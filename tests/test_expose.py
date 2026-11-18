"""McsCmd.expose end to end, with the camera, the database and the hub replaced by fakes: the PFSC file is
written whatever happens to the analysis of its frame, with the spots of that frame only."""
import logging
import pathlib
import tempfile
import time
import types
import unittest
import unittest.mock

import astropy.io.fits as pyfits
import numpy as np

import mcsActor.Commands.McsCmd as McsCmd
import mcsActor.mcsRoutines.saturation as saturation
import test_bookkeeping
from test_saturation import CONFIG
from test_spotFinder import cameraFrame, gridPositions
from test_spotMeasurement import CENT_PARMS


class Keywords(dict):
    """Keywords of a command, also indexable by position."""

    def __init__(self, keywords):
        super().__init__({name: types.SimpleNamespace(name=name, values=[value]) for name, value in keywords.items()})
        self.order = list(keywords)

    def __getitem__(self, key):
        return dict.__getitem__(self, self.order[key] if isinstance(key, int) else key)


class Cmd(test_bookkeeping.FakeCmd):
    """A command that knows whether it is still alive."""

    def __init__(self, **keywords):
        super().__init__()
        self.cmd = types.SimpleNamespace(keywords=Keywords(keywords))
        self.alive = True

    def fail(self, text):
        super().fail(text)
        self.alive = False

    def finish(self, text=''):
        super().finish(text)
        self.alive = False

    def isAlive(self):
        return self.alive


class ExposeTestCase(unittest.TestCase):

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = pathlib.Path(directory.name)
        self.bcast = Cmd()
        self.inserted = []
        config = dict(centroidParams=dict(CENT_PARMS, centSigma=25, threshSigma=4, threshFact=4),
                      detector=dict(name='MCS cam#1', gain=2.24), saturation=CONFIG, cMethod='spotFinder')

        owner = McsCmd.McsCmd.__new__(McsCmd.McsCmd)
        owner.actor = types.SimpleNamespace(actorConfig=config, models={}, bcast=self.bcast, cameraName='canon_50m')
        for name, value in dict(expTime=4800.0, prevExpTime=4800.0, newTable=None, centOverrides={},
                                dbOverride=None, findThresh=None, centThresh=None, spotFlags=None,
                                spotSaturated=None, detectorProblem='').items():
            setattr(owner, name, value)
        owner.centParms = McsCmd.centroidParameters(config['centroidParams'], owner.centOverrides)
        owner.detector, _ = McsCmd.configuredDetector(config)
        owner.cMethod = 'spotFinder'
        owner.softwareVersions = dict(version_actor='1.10.11', version_instdata='1.8.120', version_pfs_utils='7.4.47')
        owner.saturated = saturation.Saturation(None, 'not measured yet')
        owner.logger = logging.getLogger('test')
        owner.connectToDB = lambda cmd=None: types.SimpleNamespace(
            insert_dataframe=lambda table, df: self.inserted.append((table, df)))
        owner.getGeometry = lambda cmd: None
        owner.handleTelescopeGeometry = lambda *args: None
        owner._doExpose = self.doExpose
        self.owner = owner
        patcher = unittest.mock.patch.object(McsCmd.dbTools, 'loadTelescopeParametersFromDB',
                                             side_effect=lambda db, frameId: (30.0, 0.0))
        self.loadTelescope = patcher.start()
        self.addCleanup(patcher.stop)

    def doExpose(self, cmd, expTime, expType, frameId):
        fileIds = dict(pfsDay='2026-11-18', visit=self.frameId // 100, frame=self.frameId % 100,
                       frameId=self.frameId, filename=self.directory / f'PFSC{self.frameId:08d}.fits')
        return fileIds, pyfits.Header([('OBS-MOD', 'SPEC_IMAG', '')]), self.image

    def expose(self, frameId, image, **keywords):
        """Run expose on image: the command, what it raised, and the HDUs of the file written, or None."""
        self.frameId, self.image = frameId, image
        cmd = Cmd(object=None, expTime=4.8, **keywords)
        try:
            McsCmd.McsCmd.expose(self.owner, cmd)
            raised = None
        except Exception as e:
            raised = e
        path = self.directory / f'PFSC{frameId:08d}.fits'
        for _ in range(200):
            if any(text.startswith('mcsFileIds') and f'frameId={frameId};' in text for _, text in self.bcast.replies):
                break
            time.sleep(0.05)
        if not path.exists():
            return cmd, raised, None
        with pyfits.open(path) as hdus:
            hdus.verify('exception')
            return cmd, raised, {hdu.name: (hdu.header.copy(), None if hdu.data is None else hdu.data.copy())
                                 for hdu in hdus}

    def testFrameWithSpots(self):
        """The file holds the spots and the mask; mcs_data gets them; the command finishes."""
        positions = gridPositions(seed=1)
        cmd, raised, hdus = self.expose(15272400, cameraFrame(positions, seed=1))

        self.assertIsNone(raised)
        self.assertEqual(cmd.failed(), [])
        self.assertFalse(cmd.alive)
        self.assertEqual(len(hdus['SPOTS'][1]), len(positions))
        self.assertIn('MASK', hdus)
        self.assertEqual(len(self.inserted[0][1]), len(positions) + 1)

    def testNoCentroidFrameHasNoSpots(self):
        """A frame taken without centroids does not carry the spots of the frame before."""
        self.expose(15272400, cameraFrame(gridPositions(seed=1), seed=1))
        self.assertGreater(len(self.owner.centroids), 0)
        cmd, raised, hdus = self.expose(15272401, cameraFrame(gridPositions(seed=2), seed=2), noCentroid=None)

        self.assertIsNone(raised)
        self.assertEqual(len(hdus['SPOTS'][1]), 0)

    def testFrameWithoutSpots(self):
        """No spot fails the command; the file is written, without spots."""
        cmd, raised, hdus = self.expose(15272402, cameraFrame(np.zeros((0, 2)), seed=3))

        self.assertEqual(len(cmd.failed()), 1)
        self.assertIsNotNone(hdus)
        self.assertEqual(len(hdus['SPOTS'][1]), 0)

    def testDatabaseFailure(self):
        """A failure of the database after readout still writes the file, and is raised."""
        self.loadTelescope.side_effect = RuntimeError('opdb down')
        cmd, raised, hdus = self.expose(15272403, cameraFrame(gridPositions(seed=4), seed=4))

        self.assertIsInstance(raised, RuntimeError)
        self.assertIsNotNone(hdus)

    def testMissingExtension(self):
        """Without the compiled extension the command fails and the file is written."""
        with unittest.mock.patch.object(McsCmd, 'spotMeasurementModule', side_effect=ImportError('not built')):
            cmd, raised, hdus = self.expose(15272404, cameraFrame(gridPositions(seed=5), seed=5))

        self.assertIsNone(raised)
        self.assertEqual(len(cmd.failed()), 1)
        self.assertIn('not built', cmd.failed()[0])
        self.assertIsNotNone(hdus)

    def testMaskThatCannotBeWritten(self):
        """A MASK HDU that fails when written is dropped, the rest of the file written, and a warning sent."""
        class FailingHdu(pyfits.ImageHDU):
            def _writeto(self, *args, **kwargs):
                raise OSError('mask write failed')

        with unittest.mock.patch.object(saturation, 'maskHdu',
                                        side_effect=lambda shape, saturated: FailingHdu(name='MASK')):
            cmd, raised, hdus = self.expose(15272405, cameraFrame(gridPositions(seed=6), seed=6))

        self.assertIsNotNone(hdus)
        self.assertNotIn('MASK', hdus)
        self.assertIn('SPOTS', hdus)
        self.assertTrue(any(level == 'warn' and 'mask write failed' in text for level, text in self.bcast.replies))


if __name__ == '__main__':
    unittest.main()
