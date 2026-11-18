"""What the actor does with the saturation of a frame: spot flags and flux errors, the SPOTS table and its
header, the rows of mcs_data and the reports of matched fiducials and cobras."""
import os
import tempfile
import types
import unittest

import astropy.io.fits as pyfits
import numpy as np

import mcsActor.Commands.McsCmd as McsCmd
import mcsActor.mcsRoutines.saturation as saturation
from pfs.utils.datamodel.mcs import MCSSpotFlags
from test_bookkeeping import FakeCmd, commandSet
from test_saturation import CONFIG, THRESHOLDS

UNKNOWN = int(MCSSpotFlags.NO_SAT_CHECK)
PHOT, CENT, SAT = int(MCSSpotFlags.BAD_PHOTOMETRY), int(MCSSpotFlags.BAD_CENTROID), int(MCSSpotFlags.SATURATED)


def centroids(x, y, flux=100000.0, fluxerr=300.0):
    """The actor's array of centroids (spotMeasurement.SPOT_COLUMNS) for spots at x, y."""
    n = len(x)
    points = np.zeros((n, 10))
    points[:, 0] = np.arange(n) + 1
    points[:, 1], points[:, 2] = x, y
    points[:, 8], points[:, 9] = flux, fluxerr
    return points


def maskWith(counts, shape=(60, 120)):
    """A mask and spot positions, each spot with counts[k] masked pixels next to its centroid."""
    mask = np.zeros(shape, dtype=bool)
    x = 12.0 + 24.0 * np.arange(len(counts))
    y = np.full(len(counts), 30.0)
    for k, n in enumerate(counts):
        rows, cols = np.unravel_index(np.arange(n), (4, 4))
        mask[29 + rows, int(x[k]) - 1 + cols] = True
    return mask, x, y


def owner(points, **config):
    """A command set holding points, with the saturation configuration given."""
    cmdSet = commandSet()
    cmdSet.actor.actorConfig = dict(cmdSet.actor.actorConfig, saturation=dict(CONFIG, **config))
    cmdSet.centroids = points
    cmdSet.spotFlags = cmdSet.spotSaturated = None
    cmdSet.saturated = saturation.Saturation(None, 'not measured yet')
    return cmdSet


class FlagSpotsTestCase(unittest.TestCase):

    def testFlagsAndFluxErrors(self):
        """Counts, nested flags, flux errors widened by the saturated pixels, and the keyword."""
        mask, x, y = maskWith([0, 1, 2, 5, 10])
        cmdSet, cmd = owner(centroids(x, y)), FakeCmd()
        saturated = saturation.Saturation(mask, '', 2815.0, 1.0, THRESHOLDS)
        McsCmd.McsCmd.flagSaturatedSpots(cmdSet, cmd, 15272800, saturated)

        self.assertEqual(cmdSet.spotSaturated.tolist(), [0, 1, 2, 5, 10])
        self.assertEqual(cmdSet.spotFlags.tolist(), [0, 0, PHOT, PHOT | CENT, PHOT | CENT | SAT])
        np.testing.assert_allclose(cmdSet.centroids[:, 9],
                                   np.hypot(300.0, 0.008 * np.array([0, 1, 2, 5, 10]) * 100000.0))
        self.assertIn(('inform', 'saturatedSpots=15272800,3,2,1,5'), cmd.replies)

    def testWithoutAMask(self):
        """An unmeasured frame marks every spot unknown and leaves the flux errors."""
        cmdSet, cmd = owner(centroids(np.array([10.0, 20.0]), np.array([10.0, 10.0]))), FakeCmd()
        McsCmd.McsCmd.flagSaturatedSpots(cmdSet, cmd, 15272800, saturation.Saturation(None, 'flooded', 2815.0,
                                                                                         1.0, THRESHOLDS))
        self.assertEqual(cmdSet.spotSaturated.tolist(), [-1, -1])
        self.assertEqual(cmdSet.spotFlags.tolist(), [UNKNOWN, UNKNOWN])
        self.assertEqual(cmdSet.centroids[:, 9].tolist(), [300.0, 300.0])
        self.assertIn(('inform', 'saturatedSpots=15272800,-1,-1,-1,2'), cmd.replies)

    def testFailureIsContained(self):
        """A failure marks the spots unknown and warns; the exposure goes on."""
        cmdSet, cmd = owner(centroids(np.array([10.0]), np.array([10.0]))), FakeCmd()
        broken = saturation.Saturation(np.zeros((3, 3, 3), dtype=bool), '', 2815.0, 1.0, THRESHOLDS)
        McsCmd.McsCmd.flagSaturatedSpots(cmdSet, cmd, 15272800, broken)

        self.assertEqual(cmdSet.spotFlags.tolist(), [UNKNOWN])
        self.assertEqual([level for level, _ in cmd.replies], ['warn', 'inform'])


class SpotHduTestCase(unittest.TestCase):

    def written(self, cmdSet):
        """The SPOTS HDU of the command set, written to a file and read back."""
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'spots.fits')
            pyfits.HDUList([pyfits.PrimaryHDU(), McsCmd.McsCmd._makeSpotHDU(cmdSet, FakeCmd())]).writeto(
                path, output_verify='exception')
            with pyfits.open(path) as hdus:
                return hdus['SPOTS'].data.copy(), hdus['SPOTS'].header.copy()

    def testColumnsAndHeader(self):
        mask, x, y = maskWith([0, 2, 10])
        cmdSet = owner(centroids(x, y))
        cmdSet.saturated = saturation.Saturation(mask, '', 2815.0, 1.0, THRESHOLDS)
        McsCmd.McsCmd.flagSaturatedSpots(cmdSet, FakeCmd(), 15272800, cmdSet.saturated)
        data, header = self.written(cmdSet)

        self.assertEqual(data['flags'].tolist(), [0, PHOT, PHOT | CENT | SAT])
        self.assertEqual(data['nSatPix'].tolist(), [0, 2, 10])
        self.assertEqual((header['FLAG1'], header['SATPHOT'], header['SATSPOT']), ('BAD_PHOTOMETRY', 2, 10))

    def testSpotsNotJudged(self):
        """Centroids whose saturation was not judged are written unknown, with the definitions."""
        cmdSet = owner(centroids(np.array([10.0, 20.0]), np.array([10.0, 10.0])))
        data, header = self.written(cmdSet)

        self.assertEqual(data['flags'].tolist(), [UNKNOWN, UNKNOWN])
        self.assertEqual(data['nSatPix'].tolist(), [-1, -1])
        self.assertEqual(header['FLAG0'], 'NO_SAT_CHECK')
        self.assertNotIn('SATRAD', header)

    def testStaleFlagsNotWritten(self):
        """Flags of another number of spots, left by an earlier frame, are not written."""
        cmdSet = owner(centroids(np.array([10.0, 20.0]), np.array([10.0, 10.0])))
        cmdSet.spotFlags, cmdSet.spotSaturated = np.array([PHOT]), np.array([3])
        data, _ = self.written(cmdSet)
        self.assertEqual(data['flags'].tolist(), [UNKNOWN, UNKNOWN])


class DatabaseTestCase(unittest.TestCase):

    def testFlagsWrittenToMcsData(self):
        """mcs_data gets one flag per centroid, and 0 on the spot_id -1 row."""
        inserted = []
        cmdSet = owner(centroids(np.array([10.0, 20.0, 30.0]), np.array([10.0, 10.0, 10.0])))
        cmdSet.connectToDB = lambda cmd: types.SimpleNamespace(
            insert_dataframe=lambda table, df: inserted.append((table, df)))
        McsCmd.McsCmd._writeCentroids(cmdSet, cmdSet.centroids, np.array([0, PHOT | CENT, UNKNOWN]), 15272800, 1)

        table, df = inserted[0]
        self.assertEqual(table, 'mcs_data')
        self.assertEqual(dict(zip(df['spot_id'], df['flags'])), {1: 0, 2: PHOT | CENT, 3: UNKNOWN, -1: 0})


class ReportTestCase(unittest.TestCase):

    def setUp(self):
        self.cmdSet = owner(centroids(np.arange(1, 5) * 10.0, np.full(4, 10.0)), maxFiducialFraction=0.3)
        self.cmdSet.spotFlags = np.array([0, PHOT, PHOT | CENT | SAT, PHOT | CENT | SAT])

    def testFiducials(self):
        """The matched spots' counts; a warning past the fraction of SATURATED ones."""
        cmd = FakeCmd()
        McsCmd.McsCmd.reportSaturated(self.cmdSet, cmd, 'fiducials', 15272800, np.array([1, 3, 4]),
                                      limitName='maxFiducialFraction')
        self.assertEqual(cmd.replies[0], ('inform', 'saturatedFiducials=15272800,2,2,2,3'))
        self.assertEqual(cmd.replies[1][0], 'warn')
        self.assertIn('2 of the 3 matched fiducials are SATURATED', cmd.replies[1][1])

    def testCobrasWithoutLimit(self):
        """Unmatched cobras (spot id -1) are left out, and no limit means no warning."""
        cmd = FakeCmd()
        McsCmd.McsCmd.reportSaturated(self.cmdSet, cmd, 'cobras', 15272800, np.array([-1, 2, 1, -1]))
        self.assertEqual(cmd.replies, [('inform', 'saturatedCobras=15272800,1,0,0,2')])

    def testFailureIsContained(self):
        cmd = FakeCmd()
        self.cmdSet.spotFlags = None
        self.cmdSet.centroids = None
        McsCmd.McsCmd.reportSaturated(self.cmdSet, cmd, 'cobras', 15272800, np.array([1]))
        self.assertEqual([level for level, _ in cmd.replies], ['warn'])


class MeasureTestCase(unittest.TestCase):

    def testReasonReported(self):
        """A Canon coadd that cannot be measured warns; another image only informs."""
        cmdSet = owner(None)
        cmd = FakeCmd()
        McsCmd.McsCmd.measureSaturation(cmdSet, cmd, 15272800, np.zeros(saturation.CANON_SHAPE, dtype='u2'), 4.8)
        self.assertEqual(cmd.replies[0][0], 'warn')
        cmd = FakeCmd()
        McsCmd.McsCmd.measureSaturation(cmdSet, cmd, 15272800, np.zeros((100, 100), dtype='u2'), 4.8)
        self.assertEqual(cmd.replies[0][0], 'inform')


if __name__ == '__main__':
    unittest.main()
