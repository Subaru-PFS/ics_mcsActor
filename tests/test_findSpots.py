"""The findSpots step of the actor, and the switch between it and the windowed centroid."""
import types
import unittest

import numpy as np

import mcsActor.Commands.McsCmd as McsCmd
import mcsActor.mcsRoutines.spotMeasurement as spotMeasurement
from test_spotFinder import cameraFrame, gridPositions, noisyFrame
from test_spotMeasurement import CENT_PARMS


class FakeCmd:
    """Keeps the replies of a command."""

    def __init__(self, **keywords):
        self.cmd = types.SimpleNamespace(keywords={name: types.SimpleNamespace(values=[value])
                                                   for name, value in keywords.items()})
        self.replies = []

    def debug(self, text):
        self.replies.append(('debug', text))

    def inform(self, text):
        self.replies.append(('inform', text))

    def warn(self, text):
        self.replies.append(('warn', text))

    def fail(self, text):
        self.replies.append(('fail', text))

    def finish(self, text):
        self.replies.append(('finish', text))

    def failed(self):
        return [text for level, text in self.replies if level == 'fail']


def commandSet(image):
    """Stand-in for the command set, with thresholds as the windowed centroid would have left them."""
    return types.SimpleNamespace(actor=types.SimpleNamespace(image=image, cameraName='canon_50m'), findThresh=2010.0,
                                 centThresh=2010.0, avBack=1795.0, centroids=None, nCentroid=None, cMethod='win')


class FindSpotsTestCase(unittest.TestCase):

    def testSpotsKept(self):
        """The spots of the frame are kept, and the thresholds of the windowed centroid are left alone."""
        positions = gridPositions(seed=1)
        image = cameraFrame(positions, seed=1)
        owner, cmd = commandSet(image), FakeCmd()
        McsCmd.McsCmd.findSpots(owner, cmd, CENT_PARMS)

        self.assertEqual(cmd.failed(), [])
        np.testing.assert_array_equal(owner.centroids, spotMeasurement.measureSpots(image, CENT_PARMS)[0])
        self.assertEqual(owner.nCentroid, len(positions))
        self.assertEqual((owner.findThresh, owner.centThresh, owner.avBack), (2010.0, 2010.0, 1795.0))

    def testNoSpot(self):
        """A frame without a spot fails the command and keeps no centroid."""
        owner, cmd = commandSet(cameraFrame(np.zeros((0, 2)), seed=2)), FakeCmd()
        McsCmd.McsCmd.findSpots(owner, cmd, CENT_PARMS)

        self.assertEqual(len(cmd.failed()), 1)
        self.assertIsNone(owner.centroids)

    def testFloodedFrame(self):
        """A frame the library refuses fails the command, with its reason."""
        rng = np.random.default_rng(3)
        image = noisyFrame((2600, 3000), seed=3)
        image[rng.random(image.shape) < 0.45] += 3000
        owner, cmd = commandSet(image), FakeCmd()
        McsCmd.McsCmd.findSpots(owner, cmd, CENT_PARMS)

        self.assertEqual(len(cmd.failed()), 1)
        self.assertIn('too many pixels above the threshold', cmd.failed()[0])
        self.assertIsNone(owner.centroids)


class FakeDB:
    """Keeps the keywords of what is inserted."""

    def __init__(self):
        self.inserted = []

    def insert_kw(self, table, **values):
        self.inserted.append((table, values))


class MethodTestCase(unittest.TestCase):

    def testConfiguredMethod(self):
        """The config names the method, the windowed centroid when it names none; any other name is refused."""
        self.assertEqual(McsCmd.configuredCMethod({}), 'win')
        for cMethod in McsCmd.C_METHODS:
            self.assertEqual(McsCmd.configuredCMethod({'cMethod': cMethod}), cMethod)
        for cMethod in ('spotfinder', 'sep', 'cent', ''):
            with self.assertRaises(ValueError):
                McsCmd.configuredCMethod({'cMethod': cMethod})

    def testSwitchForgetsTheThresholds(self):
        """Changing the method sets it and drops the thresholds kept from the frames before."""
        owner, cmd = commandSet(None), FakeCmd(cMethod='spotFinder')
        McsCmd.McsCmd.switchCMethod(owner, cmd)

        self.assertEqual(owner.cMethod, 'spotFinder')
        self.assertEqual((owner.findThresh, owner.centThresh), (None, None))
        self.assertEqual(cmd.failed(), [])

    def testSwitchToTheSameMethod(self):
        """Switching to the method in use keeps the thresholds."""
        owner, cmd = commandSet(None), FakeCmd(cMethod='win')
        McsCmd.McsCmd.switchCMethod(owner, cmd)

        self.assertEqual(owner.cMethod, 'win')
        self.assertEqual((owner.findThresh, owner.centThresh), (2010.0, 2010.0))

    def testSwitchToAnUnknownMethod(self):
        """An unknown method fails the command and changes nothing."""
        owner, cmd = commandSet(None), FakeCmd(cMethod='spotfinder')
        McsCmd.McsCmd.switchCMethod(owner, cmd)

        self.assertEqual(len(cmd.failed()), 1)
        self.assertEqual((owner.cMethod, owner.findThresh), ('win', 2010.0))

    def testMethodRecorded(self):
        """mcs_exposure records the method that measured the frame."""
        info = dict(starttime='2026-11-17 10:00:00', frameid=15503500, visitid=155035, exptime=4800.0, altitude=90.0,
                    azimuth=0.0, instrot=0.0, adc_pa=0.0, dome_temperature=0.0, dome_pressure=0.0, dome_humidity=0.0,
                    outside_temperature=0.0, outside_pressure=0.0, outside_humidity=0.0, mcs_cover_temperature=0.0,
                    mcs_m1_temperature=0.0)
        for cMethod, recorded in (('win', "'windowed'"), ('spotFinder', "'spotFinder'")):
            db = FakeDB()
            owner = commandSet(None)
            owner.cMethod = cMethod
            owner.softwareVersions = {}
            owner.connectToDB = lambda cmd: db
            McsCmd.McsCmd._writeTelescopeInfo(owner, FakeCmd(), info)

            self.assertEqual(db.inserted[0][1]['measurement_algorithm'], recorded)


if __name__ == '__main__':
    unittest.main()
