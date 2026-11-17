"""The spots of a frame as the actor keeps them."""
import unittest

import numpy as np

import mcsActor.mcsRoutines.spotMeasurement as spotMeasurement
import mcsActor.spotFinder.spotFinder as spotFinder
from test_spotFinder import PARAMS, cameraFrame, gridPositions

CENT_PARMS = {name: value for name, value in PARAMS.items() if name != 'gain'}


class MeasureSpotsTestCase(unittest.TestCase):

    def testColumns(self):
        """One row per spot, in the order of SPOT_COLUMNS, numbered from 1."""
        image = cameraFrame(gridPositions(seed=1), seed=1)
        points, levels = spotMeasurement.measureSpots(image, CENT_PARMS)
        spots, expected = spotFinder.findSpots(image, **dict(PARAMS, gain=spotMeasurement.GAIN))
        column = {name: points[:, k] for k, name in enumerate(spotMeasurement.SPOT_COLUMNS)}

        self.assertEqual(points.shape, (len(spots), len(spotMeasurement.SPOT_COLUMNS)))
        self.assertEqual(levels, expected)
        np.testing.assert_array_equal(column['spotId'], np.arange(len(spots)) + 1)
        for name in ('x', 'y', 'x2', 'y2', 'xy', 'background', 'peak', 'flux', 'fluxerr'):
            np.testing.assert_array_equal(column[name], spots[name], err_msg=name)

    def testColumnsAreDistinct(self):
        """No two columns of a frame of spots hold the same values."""
        image = cameraFrame(gridPositions(seed=2), seed=2)
        points, _ = spotMeasurement.measureSpots(image, CENT_PARMS)

        for a in range(points.shape[1]):
            for b in range(a + 1, points.shape[1]):
                self.assertFalse(np.array_equal(points[:, a], points[:, b]), msg=(a, b))

    def testParametersAreUsed(self):
        """Each centroid parameter reaches the library."""
        image = cameraFrame(gridPositions(seed=3), seed=3)
        reference, levels = spotMeasurement.measureSpots(image, CENT_PARMS)

        points, changed = spotMeasurement.measureSpots(image, dict(CENT_PARMS, findSigma=40))
        self.assertEqual(changed['threshold'], np.floor(40 * levels['noise']))
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, aperture=6))
        self.assertTrue(np.all(points[:, 8] < reference[:, 8]))
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, nmin=2000))
        self.assertEqual(len(points), 0)
        for name, value in (('boxCent', 5), ('fwhmx', 9), ('fwhmy', 4), ('maxIt', 1)):
            points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, **{name: value}))
            self.assertFalse(np.array_equal(points[:, 1:3], reference[:, 1:3]), msg=name)
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, boxFind=30))
        self.assertLess(len(points), len(reference))

    def testNoSpot(self):
        """A frame without a spot gives no row."""
        points, _ = spotMeasurement.measureSpots(cameraFrame(np.zeros((0, 2)), seed=4), CENT_PARMS)

        self.assertEqual(points.shape, (0, 10))


if __name__ == '__main__':
    unittest.main()
