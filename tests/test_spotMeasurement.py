"""The spots of a frame as the actor keeps them."""
import unittest

import numpy as np

import mcsActor.mcsRoutines.spotMeasurement as spotMeasurement
import mcsActor.spotFinder.spotFinder as spotFinder
from test_spotFinder import PARAMS, cameraFrame, gridPositions

CENT_PARMS = {name: value for name, value in PARAMS.items() if name != 'gain'}
GAIN = PARAMS['gain']


class MeasureSpotsTestCase(unittest.TestCase):

    def testColumns(self):
        """One row per spot, in the order of SPOT_COLUMNS, numbered from 1."""
        image = cameraFrame(gridPositions(seed=1), seed=1)
        points, levels = spotMeasurement.measureSpots(image, CENT_PARMS, GAIN)
        spots, expected = spotFinder.findSpots(image, **PARAMS)
        column = {name: points[:, k] for k, name in enumerate(spotMeasurement.SPOT_COLUMNS)}

        self.assertEqual(points.shape, (len(spots), len(spotMeasurement.SPOT_COLUMNS)))
        self.assertEqual(levels, expected)
        np.testing.assert_array_equal(column['spotId'], np.arange(len(spots)) + 1)
        for name in ('x', 'y', 'x2', 'y2', 'xy', 'background', 'peak', 'flux', 'fluxerr'):
            np.testing.assert_array_equal(column[name], spots[name], err_msg=name)

    def testColumnsAreDistinct(self):
        """No two columns of a frame of spots hold the same values."""
        image = cameraFrame(gridPositions(seed=2), seed=2)
        points, _ = spotMeasurement.measureSpots(image, CENT_PARMS, GAIN)

        for a in range(points.shape[1]):
            for b in range(a + 1, points.shape[1]):
                self.assertFalse(np.array_equal(points[:, a], points[:, b]), msg=(a, b))

    def testParametersAreUsed(self):
        """Each centroid parameter reaches the library."""
        image = cameraFrame(gridPositions(seed=3), seed=3)
        reference, levels = spotMeasurement.measureSpots(image, CENT_PARMS, GAIN)

        points, changed = spotMeasurement.measureSpots(image, dict(CENT_PARMS, findSigma=40), GAIN)
        self.assertEqual(changed['threshold'], np.floor(40 * levels['noise']))
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, aperture=6), GAIN)
        self.assertTrue(np.all(points[:, 8] < reference[:, 8]))
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, nmin=2000), GAIN)
        self.assertEqual(len(points), 0)
        for name, value in (('boxCent', 5), ('fwhmx', 9), ('fwhmy', 4), ('maxIt', 1)):
            points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, **{name: value}), GAIN)
            self.assertFalse(np.array_equal(points[:, 1:3], reference[:, 1:3]), msg=name)
        points, _ = spotMeasurement.measureSpots(image, dict(CENT_PARMS, boxFind=30), GAIN)
        self.assertLess(len(points), len(reference))

    def testGainSetsTheFluxErrors(self):
        """The gain given sets the photon noise of the flux errors, and nothing else."""
        image = cameraFrame(gridPositions(seed=5), seed=5)
        reference, _ = spotMeasurement.measureSpots(image, CENT_PARMS, GAIN)
        points, _ = spotMeasurement.measureSpots(image, CENT_PARMS, 4 * GAIN)

        np.testing.assert_array_equal(points[:, :9], reference[:, :9])
        self.assertTrue(np.all(points[:, 9] < reference[:, 9]))

    def testNoSpot(self):
        """A frame without a spot gives no row."""
        points, _ = spotMeasurement.measureSpots(cameraFrame(np.zeros((0, 2)), seed=4), CENT_PARMS, GAIN)

        self.assertEqual(points.shape, (0, 10))


if __name__ == '__main__':
    unittest.main()
