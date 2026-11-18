"""Saturation thresholds, saturated pixels, the flags of the spots and the mask written with an image."""
import os
import tempfile
import unittest

import astropy.io.fits as pyfits
import fitsio
import numpy as np

import mcsActor.mcsRoutines.saturation as saturation
from pfs.utils.datamodel.mcs import MCSSpotFlags

CLIP = 2815.0
EXPTIME = 4.8
BIAS = 7 * saturation.BIAS_PER_READ
NROWS, NCOLS = 60, 80
UNKNOWN = int(MCSSpotFlags.NO_SAT_CHECK)
PHOT, CENT, SAT = int(MCSSpotFlags.BAD_PHOTOMETRY), int(MCSSpotFlags.BAD_CENTROID), int(MCSSpotFlags.SATURATED)
SPOT_SETTINGS = dict(spotRadius=10, badPhotometryPixels=2, badCentroidPixels=5, saturatedPixels=10,
                     fluxErrorPerPixel=0.008)
CONFIG = dict(clipLevel=CLIP, thresholdFactor=1.0, maxFiducialFraction=0.5, **SPOT_SETTINGS)
THRESHOLDS = saturation.SpotThresholds(radius=10.0, badPhotometry=2, badCentroid=5, saturated=10,
                                       fluxErrorPerPixel=0.008)


def blank(shape=(NROWS, NCOLS), dtype='u2', exptime=EXPTIME, offset=0):
    """An image at the bias level of an exposure time, plus offset."""
    return np.full(shape, round(saturation.nReads(exptime) * saturation.BIAS_PER_READ) + offset, dtype=dtype)


def thresholds(factor=1.0, nRows=NROWS, exptime=EXPTIME):
    """Integer threshold of each row."""
    return np.ceil(saturation.rowThresholds(nRows, exptime, CLIP, factor)).astype(int)


def canonWithSpot(value=None):
    """A Canon coadd with a 3 x 3 spot at row 2000, column 4000, at its row's threshold unless given."""
    image = blank(shape=saturation.CANON_SHAPE)
    image[2000:2003, 4000:4003] = thresholds(nRows=saturation.CANON_SHAPE[0])[2000] if value is None else value
    return image


class ContractTestCase(unittest.TestCase):

    def testBitsWrittenToTheDatabaseAndTheFile(self):
        """Values readers of mcs_data.flags and of the MASK HDU rely on."""
        self.assertEqual((UNKNOWN, PHOT, CENT, SAT), (1, 2, 4, 8))
        self.assertEqual(saturation.MASK_PLANES, dict(SAT=1))
        self.assertEqual(saturation.FLAG_VERSION, 1)

    def testSettingsOfTheDetection(self):
        """Values the detection was validated with."""
        self.assertAlmostEqual(saturation.NEIGHBOUR_FRACTION, 8 / 9)
        self.assertEqual(saturation.GROW_ITERATIONS, 4)
        self.assertEqual(saturation.ISOLATION_FRACTION, 0.3)


class ReadModelTestCase(unittest.TestCase):

    def testNReads(self):
        """canonexp takes one read more than the exposure time holds loop times of 0.8 s."""
        times = (0.8, 0.81, 1.0, 1.6, 0.8 * 3, 4.8, 12.0)
        self.assertEqual([saturation.nReads(t) for t in times], [2, 3, 3, 3, 4, 7, 16])

    def testReadExposuresSumToTheExposure(self):
        """A row whose reads cover the whole exposure captures all of it."""
        exposures = saturation.readExposures(np.array([0.0, 0.4, 0.8]), EXPTIME)
        self.assertEqual(exposures.shape, (3, 7))
        np.testing.assert_allclose(exposures.sum(axis=1), EXPTIME)
        self.assertAlmostEqual(exposures.max(), saturation.READ_PERIOD)

    def testLateRowLosesTheStartOfTheExposure(self):
        """A row first read later than one period after the opening misses the light before."""
        late = saturation.READ_PERIOD + 0.3
        self.assertAlmostEqual(saturation.readExposures(np.array(late), EXPTIME).sum(), EXPTIME - 0.3)
        self.assertAlmostEqual(float(saturation.splitRatio(np.array(late), EXPTIME)),
                               (EXPTIME - 0.3) / saturation.READ_PERIOD)

    def testSplitRatioOfAShortExposure(self):
        """A 0.8 s exposure falls in one read or is shared between two."""
        ratio = saturation.splitRatio(np.array([0.0, 0.4, 0.8]), 0.8)
        np.testing.assert_allclose(ratio, [1.0, 2.0, 1.0])


class RowThresholdsTestCase(unittest.TestCase):

    def testValues(self):
        """Thresholds at the bottom, middle and top rows of a Canon coadd."""
        rows = [0, 3000, 5777]
        for exptime, expected in ((0.8, [3328, 3328, 3328]), (1.6, [5718, 4549, 3584]), (2.4, [8789, 7532, 6138]),
                                  (4.8, [17983, 16481, 15087]), (12.0, [44829, 43327, 41933])):
            values = saturation.rowThresholds(5778, exptime, CLIP, 1.0)[rows]
            np.testing.assert_allclose(values, expected, atol=1.0, err_msg=f'{exptime} s')

    def testStraightLineForLongExposures(self):
        """Past two read periods the threshold falls linearly with the row."""
        rows = np.arange(5778)
        for exptime in (4.8, 12.0):
            captured = exptime - (saturation.PHASE_MAX - saturation.READ_PERIOD) - saturation.ROW_TIME_MAX * rows
            line = saturation.nReads(exptime) * saturation.BIAS_PER_READ + CLIP * captured / saturation.READ_PERIOD
            np.testing.assert_allclose(saturation.rowThresholds(5778, exptime, CLIP, 1.0), line, atol=5.0)

    def testFactorScalesTheLevelAboveTheBias(self):
        full = saturation.rowThresholds(5778, EXPTIME, CLIP, 1.0)
        lower = saturation.rowThresholds(5778, EXPTIME, CLIP, 0.9)
        np.testing.assert_allclose(lower - BIAS, 0.9 * (full - BIAS))

    def testCachedRatioIsReadOnly(self):
        """The cached split ratios cannot be altered through a returned array."""
        with self.assertRaises(ValueError):
            saturation.rowSplitRatio(5778, EXPTIME)[0] = 0.0


class UnmeasurableTestCase(unittest.TestCase):

    def reason(self, image, exptime=EXPTIME, clipLevel=CLIP, thresholdFactor=1.0):
        return saturation.unmeasurable(image, exptime, clipLevel, thresholdFactor)

    def testCanonCoadd(self):
        self.assertEqual(self.reason(blank(shape=saturation.CANON_SHAPE)), '')

    def testOtherShape(self):
        self.assertIn('not a Canon coadd', self.reason(blank(shape=(4096, 4096))))

    def testExposureTimeAgainstTheBackground(self):
        """An exposure time in the wrong unit, or one read off, disagrees with the bias level of the image."""
        image = blank(shape=saturation.CANON_SHAPE)
        for exptime in (4800.0, 0.8, 4.0, 5.6):
            self.assertIn('reads', self.reason(image, exptime=exptime), msg=f'{exptime} s')

    def testBackgroundTolerance(self):
        """The background may differ from the bias of the reads by a fraction of a read, not more."""
        self.assertEqual(self.reason(blank(shape=saturation.CANON_SHAPE, offset=40)), '')
        self.assertEqual(self.reason(blank(shape=saturation.CANON_SHAPE, offset=-40)), '')
        self.assertIn('reads', self.reason(blank(shape=saturation.CANON_SHAPE, offset=60)))
        self.assertIn('reads', self.reason(blank(shape=saturation.CANON_SHAPE, offset=-60)))

    def testSumThatCanWrap(self):
        """The longest exposure whose reads cannot sum past 16 bits is 13.6 s at this clip level."""
        self.assertEqual(self.reason(blank(shape=saturation.CANON_SHAPE, exptime=13.6), exptime=13.6), '')
        self.assertIn('wrap', self.reason(blank(shape=saturation.CANON_SHAPE, exptime=14.4), exptime=14.4))

    def testSettings(self):
        """Settings that are not levels a read can reach."""
        image = blank(shape=saturation.CANON_SHAPE)
        for clipLevel, factor in ((0.0, 1.0), (-2815.0, 1.0), (CLIP, 0.0), (CLIP, -1.0), (np.nan, 1.0)):
            self.assertIn('must be positive', self.reason(image, clipLevel=clipLevel, thresholdFactor=factor))
        self.assertIn('cannot reach', self.reason(image, clipLevel=28150.0))
        self.assertIn('cannot reach', self.reason(image, thresholdFactor=1.5))

    def testExposureTime(self):
        image = blank(shape=saturation.CANON_SHAPE)
        for exptime in (0.0, -4.8, np.nan):
            self.assertIn('exposure time', self.reason(image, exptime=exptime))

    def testThresholdBeyondTheImageType(self):
        """A 12 s threshold cannot be reached by a signed 16-bit pixel."""
        image = blank(shape=saturation.CANON_SHAPE, dtype='i2', exptime=12.0)
        self.assertIn('range of the image type', self.reason(image, exptime=12.0))


class SaturatedPixelsTestCase(unittest.TestCase):

    def setUp(self):
        self.level = thresholds()
        self.image = blank()
        threshold = saturation.rowThresholds(NROWS, EXPTIME, CLIP, 1.0)[20]
        self.neighbour = int(np.ceil(BIAS + saturation.NEIGHBOUR_FRACTION * (threshold - BIAS)))

    def mask(self, image=None):
        return saturation.saturatedPixels(self.image if image is None else image, EXPTIME, CLIP, 1.0)

    def testNothingSaturated(self):
        self.assertFalse(self.mask().any())

    def testAtAndBelowTheThreshold(self):
        """A pixel is saturated from the threshold of its row, rounded up, not one ADU below."""
        self.image[20, 30:32] = self.level[20]
        self.image[40, 30:32] = self.level[40] - 1
        mask = self.mask()
        self.assertTrue(mask[20, 30:32].all())
        self.assertEqual(int(mask.sum()), 2)

    def testIsolatedPixelIsDropped(self):
        """A pixel above the threshold with no bright neighbour is not saturated, at the border either."""
        for row, col in ((20, 30), (0, 30), (NROWS - 1, 30), (20, 0), (20, NCOLS - 1), (0, 0)):
            image = blank()
            image[row, col] = 60000
            self.assertFalse(self.mask(image).any(), msg=f'pixel {row}, {col}')

    def testIsolationLevel(self):
        """A pixel is kept when a neighbour reaches ISOLATION_FRACTION of its value above the bias."""
        net = self.level[20] - BIAS
        for fraction, kept in ((saturation.ISOLATION_FRACTION + 0.02, True),
                               (saturation.ISOLATION_FRACTION - 0.02, False)):
            image = blank()
            image[20, 30] = self.level[20]
            image[20, 31] = round(BIAS + fraction * net)
            self.assertEqual(bool(self.mask(image)[20, 30]), kept, msg=f'neighbour at {fraction:.2f}')

    def testNeighbourLevel(self):
        """A pixel touching a saturated one is saturated from NEIGHBOUR_FRACTION of the threshold above the bias."""
        self.image[20, 30:32] = self.level[20]
        self.image[20, 32] = self.neighbour
        self.image[20, 29] = self.neighbour - 1
        mask = self.mask()
        self.assertTrue(mask[20, 32])
        self.assertFalse(mask[20, 29])

    def testNeighbourMustTouch(self):
        self.image[20, 30:32] = self.level[20]
        self.image[20, 40] = self.neighbour
        self.assertFalse(self.mask()[20, 40])

    def testGrowthStops(self):
        """Pixels at the neighbour level are added GROW_ITERATIONS steps away from a saturated pixel, no further."""
        self.image[20, 30:32] = self.level[20]
        self.image[20, 32:42] = self.neighbour
        mask = self.mask()
        last = 31 + saturation.GROW_ITERATIONS
        self.assertTrue(mask[20, 30:last + 1].all())
        self.assertFalse(mask[20, last + 1:].any())

    def testThresholdOfEachRow(self):
        """The same value saturates an upper row and not a lower one."""
        coadd = blank(shape=(5778, 400))
        level = thresholds(nRows=5778)
        value = (level[100] + level[5600]) // 2
        coadd[100, 50:52] = value
        coadd[5600, 50:52] = value
        mask = saturation.saturatedPixels(coadd, EXPTIME, CLIP, 1.0)
        self.assertFalse(mask[100].any())
        self.assertTrue(mask[5600, 50:52].all())

    def testFloodedImage(self):
        """An image with more pixels above the threshold than spots can hold is refused."""
        image = blank(shape=(600, 600))
        image[:500, :500] = 60000
        with self.assertRaises(saturation.FloodedImage):
            saturation.saturatedPixels(image, EXPTIME, CLIP, 1.0)

    def testFloatImage(self):
        image = blank(dtype='f8')
        image[20, 30:32] = self.level[20]
        self.assertEqual(int(self.mask(image).sum()), 2)

    def testImageIsNotModified(self):
        self.image[20, 30:32] = self.level[20]
        before = self.image.copy()
        self.mask()
        np.testing.assert_array_equal(self.image, before)

    def testSpotAtTheImageEdge(self):
        self.image[0, 0:2] = self.level[0]
        self.image[NROWS - 1, NCOLS - 2:] = self.level[NROWS - 1]
        self.assertEqual(int(self.mask().sum()), 4)


class MeasureTestCase(unittest.TestCase):

    def testMask(self):
        saturated = saturation.measure(canonWithSpot(), EXPTIME, CLIP, 1.0)
        self.assertEqual(saturated.reason, '')
        self.assertEqual(int(saturated.mask.sum()), 9)
        self.assertTrue(saturated.mask[2000:2003, 4000:4003].all())

    def testFactorLeavesAMargin(self):
        """A spot under the clip level is saturated once the factor is lowered."""
        image = canonWithSpot(value=thresholds(nRows=saturation.CANON_SHAPE[0])[2000] - 500)
        self.assertFalse(saturation.measure(image, EXPTIME, CLIP, 1.0).mask.any())
        self.assertEqual(int(saturation.measure(image, EXPTIME, CLIP, 0.9).mask.sum()), 9)

    def testUnmeasurableImage(self):
        image = blank(shape=(4096, 4096))
        image[2000:2003, 2000:2003] = 65535
        saturated = saturation.measure(image, EXPTIME, CLIP, 1.0)
        self.assertIsNone(saturated.mask)
        self.assertIn('not a Canon coadd', saturated.reason)

    def testFloodedImage(self):
        """A coadd lit over part of the chip keeps the background of its reads and is refused."""
        image = blank(shape=saturation.CANON_SHAPE)
        image[:1000, :1000] = 60000
        saturated = saturation.measure(image, EXPTIME, CLIP, 1.0)
        self.assertIsNone(saturated.mask)
        self.assertIn('above their threshold', saturated.reason)

    def testConfigured(self):
        """The mask, with the settings it was measured with and the thresholds of the spots."""
        saturated = saturation.measureConfigured(CONFIG, canonWithSpot(), EXPTIME)
        self.assertEqual(saturated.reason, '')
        self.assertEqual(int(saturated.mask.sum()), 9)
        self.assertEqual((saturated.clipLevel, saturated.thresholdFactor, saturated.spots), (CLIP, 1.0, THRESHOLDS))

    def testNoConfiguration(self):
        saturated = saturation.measureConfigured(None, canonWithSpot(), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertIn('no saturation configuration', saturated.reason)

    def testMissingSetting(self):
        """A configuration without one of its settings gives no mask, no thresholds, and names the setting."""
        for name in CONFIG:
            if name == 'maxFiducialFraction':
                continue
            config = {key: value for key, value in CONFIG.items() if key != name}
            saturated = saturation.measureConfigured(config, canonWithSpot(), EXPTIME)
            self.assertIsNone(saturated.mask, msg=name)
            self.assertIsNone(saturated.spots, msg=name)
            self.assertIn(name, saturated.reason)

    def testConfigurationOfAnotherType(self):
        for config in ('2815', [CLIP, 1.0], dict(CONFIG, clipLevel='high')):
            saturated = saturation.measureConfigured(config, canonWithSpot(), EXPTIME)
            self.assertIsNone(saturated.mask)
            self.assertTrue(saturated.reason)

    def testSpotThresholdsRefused(self):
        """Thresholds that cannot be applied leave the saturation unknown."""
        for change in (dict(spotRadius=0), dict(spotRadius=float('nan')), dict(badPhotometryPixels=0),
                       dict(badCentroidPixels=1.5), dict(badCentroidPixels=True), dict(saturatedPixels=4),
                       dict(badPhotometryPixels=6), dict(fluxErrorPerPixel=-0.001)):
            saturated = saturation.measureConfigured(dict(CONFIG, **change), canonWithSpot(), EXPTIME)
            self.assertIsNone(saturated.mask, msg=change)
            self.assertIsNone(saturated.spots, msg=change)
            self.assertIn('saturation configuration', saturated.reason)

    def testFailureKeepsTheSettings(self):
        """An image that cannot be measured still records the settings it was to be measured with."""
        saturated = saturation.measureConfigured(CONFIG, blank(shape=(4096, 4096)), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertEqual((saturated.clipLevel, saturated.spots), (CLIP, THRESHOLDS))

    def testSettingOutOfRange(self):
        """A clip level no read can reach leaves the saturation unknown, it does not report none."""
        config = dict(CONFIG, clipLevel=28150.0)
        saturated = saturation.measureConfigured(config, canonWithSpot(value=60000), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertIn('cannot reach', saturated.reason)


class PixelsNearSpotsTestCase(unittest.TestCase):

    def setUp(self):
        self.mask = np.zeros((NROWS, NCOLS), dtype=bool)
        self.mask[20, 30:33] = True    # three pixels in a row
        self.mask[50, 70] = True

    def testCounts(self):
        """The mask pixels within the radius, the centroid's own pixel included."""
        x, y = np.array([31.0, 31.0, 31.0, 40.0, 70.0]), np.array([20.0, 22.0, 26.0, 20.0, 50.0])
        counts = saturation.pixelsNearSpots(self.mask, x, y, radius=5.0)
        self.assertEqual(counts.tolist(), [3, 3, 0, 0, 1])
        self.assertEqual(counts.dtype, np.int32)

    def testRadiusIsInclusive(self):
        """A pixel exactly at the radius counts; one beyond does not."""
        self.assertEqual(saturation.pixelsNearSpots(self.mask, np.array([30.0]), np.array([25.0]), 5.0).tolist(), [1])
        self.assertEqual(saturation.pixelsNearSpots(self.mask, np.array([30.0]), np.array([25.5]), 5.0).tolist(), [0])

    def testSubPixelCentroid(self):
        """Distances are taken from the centroid, not from its pixel."""
        mask = np.zeros((NROWS, NCOLS), dtype=bool)
        mask[20, 36] = True
        self.assertEqual(saturation.pixelsNearSpots(mask, np.array([31.4]), np.array([20.0]), 5.0).tolist(), [1])
        self.assertEqual(saturation.pixelsNearSpots(mask, np.array([30.6]), np.array([20.0]), 5.0).tolist(), [0])

    def testNonFiniteAndOffImage(self):
        """A spot at a non-finite position, or off the image, has no saturated pixel; the edge does not wrap."""
        mask = np.zeros((NROWS, NCOLS), dtype=bool)
        mask[1, 0] = True
        mask[1, NCOLS - 1] = True
        counts = saturation.pixelsNearSpots(mask, np.array([np.nan, 1.0, -2.0, -20.0, 500.0]),
                                            np.array([1.0, np.nan, 1.0, 1.0, 1.0]), 5.0)
        self.assertEqual(counts.tolist(), [0, 0, 1, 0, 0])

    def testWithoutAMask(self):
        self.assertEqual(saturation.pixelsNearSpots(None, np.array([1.0, 2.0]), np.array([1.0, 2.0]), 5.0).tolist(),
                         [-1, -1])

    def testNoSpot(self):
        self.assertEqual(saturation.pixelsNearSpots(self.mask, np.array([]), np.array([]), 5.0).shape, (0,))


class SpotFlagsTestCase(unittest.TestCase):

    def testThresholds(self):
        """Each flag from its threshold on, inclusive, the flags nesting."""
        counts = np.array([0, 1, 2, 4, 5, 9, 10, 60])
        flags = saturation.spotFlags(counts, THRESHOLDS)
        self.assertEqual(flags.tolist(), [0, 0, PHOT, PHOT, PHOT | CENT, PHOT | CENT, PHOT | CENT | SAT,
                                          PHOT | CENT | SAT])
        self.assertEqual(flags.dtype, np.int32)

    def testUnknownCount(self):
        """A spot whose count is unknown carries NO_SAT_CHECK alone."""
        self.assertEqual(saturation.spotFlags(np.array([-1, 12]), THRESHOLDS).tolist(), [UNKNOWN, PHOT | CENT | SAT])

    def testWithoutThresholds(self):
        self.assertEqual(saturation.spotFlags(np.array([0, 12]), None).tolist(), [UNKNOWN, UNKNOWN])

    def testEqualThresholds(self):
        """Thresholds may coincide."""
        thresholds = THRESHOLDS._replace(badPhotometry=3, badCentroid=3, saturated=3)
        self.assertEqual(saturation.spotFlags(np.array([2, 3]), thresholds).tolist(), [0, PHOT | CENT | SAT])

    def testNoSpot(self):
        self.assertEqual(saturation.spotFlags(np.array([], dtype=int), THRESHOLDS).shape, (0,))


class FluxErrorsTestCase(unittest.TestCase):

    def testTermAddedInQuadrature(self):
        flux, fluxerr = np.array([100000.0, 100000.0, 50000.0]), np.array([300.0, 300.0, 200.0])
        errors = saturation.fluxErrors(flux, fluxerr, np.array([0, 2, 5]), THRESHOLDS)
        np.testing.assert_allclose(errors, [300.0, np.hypot(300.0, 1600.0), np.hypot(200.0, 2000.0)])

    def testUnknownOrNotFinite(self):
        """An unknown count, a missing threshold or a flux that is not finite leaves the error as it was."""
        flux, fluxerr = np.array([100000.0, np.nan, 100000.0]), np.array([300.0, np.nan, 300.0])
        errors = saturation.fluxErrors(flux, fluxerr, np.array([-1, 4, 3]), THRESHOLDS)
        np.testing.assert_allclose(errors, [300.0, np.nan, np.hypot(300.0, 2400.0)])
        np.testing.assert_array_equal(saturation.fluxErrors(flux, fluxerr, np.array([3, 3, 3]), None), fluxerr)

    def testInputsUntouched(self):
        flux, fluxerr = np.array([100000.0]), np.array([300.0])
        saturation.fluxErrors(flux, fluxerr, np.array([3]), THRESHOLDS)
        self.assertEqual(fluxerr.tolist(), [300.0])


class FlagsOfSpotsTestCase(unittest.TestCase):

    def setUp(self):
        self.spotIds = np.arange(1, 5)
        self.flags = np.array([0, PHOT, PHOT | CENT, PHOT | CENT | SAT])

    def testWanted(self):
        """The flags of the identifiers wanted, in their order; those that are not positive are left out."""
        flags = saturation.flagsOfSpots(self.flags, self.spotIds, np.array([4.0, -1.0, 1.0, 2.0, 0.0]))
        self.assertEqual(flags.tolist(), [PHOT | CENT | SAT, 0, PHOT])

    def testSpotIdsInAnyOrder(self):
        flags = saturation.flagsOfSpots(self.flags, np.array([40, 10, 30, 20]), np.array([10, 20, 30, 40]))
        self.assertEqual(flags.tolist(), [PHOT, PHOT | CENT | SAT, CENT | PHOT, 0])

    def testIdentifierThatIsNotASpot(self):
        """An identifier that is not a spot is unknown, not unsaturated."""
        flags = saturation.flagsOfSpots(self.flags, np.array([2, 4, 6, 8]), np.array([1, 5, 9, 4]))
        self.assertEqual(flags.tolist(), [UNKNOWN, UNKNOWN, UNKNOWN, PHOT])

    def testNoSpot(self):
        flags = saturation.flagsOfSpots(np.array([], dtype=int), np.array([]), np.array([3, -1]))
        self.assertEqual(flags.tolist(), [UNKNOWN])
        self.assertEqual(saturation.flagsOfSpots(self.flags, self.spotIds, np.array([-1.0, -1.0])).shape, (0,))


class ReportsTestCase(unittest.TestCase):

    FLAGS = np.array([0, PHOT, PHOT | CENT, PHOT | CENT | SAT, PHOT | CENT | SAT])

    def testCountFlags(self):
        self.assertEqual(saturation.countFlags(self.FLAGS), (4, 3, 2))
        self.assertEqual(saturation.countFlags(np.array([], dtype=int)), (0, 0, 0))
        self.assertIsNone(saturation.countFlags(np.array([0, UNKNOWN])))

    def testCountKeyword(self):
        self.assertEqual(saturation.countKeyword('saturatedSpots', 15272800, self.FLAGS),
                         'saturatedSpots=15272800,4,3,2,5')

    def testCountKeywordWhenUnknown(self):
        """Unknown counts are -1, the number of spots is still given."""
        self.assertEqual(saturation.countKeyword('saturatedCobras', 15272800, np.full(3, UNKNOWN)),
                         'saturatedCobras=15272800,-1,-1,-1,3')

    def testExceedsFraction(self):
        """The fraction is of SATURATED spots unless another flag is named."""
        self.assertFalse(saturation.exceedsFraction(self.FLAGS, 0.4))
        self.assertTrue(saturation.exceedsFraction(self.FLAGS, 0.39))
        self.assertTrue(saturation.exceedsFraction(self.FLAGS, '0.39'))
        self.assertTrue(saturation.exceedsFraction(self.FLAGS, 0.7, flag=MCSSpotFlags.BAD_PHOTOMETRY))
        self.assertFalse(saturation.exceedsFraction(self.FLAGS, None))
        self.assertFalse(saturation.exceedsFraction(np.array([], dtype=int), 0.0))

    def testUnknownIsNotAnExcess(self):
        self.assertFalse(saturation.exceedsFraction(np.array([SAT, UNKNOWN]), 0.0))


class SpotCardsTestCase(unittest.TestCase):

    def header(self, thresholds):
        """The cards as a header written to a file and read back."""
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'spots.fits')
            table = pyfits.BinTableHDU.from_columns([pyfits.Column(name='flags', format='J', array=[0])], name='SPOTS')
            table.header.extend(saturation.spotCards(thresholds))
            pyfits.HDUList([pyfits.PrimaryHDU(), table]).writeto(path, output_verify='exception')
            with pyfits.open(path) as hdus:
                return hdus['SPOTS'].header.copy()

    def testDefinitions(self):
        """FLAG<n> names the flag worth 2**n, each card fits without truncation."""
        header = self.header(THRESHOLDS)
        self.assertEqual(header['FLAGVER'], 1)
        for flag in MCSSpotFlags:
            if flag is MCSSpotFlags.NONE:
                continue
            self.assertEqual(MCSSpotFlags[header[f'FLAG{int(flag).bit_length() - 1}']], flag)
            self.assertTrue(header.comments[f'FLAG{int(flag).bit_length() - 1}'])
        self.assertEqual((header['SATRAD'], header['SATPHOT'], header['SATCENT'], header['SATSPOT'], header['SATFERR']),
                         (10.0, 2, 5, 10, 0.008))

    def testWithoutThresholds(self):
        """The definitions are written even when the frame had no thresholds."""
        header = self.header(None)
        self.assertEqual(header['FLAG0'], 'NO_SAT_CHECK')
        self.assertNotIn('SATRAD', header)


class MaskHduTestCase(unittest.TestCase):

    def written(self, hdu):
        """The HDU as astropy and as fitsio read it back from a file: (data, header) of each."""
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'mask.fits')
            pyfits.HDUList([pyfits.PrimaryHDU(), hdu]).writeto(path)
            with pyfits.open(path) as hdus:
                fromAstropy = hdus['MASK'].data.copy(), hdus['MASK'].header.copy()
            with fitsio.FITS(path) as hdus:
                fromFitsio = hdus['MASK'].read(), hdus['MASK'].read_header()
        return fromAstropy, fromFitsio

    def testSaturatedPixels(self):
        """A saturated pixel has the bit of the SAT plane, which the MP_SAT card names, with the settings."""
        mask = np.zeros((NROWS, NCOLS), dtype=bool)
        mask[20, 30:32] = True
        saturated = saturation.Saturation(mask, '', CLIP, 1.0, THRESHOLDS)
        for data, header in self.written(saturation.maskHdu(mask.shape, saturated)):
            self.assertEqual(data.shape, mask.shape)
            np.testing.assert_array_equal(data != 0, mask)
            self.assertEqual(int(data[20, 30]), 1 << int(header['MP_SAT']))
            self.assertTrue(header['SATKNOWN'])
            self.assertNotIn('SATWHY', header)
            self.assertEqual((header['CLIPLVL'], header['SATFACT']), (CLIP, 1.0))

    def testUnknown(self):
        """Without a measurement no bit is set and the header says so, with the reason."""
        reason = 'a 4.8 s exposure sums 7 reads but the background (1283 ADU) is the bias of 5.0'
        for data, header in self.written(saturation.maskHdu((NROWS, NCOLS), saturation.Saturation(None, reason))):
            self.assertEqual(data.shape, (NROWS, NCOLS))
            self.assertFalse(data.any())
            self.assertFalse(header['SATKNOWN'])
            self.assertEqual(header['SATWHY'], reason)
            self.assertNotIn('CLIPLVL', header)

    def testReasonACardCannotHold(self):
        """A reason with a line break or a non-ASCII character is still written."""
        saturated = saturation.Saturation(None, 'KeyError: \u00e9t\u00e9\nline two')
        for data, header in self.written(saturation.maskHdu((NROWS, NCOLS), saturated)):
            self.assertFalse(header['SATKNOWN'])
            self.assertEqual(header['SATWHY'], 'KeyError:  t  line two')


class DutyScanTestCase(unittest.TestCase):
    """The flags against the calibration of docs/saturation.md, on spots built to its numbers."""

    def testCalibratedCounts(self):
        """A spot with as many masked pixels as the calibration's classes gets their flags."""
        mask = np.zeros((NROWS, NCOLS), dtype=bool)
        positions = []
        for k, n in enumerate((1, 2, 5, 10)):
            x0, y0 = 10 + 18 * k, 30
            rows, cols = np.unravel_index(np.arange(n), (4, 4))
            mask[y0 - 1 + rows, x0 - 1 + cols] = True
            positions.append((x0, y0))
        x, y = np.array(positions, dtype=float).T
        counts = saturation.pixelsNearSpots(mask, x, y, THRESHOLDS.radius)
        self.assertEqual(counts.tolist(), [1, 2, 5, 10])
        self.assertEqual(saturation.spotFlags(counts, THRESHOLDS).tolist(),
                         [0, PHOT, PHOT | CENT, PHOT | CENT | SAT])


if __name__ == '__main__':
    unittest.main()
