"""Saturation thresholds, saturated pixels, the flags of the spots and the mask written with an image."""
import os
import tempfile
import unittest

import astropy.io.fits as pyfits
import fitsio
import numpy as np

import mcsActor.mcsRoutines.saturation as saturation
from mcsActor.mcsRoutines.saturation import SpotFlag

CLIP = 2815.0
EXPTIME = 4.8
BIAS = 7 * saturation.BIAS_PER_READ
NROWS, NCOLS = 60, 80
UNKNOWN = int(SpotFlag.SATURATION_UNKNOWN)


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
        self.assertEqual((int(SpotFlag.SATURATED), int(SpotFlag.SATURATION_UNKNOWN)), (1, 2))
        self.assertEqual(saturation.MASK_PLANES, dict(SAT=1))

    def testSettingsOfTheDetection(self):
        """Values the detection was validated with."""
        self.assertAlmostEqual(saturation.NEIGHBOUR_FRACTION, 8 / 9)
        self.assertEqual(saturation.GROW_ITERATIONS, 4)
        self.assertEqual(saturation.ISOLATION_FRACTION, 0.3)
        self.assertEqual(saturation.SPOT_RADIUS, 5.0)


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
        config = dict(clipLevel=CLIP, thresholdFactor=1.0, maxFiducialFraction=0.5)
        saturated = saturation.measureConfigured(config, canonWithSpot(), EXPTIME)
        self.assertEqual(saturated.reason, '')
        self.assertEqual(int(saturated.mask.sum()), 9)

    def testNoConfiguration(self):
        saturated = saturation.measureConfigured(None, canonWithSpot(), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertIn('no saturation configuration', saturated.reason)

    def testMissingSetting(self):
        """A configuration without one of its settings gives no mask and names the setting."""
        with self.assertLogs('saturation', level='ERROR'):
            saturated = saturation.measureConfigured(dict(clipLevel=CLIP), canonWithSpot(), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertIn('thresholdFactor', saturated.reason)

    def testConfigurationOfAnotherType(self):
        for config in ('2815', [CLIP, 1.0], dict(clipLevel='high', thresholdFactor=1.0)):
            with self.assertLogs('saturation', level='ERROR'):
                saturated = saturation.measureConfigured(config, canonWithSpot(), EXPTIME)
            self.assertIsNone(saturated.mask)
            self.assertTrue(saturated.reason)

    def testSettingOutOfRange(self):
        """A clip level no read can reach leaves the saturation unknown, it does not report none."""
        config = dict(clipLevel=28150.0, thresholdFactor=1.0)
        saturated = saturation.measureConfigured(config, canonWithSpot(value=60000), EXPTIME)
        self.assertIsNone(saturated.mask)
        self.assertIn('cannot reach', saturated.reason)


class SpotFlagsTestCase(unittest.TestCase):

    def setUp(self):
        self.mask = np.zeros((NROWS, NCOLS), dtype=bool)
        self.mask[20, 30] = True

    def testWithinTheRadius(self):
        x = np.array([30.0, 34.9, 35.2, 30.0, 33.0])
        y = np.array([20.0, 20.0, 20.0, 26.0, 24.0])
        self.assertEqual(saturation.spotsOnMask(self.mask, x, y).tolist(), [True, True, False, False, True])

    def testNonFinitePosition(self):
        flagged = saturation.spotsOnMask(self.mask, np.array([np.nan, 30.0]), np.array([20.0, np.nan]))
        self.assertEqual(flagged.tolist(), [False, False])

    def testSpotOffTheImage(self):
        self.mask[0, 0] = True
        flagged = saturation.spotsOnMask(self.mask, np.array([-2.0, -20.0, 500.0]), np.array([1.0, 1.0, 1.0]))
        self.assertEqual(flagged.tolist(), [True, False, False])

    def testNoSpot(self):
        self.assertEqual(saturation.spotsOnMask(self.mask, np.array([]), np.array([])).shape, (0,))
        self.assertEqual(saturation.spotFlags(self.mask, np.array([]), np.array([])).shape, (0,))

    def testFlags(self):
        flags = saturation.spotFlags(self.mask, np.array([31.0, 60.0]), np.array([21.0, 50.0]))
        self.assertEqual(flags.tolist(), [SpotFlag.SATURATED, 0])

    def testWithoutAMask(self):
        """Every spot of an image that could not be measured is unknown, none is clean."""
        flags = saturation.spotFlags(None, np.array([31.0, 60.0]), np.array([21.0, 50.0]))
        self.assertEqual(flags.tolist(), [UNKNOWN, UNKNOWN])


class FlagsOfSpotsTestCase(unittest.TestCase):

    def setUp(self):
        self.flags = np.array([0, 1, 0, 1], dtype=np.int32)
        self.spotIds = np.array([1.0, 2.0, 3.0, 4.0])

    def testWanted(self):
        """Flags follow the wanted identifiers, in their order; those that are not positive are left out."""
        flags = saturation.flagsOfSpots(self.flags, self.spotIds, np.array([4.0, -1.0, 1.0, 2.0, 0.0]))
        self.assertEqual(flags.tolist(), [1, 0, 1])

    def testSpotIdsInAnyOrder(self):
        flags = saturation.flagsOfSpots(self.flags, np.array([40, 10, 30, 20]), np.array([10, 20, 30, 40]))
        self.assertEqual(flags.tolist(), [1, 1, 0, 0])

    def testIdentifierThatIsNotASpot(self):
        """A wanted spot that is not among the spots is unknown, whether below, between or above them."""
        flags = saturation.flagsOfSpots(self.flags, np.array([2, 4, 6, 8]), np.array([1, 5, 9, 4]))
        self.assertEqual(flags.tolist(), [UNKNOWN, UNKNOWN, UNKNOWN, 1])

    def testNoSpot(self):
        flags = saturation.flagsOfSpots(np.array([], dtype=int), np.array([]), np.array([3, -1]))
        self.assertEqual(flags.tolist(), [UNKNOWN])
        self.assertEqual(saturation.flagsOfSpots(self.flags, self.spotIds, np.array([-1.0, -1.0])).shape, (0,))


class ReportsTestCase(unittest.TestCase):

    def testCountSaturated(self):
        self.assertEqual(saturation.countSaturated(np.array([0, 1, 1, 0])), 2)
        self.assertEqual(saturation.countSaturated(np.array([], dtype=int)), 0)
        self.assertIsNone(saturation.countSaturated(np.array([0, UNKNOWN])))

    def testCountKeyword(self):
        self.assertEqual(saturation.countKeyword('saturatedSpots', 15272800, np.array([0, 1, 1, 0])),
                         'saturatedSpots=15272800,2,4')

    def testCountKeywordWhenUnknown(self):
        """An unknown saturation is reported as -1, never as no saturated spot."""
        self.assertEqual(saturation.countKeyword('saturatedCobras', 15272800, np.full(3, UNKNOWN)),
                         'saturatedCobras=15272800,-1,3')

    def testExceedsFraction(self):
        flags = np.array([1, 1, 0, 0])
        self.assertFalse(saturation.exceedsFraction(flags, 0.5))
        self.assertTrue(saturation.exceedsFraction(flags, 0.4))
        self.assertTrue(saturation.exceedsFraction(flags, '0.4'))
        self.assertFalse(saturation.exceedsFraction(flags, None))
        self.assertFalse(saturation.exceedsFraction(np.array([], dtype=int), 0.0))

    def testUnknownIsNotAnExcess(self):
        self.assertFalse(saturation.exceedsFraction(np.array([1, UNKNOWN]), 0.0))


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
        """A saturated pixel has the bit of the SAT plane, which the MP_SAT card names."""
        mask = np.zeros((NROWS, NCOLS), dtype=bool)
        mask[20, 30:32] = True
        for data, header in self.written(saturation.maskHdu(mask.shape, mask, '')):
            self.assertEqual(data.shape, mask.shape)
            np.testing.assert_array_equal(data != 0, mask)
            self.assertEqual(int(data[20, 30]), 1 << int(header['MP_SAT']))
            self.assertTrue(header['SATKNOWN'])
            self.assertNotIn('SATWHY', header)

    def testUnknown(self):
        """Without a measurement no bit is set and the header says so, with the reason."""
        reason = 'a 4.8 s exposure sums 7 reads but the background (1283 ADU) is the bias of 5.0'
        for data, header in self.written(saturation.maskHdu((NROWS, NCOLS), None, reason)):
            self.assertEqual(data.shape, (NROWS, NCOLS))
            self.assertFalse(data.any())
            self.assertFalse(header['SATKNOWN'])
            self.assertEqual(header['SATWHY'], reason)

    def testReasonACardCannotHold(self):
        """A reason with a line break or a non-ASCII character is still written."""
        for data, header in self.written(saturation.maskHdu((NROWS, NCOLS), None, 'KeyError: \u00e9t\u00e9\nline two')):
            self.assertFalse(header['SATKNOWN'])
            self.assertEqual(header['SATWHY'], 'KeyError:  t  line two')


if __name__ == '__main__':
    unittest.main()
