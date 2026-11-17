"""The spotFinder library: column levels and noise, detection, windowed centroids, fluxes."""
import unittest

import numpy as np

import mcsActor.spotFinder.spotFinder as spotFinder

SHAPE = (1400, 2000)
LOWER_PERCENTILE = 15.865525393145707
PARAMS = dict(findSigma=25, nmin=10, boxFind=14, boxCent=14, fwhmx=4, fwhmy=9, maxIt=20, aperture=10, gain=2.24)


def noisyFrame(shape, seed, *, noise=6.0, level=1800):
    """int32 image of Gaussian noise about a level that differs from column to column."""
    rng = np.random.default_rng(seed)
    columns = level + rng.integers(-12, high=13, size=shape[1])
    return np.rint(columns + rng.normal(loc=0, scale=noise, size=shape)).astype('<i4')


def expectedNoise(image):
    """Root mean square over the columns of the median minus the lower percentile, with numpy."""
    lower = np.median(image, axis=0) - np.percentile(image, LOWER_PERCENTILE, axis=0)
    return np.sqrt(np.mean(lower ** 2))


def subtracted(image):
    """The image less floor(level + 0.5) of each column, with numpy."""
    return image.astype('<i4') - np.floor(np.median(image, axis=0) + 0.5).astype('<i4')


def addSpots(image, positions, *, amplitude=5000.0, sigma=2.0, lobe=0.0, lobeOffset=(0.0, -4.0)):
    """Add Gaussian spots to an image, in place, each with an optional second, displaced Gaussian.

    Parameters
    ----------
    positions : `numpy.ndarray`, shape (n, 2)
        x, y of the spots; a pixel is centred on its indices.
    amplitude : `float` or `numpy.ndarray`, shape (n,)
    lobe : `float`
        Amplitude of the second Gaussian, as a fraction of the first.
    lobeOffset : `tuple` of `float`
        Its displacement along x and y (pix).
    """
    half = 20
    amplitudes = np.broadcast_to(amplitude, len(positions))
    for (x, y), height in zip(positions, amplitudes):
        i0, j0 = int(round(x)), int(round(y))
        rows, cols = np.mgrid[j0 - half:j0 + half + 1, i0 - half:i0 + half + 1]
        spot = np.exp(-((cols - x) ** 2 + (rows - y) ** 2) / (2 * sigma ** 2))
        spot += lobe * np.exp(-((cols - x - lobeOffset[0]) ** 2 + (rows - y - lobeOffset[1]) ** 2) / (2 * sigma ** 2))
        image[j0 - half:j0 + half + 1, i0 - half:i0 + half + 1] += np.rint(height * spot).astype(image.dtype)


def gridPositions(seed, *, jitterX=0.5, jitterY=0.5):
    """Spot positions 100 pixels apart or more, at random sub-pixel offsets."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(60, 1950, step=120), np.arange(60, 1350, step=100))
    positions = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    return positions + np.column_stack([rng.uniform(low=-jitterX, high=jitterX, size=len(positions)),
                                        rng.uniform(low=-jitterY, high=jitterY, size=len(positions))])


def spotFrame(positions, **kwargs):
    """int32 frame of zeros with spots at positions."""
    image = np.zeros(SHAPE, dtype='<i4')
    addSpots(image, positions, **kwargs)
    return image


def cameraFrame(positions, seed, *, level=1795, noise=5.0, amplitude=4000.0):
    """uint16 frame like the camera's: a level per column, every eighth one lower and noisier, noise, spots."""
    rng = np.random.default_rng(seed)
    columns = level + rng.integers(-6, high=7, size=SHAPE[1])
    columns[::8] -= 8
    columns[::5] += 40
    pixels = rng.normal(loc=0, scale=noise, size=SHAPE)
    pixels[:, ::8] *= 3
    image = np.rint(columns + pixels).astype('<i4')
    addSpots(image, positions, amplitude=amplitude)
    return image.astype('<u2')


def find(image, **changes):
    """findSpots with PARAMS, changed as given."""
    return spotFinder.findSpots(image, **{**PARAMS, **changes})


def distances(spots, positions):
    """Distance (pix) from each position to the nearest spot."""
    return np.hypot(spots['x'][None, :] - positions[:, 0][:, None],
                    spots['y'][None, :] - positions[:, 1][:, None]).min(axis=1)


def windowedStep(image, x, y, *, box=14, fwhmx=4.0, fwhmy=9.0):
    """Twice the Gaussian-weighted mean offset of the pixels from (x, y): zero where the windowed centroid settles."""
    rows, cols = np.mgrid[:image.shape[0], :image.shape[1]]
    radius2 = (cols - x) ** 2 + (rows - y) ** 2
    inside = radius2 < box ** 2
    values, radius2 = image[inside].astype(float), radius2[inside]
    weightX = np.exp(-radius2 / (2 * (fwhmx / np.sqrt(8 * np.log(2))) ** 2))
    weightY = np.exp(-radius2 / (2 * (fwhmy / np.sqrt(8 * np.log(2))) ** 2))
    return (2 * np.sum(weightX * values * (cols[inside] - x)) / np.sum(weightX * values),
            2 * np.sum(weightY * values * (rows[inside] - y)) / np.sum(weightY * values))


class ColumnLevelsTestCase(unittest.TestCase):

    def testLevelsAreColumnMedians(self):
        """The levels are numpy's column medians, for an odd and an even number of rows."""
        for shape in ((501, 300), (500, 300)):
            image = noisyFrame(shape, seed=1)
            levels, _ = spotFinder.columnLevels(image)

            np.testing.assert_array_equal(levels, np.median(image, axis=0))

        self.assertTrue(np.any(levels % 1 != 0), msg='an even number of rows should give some half-integer medians')

    def testNoiseIsLowerHalfOfColumns(self):
        """The noise is the rms over columns of median minus 15.87th percentile, and close to the true sigma."""
        image = noisyFrame((2000, 400), seed=2, noise=6.0)
        _, noise = spotFinder.columnLevels(image)

        self.assertAlmostEqual(noise, expectedNoise(image), places=9)
        self.assertAlmostEqual(noise, 6.0, delta=0.3)

    def testColumnsOutsideTheHistogram(self):
        """Columns far from the level of the image, or split between two levels, still get their median."""
        image = noisyFrame((600, 200), seed=4)
        image[:, 10] = 40000
        image[:, 20] = 0
        image[::2, 30] += 5000
        image[:400, 40] -= 3000
        image[:50, 50] += 20000
        levels, noise = spotFinder.columnLevels(image)

        np.testing.assert_array_equal(levels, np.median(image, axis=0))
        self.assertAlmostEqual(noise, expectedNoise(image), places=9)

    def testExtremeValues(self):
        """Pixels at the ends of the range handled do not disturb the columns they are in."""
        image = noisyFrame((400, 100), seed=5)
        image[0, 5] = spotFinder.VALUE_LIMIT
        image[1, 5] = -spotFinder.VALUE_LIMIT
        levels, _ = spotFinder.columnLevels(image)

        np.testing.assert_array_equal(levels, np.median(image, axis=0))

    def testFramesRefused(self):
        """A frame with a pixel beyond the range handled, or not of integers, is refused."""
        for dtype, value in (('<i4', spotFinder.VALUE_LIMIT + 1), ('<i8', 2 ** 32 + 3000), ('<i8', -2 ** 31 - 5),
                             ('<u4', 2 ** 31 + 7)):
            image = noisyFrame((400, 100), seed=5).astype(dtype)
            image[3, 7] = value
            for function in (spotFinder.columnLevels, find, lambda frame: spotFinder.apertureFlux(frame, [50], [200], 5)):
                with self.assertRaises(ValueError):
                    function(image)
        for function in (spotFinder.columnLevels, find):
            with self.assertRaises(TypeError):
                function(noisyFrame((400, 100), seed=5).astype('<f8'))

    def testLightLeavesTheNoise(self):
        """Light on a twentieth of the pixels moves the noise by a few percent, their standard deviation far more."""
        rng = np.random.default_rng(6)
        dark = noisyFrame((2000, 400), seed=6, noise=6.0)
        lit = dark.copy()
        lit[rng.random(dark.shape) < 0.05] += 3000
        _, darkNoise = spotFinder.columnLevels(dark)
        _, litNoise = spotFinder.columnLevels(lit)

        self.assertAlmostEqual(litNoise / darkNoise, 1.0, delta=0.05)
        self.assertGreater((lit - np.median(lit, axis=0)).std(), 50 * darkNoise)

    def testFewRows(self):
        """One, two and three rows give numpy's medians and noise."""
        for nRows in (1, 2, 3):
            image = noisyFrame((nRows, 40), seed=8)
            levels, noise = spotFinder.columnLevels(image)

            np.testing.assert_array_equal(levels, np.median(image, axis=0))
            self.assertAlmostEqual(noise, expectedNoise(image), places=9)

    def testConstantImage(self):
        """An image of one value has that level and no noise."""
        levels, noise = spotFinder.columnLevels(np.full((300, 50), 1795, dtype='<u2'))

        np.testing.assert_array_equal(levels, 1795)
        self.assertEqual(noise, 0.0)

    def testThreadsDoNotChangeTheLevels(self):
        """Whatever the number of threads, the same levels and noise; more threads than columns too."""
        image = noisyFrame((700, 37), seed=9)
        reference = spotFinder.columnLevels(image, nThreads=1)
        for nThreads in (2, 3, 8, 64, 1000):
            levels, noise = spotFinder.columnLevels(image, nThreads=nThreads)

            np.testing.assert_array_equal(levels, reference[0])
            self.assertEqual(noise, reference[1])

    def testSizesRefused(self):
        """An image with no pixel, or with more rows than a histogram bin can count, is refused."""
        for shape in ((65536, 3), (0, 5), (5, 0)):
            with self.assertRaises(ValueError):
                spotFinder.columnLevels(np.zeros(shape, dtype='<i4'))
            with self.assertRaises(ValueError):
                find(np.zeros(shape, dtype='<i4'))


class ApertureFluxTestCase(unittest.TestCase):

    def expected(self, image, x, y, radius):
        """Sum of the pixels whose centre is within radius of (x, y), with numpy."""
        rows, cols = np.mgrid[:image.shape[0], :image.shape[1]]
        return float(image[(cols - x) ** 2 + (rows - y) ** 2 <= radius ** 2].sum())

    def testSumsPixelsWithinRadius(self):
        """The flux is the sum of the pixels whose centre lies within the radius, the limit included."""
        rng = np.random.default_rng(7)
        image = rng.integers(-50, high=500, size=(80, 120)).astype('<i4')
        x = np.array([30.0, 41.37, 60.5, 90.49, 25.0])
        y = np.array([40.0, 20.81, 33.5, 50.02, 25.0])
        for radius in (3.0, 8.0, 10.0, 5.0):
            flux = spotFinder.apertureFlux(image, x, y, radius)

            np.testing.assert_allclose(flux, [self.expected(image, a, b, radius) for a, b in zip(x, y)])

    def testNoFluxOffTheImage(self):
        """A position that is not on the image, or whose aperture leaves it, gets NaN."""
        image = np.ones((60, 80), dtype='<i4')
        x = np.array([np.nan, 40.0, 4.0, 76.0, 40.0, 40.0, 1e12, -1e12, 40.0, 5.0])
        y = np.array([30.0, np.inf, 30.0, 30.0, 4.0, 56.0, 30.0, 30.0, -1e300, 5.0])
        flux = spotFinder.apertureFlux(image, x, y, 5.0)

        self.assertTrue(np.all(np.isnan(flux[:9])))
        self.assertEqual(flux[9], self.expected(image, 5.0, 5.0, 5.0))

    def testPositionsMustPair(self):
        """As many x as y."""
        with self.assertRaises(ValueError):
            spotFinder.apertureFlux(np.ones((60, 80), dtype='<i4'), np.zeros(3), np.zeros(2), 5.0)


class FindSpotsTestCase(unittest.TestCase):

    def testPositions(self):
        """Every spot is found once, at its position, in the order of the frame."""
        positions = gridPositions(seed=1)
        spots, _ = find(spotFrame(positions))

        self.assertEqual(len(spots), len(positions))
        self.assertLess(distances(spots, positions).max(), 0.01)
        self.assertTrue(np.all(np.diff(np.rint(spots['y'] / 50)) >= 0))

    def testAnyNumberOfSpots(self):
        """No spot, one, a few, wherever they are on the frame."""
        corners = np.array([(60.3, 70.2), (1930.7, 61.4), (64.1, 1330.6), (1935.2, 1333.9), (1000.5, 700.5)])
        for nSpots in range(len(corners) + 1):
            for first in range(len(corners)):
                positions = np.roll(corners, -first, axis=0)[:nSpots]
                with self.subTest(positions=positions.tolist()):
                    spots, _ = find(spotFrame(positions))

                    self.assertEqual(len(spots), nSpots)
                    if nSpots:
                        self.assertLess(distances(spots, positions).max(), 0.01)

    def testLevels(self):
        """Background, noise and threshold are those of the columns of the frame."""
        image = cameraFrame(gridPositions(seed=2), seed=2)
        _, levels = find(image, findSigma=20)

        self.assertEqual(levels['background'], np.median(np.median(image, axis=0)))
        self.assertAlmostEqual(levels['noise'], expectedNoise(image.astype('<i4')), places=9)
        self.assertAlmostEqual(levels['noise'], np.sqrt((7 * 5.0 ** 2 + 15.0 ** 2) / 8), delta=0.7)
        self.assertEqual(levels['threshold'], np.floor(20 * levels['noise']))
        self.assertNotEqual(levels['threshold'], 20 * levels['noise'])

    def testCameraFrame(self):
        """On a frame like the camera's, every spot is found, close to its position."""
        positions = gridPositions(seed=3)
        spots, _ = find(cameraFrame(positions, seed=3))

        self.assertEqual(len(spots), len(positions))
        self.assertLess(distances(spots, positions).max(), 0.05)

    def testLevelDoesNotMoveTheSpots(self):
        """The same frame 2,000 ADU higher gives the same spots, on backgrounds 2,000 ADU higher."""
        low = cameraFrame(gridPositions(seed=4), seed=4)
        high = (low.astype('<i4') + 2000).astype('<u2')
        spotsLow, levelsLow = find(low)
        spotsHigh, levelsHigh = find(high)
        for name in spotsLow.dtype.names:
            if name != 'background':
                np.testing.assert_array_equal(spotsHigh[name], spotsLow[name], err_msg=name)
        np.testing.assert_array_equal(spotsHigh['background'], spotsLow['background'] + 2000)
        self.assertEqual(levelsHigh, dict(levelsLow, background=levelsLow['background'] + 2000))

    def testBackgroundOfEachSpot(self):
        """The background of a spot is the median of the column of its centroid, rounded; with the peak, the raw peak."""
        image = cameraFrame(gridPositions(seed=19), seed=19)
        spots, levels = find(image)
        columns = np.floor(np.median(image, axis=0) + 0.5)
        cols, rows = np.rint(spots['x']).astype(int), np.rint(spots['y']).astype(int)
        rawPeaks = np.array([image[j - 3:j + 4, i - 3:i + 4].max() for i, j in zip(cols, rows)])

        np.testing.assert_array_equal(spots['background'], columns[cols])
        self.assertGreater(len(np.unique(spots['background'])), 3)
        np.testing.assert_allclose(spots['peak'] + spots['background'], rawPeaks, atol=60)

    def testFrameUntouched(self):
        """The frame given is not modified, whatever its type."""
        positions = gridPositions(seed=5)
        for dtype in ('<u2', '<i4', '>u2', '<i8'):
            image = cameraFrame(positions, seed=5).astype(dtype)
            original = image.copy()
            spots, _ = find(image)

            np.testing.assert_array_equal(image, original)
            self.assertEqual(len(spots), len(positions))

    def testFluxIsApertureSum(self):
        """The flux is the sum, within the aperture of the centroid, of the frame less its column levels."""
        image = cameraFrame(gridPositions(seed=6), seed=6)
        for aperture in (10, 6.5):
            spots, _ = find(image, aperture=aperture)

            np.testing.assert_array_equal(spots['flux'],
                                          spotFinder.apertureFlux(subtracted(image), spots['x'], spots['y'], aperture))
        np.testing.assert_allclose(find(image)[0]['flux'], 2 * np.pi * 2.0 ** 2 * 4000.0, rtol=0.03)

    def testFluxError(self):
        """The error is the noise of the pixels summed and the photon noise of the flux, in quadrature."""
        image = cameraFrame(gridPositions(seed=7), seed=7)
        spots, levels = find(image, aperture=10, gain=3.0)
        nPixels = spotFinder.apertureFlux(np.ones(SHAPE, dtype='<i4'), spots['x'], spots['y'], 10)

        self.assertTrue(np.all((nPixels > 300) & (nPixels < 330)))
        np.testing.assert_allclose(spots['fluxerr'], np.sqrt(nPixels * levels['noise'] ** 2 + spots['flux'] / 3.0))

    def testNoFluxWithoutAperture(self):
        """No aperture, no flux and no error."""
        spots, _ = find(spotFrame(gridPositions(seed=8)), aperture=0)

        self.assertTrue(np.all(np.isnan(spots['flux'])) and np.all(np.isnan(spots['fluxerr'])))

    def testPeak(self):
        """The peak is the highest pixel of the spot, above the level of its column."""
        image = cameraFrame(gridPositions(seed=9), seed=9)
        spots, _ = find(image)
        flat = subtracted(image)
        cols, rows = np.rint(spots['x']).astype(int), np.rint(spots['y']).astype(int)

        np.testing.assert_array_equal(spots['peak'], [flat[j - 3:j + 4, i - 3:i + 4].max() for i, j in zip(cols, rows)])
        self.assertTrue(np.all((spots['peak'] > 3200) & (spots['peak'] < 4100)))

    def testMoments(self):
        """A spot stretched along y has the larger second moment there, and no cross moment."""
        spots, _ = find(spotFrame(gridPositions(seed=10), lobe=0.6, lobeOffset=(0.0, -3.0)))

        np.testing.assert_allclose(spots['x2'], 2.0 ** 2, rtol=0.1)
        self.assertTrue(np.all(spots['y2'] > 1.4 * spots['x2']))
        self.assertLess(np.abs(spots['xy']).max(), 0.1)

    def testCrossMoment(self):
        """A spot stretched along a diagonal has a cross moment of the sign of its slope."""
        rising, _ = find(spotFrame(gridPositions(seed=10), lobe=0.6, lobeOffset=(3.0, 3.0)))
        falling, _ = find(spotFrame(gridPositions(seed=10), lobe=0.6, lobeOffset=(3.0, -3.0)))

        self.assertTrue(np.all(rising['xy'] > 0.5))
        self.assertTrue(np.all(falling['xy'] < -0.5))
        np.testing.assert_allclose(rising['x2'], rising['y2'], rtol=0.05)

    def testBothAxesConverge(self):
        """A spot symmetric along one axis and lopsided along the other is iterated until both have settled."""
        for lobeOffset, jitter in (((0.0, -4.0), dict(jitterX=0.0)), ((-4.0, 0.0), dict(jitterY=0.0))):
            image = spotFrame(gridPositions(seed=11, **jitter), lobe=0.3, lobeOffset=lobeOffset)
            spots, _ = find(image)

            for spot in spots:
                i0, j0 = int(spot['x']) - 40, int(spot['y']) - 40
                dx, dy = windowedStep(image[j0:j0 + 81, i0:i0 + 81], spot['x'] - i0, spot['y'] - j0)

                self.assertLess(abs(dx), 2e-4)
                self.assertLess(abs(dy), 2e-4)

    def testWindowIsCircular(self):
        """Light beyond boxCent of the centroid does not move it, even within boxCent along each axis."""
        positions = gridPositions(seed=16)
        alone = spotFrame(positions)
        company = alone.copy()
        for x, y in positions:
            i, j = int(round(x)) + 11, int(round(y)) + 11
            company[j:j + 2, i:i + 2] = 3000
        spotsAlone, _ = find(alone)
        spotsCompany, _ = find(company)

        np.testing.assert_array_equal(spotsCompany['x'], spotsAlone['x'])
        np.testing.assert_array_equal(spotsCompany['y'], spotsAlone['y'])

    def testIterationLimitNotReached(self):
        """Spots two pixels wide settle within the 20 iterations of the actor."""
        image = spotFrame(gridPositions(seed=12), lobe=0.3)
        short, _ = find(image, maxIt=20)
        long, _ = find(image, maxIt=200)

        self.assertEqual(short.tobytes(), long.tobytes())

    def testWideWindowSettlesQuickly(self):
        """With a window three times wider than the spots, a few iterations reach the converged position."""
        image = spotFrame(gridPositions(seed=12), lobe=0.3)
        wide = dict(fwhmx=14.0, fwhmy=14.0, boxCent=20)
        short, _ = find(image, maxIt=6, **wide)
        long, _ = find(image, maxIt=500, **wide)

        np.testing.assert_allclose(short['x'], long['x'], rtol=0, atol=1e-3)
        np.testing.assert_allclose(short['y'], long['y'], rtol=0, atol=1e-3)

    def testThreadsDoNotChangeTheSpots(self):
        """Whatever the number of threads, the same spots."""
        image = cameraFrame(gridPositions(seed=13), seed=13)
        reference, levels = find(image, nThreads=1)
        for nThreads in (2, 5, 16, 1000):
            spots, other = find(image, nThreads=nThreads)

            self.assertEqual(spots.tobytes(), reference.tobytes())
            self.assertEqual(other, levels)

    def testThresholdOfAFrameWithoutNoise(self):
        """A frame without noise is searched at 1 ADU, not at zero."""
        positions = gridPositions(seed=14)
        image = np.full(SHAPE, 1795, dtype='<u2')
        image += spotFrame(positions).astype('<u2')
        spots, levels = find(image)

        self.assertEqual((levels['noise'], levels['threshold']), (0.0, 1.0))
        self.assertEqual(len(spots), len(positions))

    def testMinimumSize(self):
        """A group of pixels above the threshold smaller than nmin is not a spot."""
        image = np.zeros(SHAPE, dtype='<i4')
        image[700:703, 1000:1003] = 500
        image[400:404, 600:604] = 500

        self.assertEqual(len(find(image, nmin=10)[0]), 1)
        self.assertEqual(len(find(image, nmin=9)[0]), 2)
        self.assertEqual(len(find(image, nmin=17)[0]), 0)

    def testSpotStartsAboveTheThreshold(self):
        """A spot needs a pixel above the threshold, and then takes the pixels at the threshold beside it."""
        image = np.zeros(SHAPE, dtype='<i4')
        image[700:704, 1000:1004] = 1
        image[400:404, 600:604] = 1
        image[401, 601] = 2
        spots, levels = find(image, nmin=16, boxCent=0, aperture=0)

        self.assertEqual(levels['threshold'], 1.0)
        self.assertEqual(len(spots), 1)
        np.testing.assert_allclose([spots['x'][0], spots['y'][0]], [601.5 - 0.5 / 17, 401.5 - 0.5 / 17])

    def testSidesConnectNotCorners(self):
        """Pixels touching by a corner only are two spots; by a side, one."""
        byCorner = np.zeros(SHAPE, dtype='<i4')
        byCorner[700:704, 1000:1004] = 500
        byCorner[704:708, 1004:1008] = 500
        bySide = np.zeros(SHAPE, dtype='<i4')
        bySide[700:704, 1000:1004] = 500
        bySide[704:708, 1003:1007] = 500

        self.assertEqual(len(find(byCorner, boxFind=1, boxCent=1, aperture=0)[0]), 2)
        self.assertEqual(len(find(bySide, boxFind=1, boxCent=1, aperture=0)[0]), 1)

    def testEndsOfRowsDoNotConnect(self):
        """The last pixel of a row and the first of the next are not neighbours, whichever is looked at first."""
        nx = SHAPE[1]
        rightFirst = np.zeros(SHAPE, dtype='<i4')
        rightFirst[100:104, nx - 3:] = 5
        rightFirst[101:105, :3] = 5
        leftFirst = np.zeros(SHAPE, dtype='<i4')
        leftFirst[99, nx - 3:] = 1
        leftFirst[100:103, nx - 3:] = 5
        leftFirst[100:104, :3] = 5

        for image in (rightFirst, leftFirst):
            spots, _ = find(image, boxFind=0, boxCent=0, aperture=0)

            self.assertEqual(len(spots), 2)
            np.testing.assert_allclose(np.sort(spots['x']), [1, nx - 2])

    def testEdges(self):
        """No spot is kept closer to an edge than three times the larger of the two boxes."""
        inside = np.array([(50.3, 700.2), (1000.4, 49.6), (1949.2, 700.3), (1000.1, 1349.7),
                           (43.3, 500.2), (800.4, 43.3), (1955.8, 500.3), (800.1, 1355.8)])
        outside = np.array([(40.6, 300.2), (500.4, 40.6), (1958.5, 300.3), (500.1, 1358.5)])
        for boxes in (dict(boxFind=14, boxCent=14), dict(boxFind=5, boxCent=14), dict(boxFind=14, boxCent=5)):
            spots, _ = find(spotFrame(np.vstack([inside, outside])), **boxes)

            self.assertEqual(len(spots), len(inside), msg=boxes)
            self.assertLess(distances(spots, inside).max(), 0.05, msg=boxes)

    def testRegionsSettlingOnOneSpot(self):
        """A faint group of pixels whose centroid settles on a bright spot nearby is that spot, counted once."""
        alone = (1400.6, 700.2)
        for bright, faint in (((1000.3, 700.4), (1000.3, 712.4)), ((1000.3, 712.4), (1000.3, 700.4))):
            image = spotFrame(np.array([bright, faint, alone]), amplitude=np.array([5000.0, 400.0, 400.0]), sigma=1.2)
            spots, _ = find(image)
            brightPixel = image[int(round(bright[1])), int(round(bright[0]))]

            self.assertEqual(len(spots), 2)
            self.assertLess(distances(spots, np.array([bright, alone])).max(), 0.02)
            np.testing.assert_array_equal(np.sort(spots['peak']), [image[700, 1401], brightPixel])

    def testClosestSpots(self):
        """Two spots are one when their centroids are closer than boxFind + 1 pixels, the brighter one."""
        for separation, expected in ((14.5, 1), (15.5, 2)):
            for brighter in (0, 1):
                positions = np.array([(1000.0, 700.0), (1000.0 + separation, 700.0)])
                amplitude = np.where(np.arange(2) == brighter, 5000.0, 3000.0)
                image = spotFrame(positions, amplitude=amplitude, sigma=1.2)
                spots, _ = find(image, boxCent=3)
                brighterPeak = image[690:711, int(positions[brighter, 0]) - 2:int(positions[brighter, 0]) + 3].max()

                self.assertEqual(len(spots), expected, msg=(separation, brighter))
                self.assertEqual(spots['peak'].max(), brighterPeak)
                if expected == 1:
                    self.assertLess(abs(spots['x'][0] - positions[brighter, 0]), 0.01)

    def testFloodedFrame(self):
        """A frame with millions of pixels above the threshold is not searched."""
        rng = np.random.default_rng(15)
        image = noisyFrame((2600, 3000), seed=15)
        image[rng.random(image.shape) < 0.45] += 3000

        with self.assertRaises(RuntimeError):
            find(image)

    def testTooManySpots(self):
        """A frame with spots by the tens of thousands is given up."""
        image = np.zeros(SHAPE, dtype='<i4')
        for dy in range(4):
            for dx in range(4):
                image[50 + dy:1350:8, 50 + dx:1950:8] = 500

        with self.assertRaises(RuntimeError):
            find(image)


class WindowedCentroidTestCase(unittest.TestCase):

    def matched(self, image):
        """Spots of findSpots and of the windowed centroid on a frame without background, in the same order."""
        import mcsActor.windowedCentroid.centroid as centroid

        spots, _ = find(image, findSigma=0, aperture=0)
        old = np.frombuffer(centroid.centroid_only(image, 4, 9, 1, 1, 14, 14, 10, 20, 0), dtype='<f8').reshape(-1, 7)
        nearest = np.hypot(old[:, 0][None, :] - spots['x'][:, None], old[:, 1][None, :] - spots['y'][:, None]).argmin(axis=1)

        self.assertEqual(len(old), len(spots))
        self.assertEqual(len(set(nearest)), len(spots))
        return spots, old[nearest]

    def testSameSpotsAsTheWindowedCentroid(self):
        """The spots, positions, moments and peaks of symmetric spots are those of the windowed centroid."""
        spots, old = self.matched(spotFrame(gridPositions(seed=17)))

        np.testing.assert_allclose(spots['x'], old[:, 0], atol=1e-3)
        np.testing.assert_allclose(spots['y'], old[:, 1], atol=1e-3)
        np.testing.assert_allclose(spots['x2'], old[:, 2], rtol=1e-9)
        np.testing.assert_allclose(spots['y2'], old[:, 3], rtol=1e-9)
        np.testing.assert_array_equal(spots['peak'], old[:, 4])
        np.testing.assert_allclose(spots['xy'], old[:, 5], atol=1e-9)

    def testLopsidedSpotsSettleFurther(self):
        """On spots lopsided along y, y settles where the windowed centroid, stopping with x, leaves it short."""
        image = spotFrame(gridPositions(seed=18, jitterX=0.0), lobe=0.3)
        spots, old = self.matched(image)

        np.testing.assert_allclose(spots['x'], old[:, 0], atol=1e-3)
        steps = []
        for (x, y), (xOld, yOld) in zip(np.column_stack([spots['x'], spots['y']]), old[:, :2]):
            i0, j0 = int(x) - 40, int(y) - 40
            stamp = image[j0:j0 + 81, i0:i0 + 81]
            steps.append((windowedStep(stamp, x - i0, y - j0)[1], windowedStep(stamp, xOld - i0, yOld - j0)[1]))
        steps = np.abs(np.array(steps))

        self.assertLess(steps[:, 0].max(), 2e-4)
        self.assertGreater(np.median(steps[:, 1]), 10 * steps[:, 0].max())


if __name__ == '__main__':
    unittest.main()
