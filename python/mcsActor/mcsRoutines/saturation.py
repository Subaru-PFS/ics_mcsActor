"""Saturated pixels and spots of a Canon MCS coadd.

``canonexp`` sums n = ceil(T / 0.8 s) + 1 reads into an unsigned 16-bit coadd.  Each read starts from
an empty pixel, carries a bias and clips at a level above it, so a pixel saturates when its most exposed
read reaches the clip level, at a coadd value of

    n bias + clip (seconds of the exposure the row captured) / (seconds in its most exposed read).

The camera free-runs with a rolling readout that is not synchronised with the shutter: row y is read
ROW_TIME y after row 0, and how the exposure is shared between the reads of a row depends on the delay
(the phase) between the shutter opening and the first read of row 0.  The phase changes from frame to
frame and is not recorded, so the threshold of a row is the lowest value the expression takes over the
phases a frame can have.

The clip level drifts with time and is configuration, as is the factor applied to it.
"""
import enum
import functools
import logging
import math
import typing

import astropy.io.fits as pyfits
import numpy as np

CANON_SHAPE = (5778, 8960)
"""(rows, columns) of a Canon coadd."""

BIAS_PER_READ = 256.5
"""Bias of one read (ADU)."""

LOOP_TIME = 0.8
"""Divisor ``canonexp`` uses to count the reads of an exposure (s)."""

READ_PERIOD = 0.826
"""Time between two reads of the same row (s)."""

ROW_TIME = 1.43e-4
"""Delay of the rolling readout between consecutive rows (s); rows are read in increasing y."""

ROW_TIME_MAX = 1.47e-4
"""Upper bound of the row time (s); its lower bound is taken symmetric about ROW_TIME."""

PHASE_MIN = -0.20
"""Lowest phase of a frame (s): minus the longest delay from the start of the acquisition to the
shutter opening."""

PHASE_MAX = READ_PERIOD + 0.05
"""Highest phase of a frame (s): one read period, plus an allowance for the delay between the start
of the acquisition and the first delivered read."""

PHASE_STEP = 1e-3
"""Step (s) of the phases over which the lowest split ratio of a row is searched."""

MAX_CLIP_FRACTION = 1.14
"""Highest clip level of a pixel, as a fraction of the clip level of the detector."""

NEIGHBOUR_FRACTION = 8.0 / 9.0
"""Fraction of its row's threshold (above the bias) from which a pixel touching a saturated pixel is
saturated too."""

ISOLATION_FRACTION = 0.3
"""A pixel above the threshold whose four neighbours all stay below this fraction of its own value
(above the bias) is a hot pixel or a particle hit, not part of a spot."""

GROW_ITERATIONS = 4
"""Number of times the saturated pixels are extended to their neighbours."""

SPOT_RADIUS = 5.0
"""Distance (pix) from its centroid within which a saturated pixel makes a spot saturated."""

ADC_MAX = 4095
"""Largest value of one read (12 bit)."""

SUM_MAX = 65535
"""Largest value of the coadd (unsigned 16 bit)."""

BACKGROUND_STRIDE = 16
"""Stride, along both axes, of the pixels whose median gives the background of an image."""

MAX_SATURATED_PIXELS = 200000
"""Largest number of pixels above their threshold in an image of spots; an image with more is flooded."""

MASK_PLANES = dict(SAT=1)
"""Bit number of each plane of the mask written with an image."""


class SpotFlag(enum.IntFlag):
    """Bits of ``mcs_data.flags``."""

    SATURATED = 1
    """A saturated pixel lies within SPOT_RADIUS of the centroid."""

    SATURATION_UNKNOWN = 2
    """The saturation of the frame could not be measured."""


class Saturation(typing.NamedTuple):
    """Saturated pixels of an image."""

    mask: typing.Optional[np.ndarray]
    """`bool`, shape of the image; None when the saturation could not be measured."""

    reason: str
    """Why it could not be measured, else empty."""


class FloodedImage(Exception):
    """More pixels of an image exceed their threshold than a field of spots can hold."""


def nReads(exptime):
    """Number of reads ``canonexp`` sums for an exposure time.

    Parameters
    ----------
    exptime : `float`
        Exposure time (s).

    Returns
    -------
    `int`
    """
    return int(math.ceil(round(exptime * 1000) / round(LOOP_TIME * 1000))) + 1


def readExposures(phase, exptime):
    """Seconds of the exposure each read of a row collects.

    Read k = 1..n of a row ends ``phase + (k - 1) READ_PERIOD`` after the shutter opens and holds the
    light of the READ_PERIOD before.

    Parameters
    ----------
    phase : `numpy.ndarray`
        Delay (s) between the shutter opening and the end of the first read of the row.
    exptime : `float`
        Exposure time (s).

    Returns
    -------
    `numpy.ndarray`, ``phase.shape + (n,)``
    """
    phase = np.asarray(phase, dtype=float)[..., None]
    start = phase + (np.arange(1, nReads(exptime) + 1) - 2) * READ_PERIOD
    return np.clip(np.minimum(start + READ_PERIOD, exptime) - np.maximum(start, 0.0), 0.0, None)


def splitRatio(phase, exptime):
    """Seconds of the exposure a row captures over the seconds in its most exposed read.

    Parameters
    ----------
    phase : `numpy.ndarray`
        Delay (s) between the shutter opening and the end of the first read of the row.
    exptime : `float`
        Exposure time (s).

    Returns
    -------
    `numpy.ndarray`, ``phase.shape``
        1 where the row captures nothing.
    """
    exposures = readExposures(phase, exptime)
    longest = exposures.max(axis=-1)
    lit = longest > 0
    return np.where(lit, exposures.sum(axis=-1) / np.where(lit, longest, 1.0), 1.0)


@functools.lru_cache(maxsize=16)
def rowSplitRatio(nRows, exptime):
    """Lowest split ratio each row can have, over the phases of a frame and the bounds of the row time.

    Parameters
    ----------
    nRows : `int`
        Number of image rows.
    exptime : `float`
        Exposure time (s).

    Returns
    -------
    `numpy.ndarray`, (nRows,)
        Read-only.
    """
    rows = np.arange(nRows)
    low = PHASE_MIN + (2 * ROW_TIME - ROW_TIME_MAX) * rows
    high = PHASE_MAX + ROW_TIME_MAX * rows
    origin = low.min() - PHASE_STEP
    phases = np.arange(origin, high.max() + 2 * PHASE_STEP, PHASE_STEP)
    ratio = splitRatio(phases, exptime)
    first = np.clip(np.floor((low - origin) / PHASE_STEP).astype(int), 0, len(phases) - 1)
    last = np.clip(np.ceil((high - origin) / PHASE_STEP).astype(int), 0, len(phases) - 1)
    lowest = np.array([ratio[a:b + 1].min() for a, b in zip(first, last)])
    lowest.setflags(write=False)
    return lowest


def rowThresholds(nRows, exptime, clipLevel, thresholdFactor):
    """Coadd value from which a pixel of each row is saturated.

    Parameters
    ----------
    nRows : `int`
        Number of image rows.
    exptime : `float`
        Exposure time (s).
    clipLevel : `float`
        Clip level of one read, above the bias (ADU).
    thresholdFactor : `float`
        Factor applied to the clip level: 1 puts the threshold at the clip level, a lower value
        leaves a margin.

    Returns
    -------
    `numpy.ndarray`, (nRows,)
        Raw ADU, bias included.
    """
    return nReads(exptime) * BIAS_PER_READ + thresholdFactor * clipLevel * rowSplitRatio(nRows, float(exptime))


def unmeasurable(image, exptime, clipLevel, thresholdFactor):
    """Why the saturation of an image cannot be measured.

    Parameters
    ----------
    image : `numpy.ndarray`, (ny, nx)
        Raw image (ADU).
    exptime : `float`
        Exposure time (s).
    clipLevel : `float`
        Clip level of one read, above the bias (ADU).
    thresholdFactor : `float`
        Factor applied to the clip level, see `rowThresholds`.

    Returns
    -------
    `str`
        Empty when it can be measured: the settings are levels a read can reach, the image has the shape
        of a Canon coadd and a background that is the bias of the reads of the exposure time, the sum of
        the reads cannot exceed 16 bits and every threshold is within the range of the image type.
    """
    if not (np.isfinite(clipLevel) and clipLevel > 0 and np.isfinite(thresholdFactor) and thresholdFactor > 0):
        return f'clipLevel ({clipLevel}) and thresholdFactor ({thresholdFactor}) must be positive'
    if BIAS_PER_READ + thresholdFactor * clipLevel > ADC_MAX:
        return f'a read cannot reach {thresholdFactor:g} x {clipLevel:g} ADU above its bias'

    if image.shape != CANON_SHAPE:
        return f'an image of shape {image.shape} is not a Canon coadd'
    if not (np.isfinite(exptime) and exptime > 0):
        return f'an exposure time of {exptime} s'

    n = nReads(exptime)
    background = float(np.median(image[::BACKGROUND_STRIDE, ::BACKGROUND_STRIDE]))
    reads = background / BIAS_PER_READ
    if abs(reads - n) > min(0.15 + 0.004 * n, 0.35):
        return (f'a {exptime:g} s exposure sums {n} reads but the background ({background:.0f} ADU) '
                f'is the bias of {reads:.1f}')

    if n * min(MAX_CLIP_FRACTION * clipLevel + BIAS_PER_READ, ADC_MAX) > SUM_MAX:
        return f'the 16-bit sum of the {n} reads of a {exptime:g} s exposure can wrap'

    _, unreachable = _inImageType(rowThresholds(image.shape[0], exptime, clipLevel, thresholdFactor), image.dtype)
    if unreachable.any():
        return f'the threshold exceeds the range of the image type {image.dtype}'

    return ''


def _inImageType(levels, dtype):
    """Levels a pixel of an integer type must reach, and the levels none can reach.

    Returns
    -------
    (`numpy.ndarray`, `numpy.ndarray` of `bool`)
        The levels rounded up in ``dtype`` (unchanged for a floating type), and which exceed its range.
    """
    if not np.issubdtype(dtype, np.integer):
        return levels, np.zeros(levels.shape, dtype=bool)

    top = np.iinfo(dtype).max
    return np.ceil(np.clip(levels, 0, top)).astype(dtype), levels > top


def _hasBrightNeighbour(image, rows, cols, offset):
    """Whether a 4-neighbour of each pixel, inside the image, reaches ISOLATION_FRACTION of the pixel's
    value above offset."""
    ny, nx = image.shape
    own = image[rows, cols].astype(np.float64) - offset
    brightest = np.full(own.shape, -np.inf)
    for dj, di in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        r, c = rows + dj, cols + di
        inside = (r >= 0) & (r < ny) & (c >= 0) & (c < nx)
        value = image[np.clip(r, 0, ny - 1), np.clip(c, 0, nx - 1)].astype(np.float64) - offset
        brightest = np.maximum(brightest, np.where(inside, value, -np.inf))
    return brightest >= ISOLATION_FRACTION * own


def _grow(image, mask, rows, cols, levels):
    """Add to mask, in place, the pixels at or above their row's level that touch a pixel of the mask.

    Parameters
    ----------
    rows, cols : `numpy.ndarray` of `int`
        The pixels of the mask.
    levels : `numpy.ndarray`, (ny,)
        Level of each row, in the type of the image.
    """
    ny, nx = image.shape
    offsets = [(dj, di) for dj in (-1, 0, 1) for di in (-1, 0, 1) if dj or di]
    for _ in range(GROW_ITERATIONS):
        if rows.size == 0:
            break
        r = np.concatenate([rows + dj for dj, _ in offsets])
        c = np.concatenate([cols + di for _, di in offsets])
        inside = (r >= 0) & (r < ny) & (c >= 0) & (c < nx)
        r, c = r[inside], c[inside]
        new = (image[r, c] >= levels[r]) & ~mask[r, c]
        if not new.any():
            break
        flat = np.unique(r[new] * nx + c[new])
        rows, cols = flat // nx, flat % nx
        mask[rows, cols] = True


def saturatedPixels(image, exptime, clipLevel, thresholdFactor):
    """Pixels of a Canon coadd that hold a clipped read.

    A pixel at or above the threshold of its row is saturated unless it is isolated (see
    ISOLATION_FRACTION).  So is a pixel reaching NEIGHBOUR_FRACTION of the threshold, above the bias,
    that touches a saturated pixel; this is applied GROW_ITERATIONS times.

    Parameters
    ----------
    image : `numpy.ndarray`, (ny, nx)
        Raw coadd (ADU).  Not modified.
    exptime : `float`
        Exposure time (s).
    clipLevel : `float`
        Clip level of one read, above the bias (ADU).
    thresholdFactor : `float`
        Factor applied to the clip level, see `rowThresholds`.

    Returns
    -------
    `numpy.ndarray` of `bool`, (ny, nx)

    Raises
    ------
    FloodedImage
        More than MAX_SATURATED_PIXELS pixels are at or above their threshold.
    """
    image = np.asarray(image)
    offset = nReads(exptime) * BIAS_PER_READ
    thresholds = rowThresholds(image.shape[0], exptime, clipLevel, thresholdFactor)
    levels, unreachable = _inImageType(thresholds, image.dtype)

    mask = image >= levels[:, None]
    mask[unreachable] = False
    nAbove = int(np.count_nonzero(mask))
    if nAbove > MAX_SATURATED_PIXELS:
        raise FloodedImage(f'{nAbove} pixels are above their threshold')
    if nAbove == 0:
        return mask

    rowsWith = np.flatnonzero(mask.any(axis=1))
    sub, cols = np.nonzero(mask[rowsWith])
    rows = rowsWith[sub]
    keep = _hasBrightNeighbour(image, rows, cols, offset)
    mask[rows[~keep], cols[~keep]] = False
    rows, cols = rows[keep], cols[keep]

    neighbourLevels, _ = _inImageType(offset + NEIGHBOUR_FRACTION * (thresholds - offset), image.dtype)
    _grow(image, mask, rows, cols, neighbourLevels)
    return mask


def measure(image, exptime, clipLevel, thresholdFactor):
    """Saturated pixels of an image, or why they cannot be told.

    Parameters
    ----------
    image : `numpy.ndarray`, (ny, nx)
        Raw image (ADU).  Not modified.
    exptime : `float`
        Exposure time (s).
    clipLevel : `float`
        Clip level of one read, above the bias (ADU).
    thresholdFactor : `float`
        Factor applied to the clip level, see `rowThresholds`.

    Returns
    -------
    `Saturation`
        Without a mask when `unmeasurable` gives a reason or the image is flooded.
    """
    image = np.asarray(image)
    reason = unmeasurable(image, exptime, clipLevel, thresholdFactor)
    if reason:
        return Saturation(None, reason)

    try:
        return Saturation(saturatedPixels(image, exptime, clipLevel, thresholdFactor), '')
    except FloodedImage as e:
        return Saturation(None, str(e))


def measureConfigured(config, image, exptime):
    """Saturated pixels of an image, with the settings of a configuration section.

    Parameters
    ----------
    config : `dict` or None
        ``clipLevel`` and ``thresholdFactor``, see `rowThresholds`.
    image : `numpy.ndarray`, (ny, nx)
        Raw image (ADU).  Not modified.
    exptime : `float`
        Exposure time (s).

    Returns
    -------
    `Saturation`
        Nothing is raised: without a configuration, without one of its settings or when the measurement
        fails, there is no mask and the reason says why.
    """
    if config is None:
        return Saturation(None, 'no saturation configuration')

    try:
        return measure(image, exptime, config['clipLevel'], config['thresholdFactor'])
    except Exception as e:
        logging.getLogger('saturation').exception('failed to measure the saturation')
        return Saturation(None, f'{type(e).__name__}: {e}')


def spotsOnMask(mask, x, y, radius=SPOT_RADIUS):
    """Spots with a pixel of the mask within radius of their centroid.

    Parameters
    ----------
    mask : `numpy.ndarray` of `bool`, (ny, nx)
    x, y : `numpy.ndarray`
        Centroids (pix); a spot at a non-finite position is not on the mask.
    radius : `float`
        Distance (pix) from the centroid.

    Returns
    -------
    `numpy.ndarray` of `bool`, one per spot
    """
    x, y = np.atleast_1d(np.asarray(x, dtype=float)), np.atleast_1d(np.asarray(y, dtype=float))
    ny, nx = mask.shape
    finite = np.isfinite(x) & np.isfinite(y)
    xc, yc = np.where(finite, x, 0.0), np.where(finite, y, 0.0)

    half = int(math.ceil(radius))
    dy, dx = np.mgrid[-half:half + 1, -half:half + 1]
    rows = np.rint(yc).astype(int)[:, None] + dy.ravel()[None, :]
    cols = np.rint(xc).astype(int)[:, None] + dx.ravel()[None, :]

    valid = (np.hypot(rows - yc[:, None], cols - xc[:, None]) <= radius) & finite[:, None]
    valid &= (rows >= 0) & (rows < ny) & (cols >= 0) & (cols < nx)
    return (mask[np.clip(rows, 0, ny - 1), np.clip(cols, 0, nx - 1)] & valid).any(axis=1)


def spotFlags(mask, x, y):
    """Saturation bits of spots.

    Parameters
    ----------
    mask : `numpy.ndarray` of `bool`, (ny, nx), or None
        Saturated pixels of the image; None when its saturation could not be measured.
    x, y : `numpy.ndarray`
        Centroids of the spots (pix).

    Returns
    -------
    `numpy.ndarray` of `int`, one per spot
        `SpotFlag.SATURATED` for a spot with a saturated pixel within SPOT_RADIUS of its centroid;
        `SpotFlag.SATURATION_UNKNOWN` for every spot without a mask.
    """
    if mask is None:
        return np.full(len(np.atleast_1d(x)), int(SpotFlag.SATURATION_UNKNOWN), dtype=np.int32)

    return np.where(spotsOnMask(mask, x, y), int(SpotFlag.SATURATED), 0).astype(np.int32)


def flagsOfSpots(flags, spotIds, wanted):
    """Flags of the spots with the wanted identifiers.

    Parameters
    ----------
    flags : `numpy.ndarray` of `int`
        `SpotFlag` bits, one per spot.
    spotIds : `numpy.ndarray`
        Identifier of each of those spots.
    wanted : `numpy.ndarray`
        Identifiers; those that are not positive stand for no spot and are left out.

    Returns
    -------
    `numpy.ndarray` of `int`, one per positive identifier of ``wanted``
        `SpotFlag.SATURATION_UNKNOWN` for an identifier that is not among ``spotIds``.
    """
    flags = np.asarray(flags, dtype=int)
    spotIds = np.asarray(spotIds).astype(int)
    wanted = np.asarray(wanted).astype(int)
    wanted = wanted[wanted > 0]
    if spotIds.size == 0:
        return np.full(wanted.shape, int(SpotFlag.SATURATION_UNKNOWN))

    order = np.argsort(spotIds, kind='stable')
    position = np.clip(np.searchsorted(spotIds[order], wanted), 0, spotIds.size - 1)
    found = spotIds[order][position] == wanted
    return np.where(found, flags[order][position], int(SpotFlag.SATURATION_UNKNOWN))


def countSaturated(flags):
    """Number of saturated spots among flags.

    Parameters
    ----------
    flags : `numpy.ndarray` of `int`
        `SpotFlag` bits, one per spot.

    Returns
    -------
    `int` or None
        None when the saturation of any of the spots is unknown.
    """
    flags = np.asarray(flags, dtype=int)
    if (flags & SpotFlag.SATURATION_UNKNOWN).any():
        return None

    return int(((flags & SpotFlag.SATURATED) != 0).sum())


def countKeyword(name, frameId, flags):
    """Keyword giving how many of a frame's spots are saturated.

    Parameters
    ----------
    name : `str`
        Name of the keyword.
    frameId : `int`
        MCS frame.
    flags : `numpy.ndarray` of `int`
        `SpotFlag` bits, one per spot.

    Returns
    -------
    `str`
        ``name=frameId,nSaturated,nSpots``, nSaturated being -1 when the saturation is unknown.
    """
    nSaturated = countSaturated(flags)
    return f'{name}={frameId},{-1 if nSaturated is None else nSaturated},{len(flags)}'


def exceedsFraction(flags, maxFraction):
    """Whether more than a fraction of the spots are saturated.

    Parameters
    ----------
    flags : `numpy.ndarray` of `int`
        `SpotFlag` bits, one per spot.
    maxFraction : `float` or None
        Largest fraction of saturated spots that is not an excess; None for no limit.

    Returns
    -------
    `bool`
        False when the saturation is unknown.
    """
    nSaturated = countSaturated(flags)
    if nSaturated is None or maxFraction is None:
        return False

    return nSaturated > float(maxFraction) * len(flags)


def maskHdu(shape, mask, reason):
    """Mask of an image, with its saturated pixels in the SAT plane.

    Parameters
    ----------
    shape : (`int`, `int`)
        Shape of the image.
    mask : `numpy.ndarray` of `bool`, ``shape``, or None
        Saturated pixels; None when the saturation could not be measured.
    reason : `str`
        Why it could not be measured, else empty.

    Returns
    -------
    `astropy.io.fits.CompImageHDU`
        MASK, unsigned 8 bit: a saturated pixel has the bit ``MASK_PLANES['SAT']`` set, and each plane's
        bit number is in a card ``MP_<plane>``.  SATKNOWN says whether the saturation was measured,
        SATWHY why it was not, any character a card cannot hold being replaced by a space; no bit is
        set when it was not.
    """
    if mask is None:
        data = np.zeros(shape, dtype=np.uint8)
    else:
        data = mask.astype(np.uint8) << MASK_PLANES['SAT']

    hdu = pyfits.CompImageHDU(data, name='MASK', compression_type='RICE_1')
    for plane, bit in MASK_PLANES.items():
        hdu.header[f'MP_{plane}'] = (bit, 'Bit number of the mask plane')
    hdu.header['SATKNOWN'] = (mask is not None, 'Saturation of the image measured')
    if reason:
        hdu.header['SATWHY'] = ''.join(c if ' ' <= c <= '~' else ' ' for c in reason)
    return hdu
