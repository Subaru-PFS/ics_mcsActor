#cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stdlib cimport free

cdef extern from "spotFinder.h" nogil:
    int SPOTS_BAD_IMAGE
    int SPOTS_FLOODED
    int SPOTS_TOO_MANY
    int SPOTS_NO_MEMORY

    struct spotParams:
        double findSigma
        int nmin
        int boxFind
        int boxCent
        double fwhmx
        double fwhmy
        int maxIt
        double aperture
        double gain
        int nThreads

    struct spot:
        double x
        double y
        double x2
        double y2
        double xy
        double peak
        double flux
        double fluxerr
        double background

    struct frameLevels:
        double background
        double noise
        double threshold

    int c_columnLevels "columnLevels"(const int *image, int n_x, int n_y, int nThreads, double *levels, double *noise)
    double c_apertureSum "apertureSum"(const int *image, int n_x, int n_y, double x, double y, double radius, int *nPixels)
    spot *c_findSpots "findSpots"(int *image, int n_x, int n_y, const spotParams *params, frameLevels *levels, int *nSpots)


SPOT_DTYPE = [(name, '<f8') for name in ('x', 'y', 'x2', 'y2', 'xy', 'peak', 'flux', 'fluxerr', 'background')]
"""Fields of a spot, as spotFinder.h declares them."""

VALUE_LIMIT = 2 ** 30
"""Largest pixel value, of either sign, the library reduces."""


def integerFrame(image):
    """The frame as an array, checked to hold integers within `VALUE_LIMIT`.

    Raises
    ------
    TypeError
        The frame is not of an integer type.
    ValueError
        A pixel is beyond `VALUE_LIMIT`.
    """
    frame = np.asarray(image)
    if frame.dtype.kind not in 'iu':
        raise TypeError(f'a frame of integers is needed, not of {frame.dtype}')
    if frame.dtype.itemsize > 2 and frame.size and (int(frame.min()) < -VALUE_LIMIT or int(frame.max()) > VALUE_LIMIT):
        raise ValueError(f'a frame with pixels beyond {VALUE_LIMIT} cannot be reduced')
    return frame


def findSpots(image, *, double findSigma, int nmin, int boxFind, int boxCent, double fwhmx, double fwhmy, int maxIt,
              double aperture, double gain, int nThreads=8):
    """Spots of a frame, found and measured above the level of each of its columns.

    The level of a column is its median. The noise of the frame is the root mean square, over the columns,
    of the median minus the 15.87th percentile of the column: the standard deviation of the pixels,
    whatever light is on them. A spot is a set of connected pixels above findSigma times the noise.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
        Frame, of any integer type, with pixels within 2^30 of zero; left untouched.
    findSigma : `float`
        Detection threshold, in units of the noise. The threshold is truncated to an integer, 1 ADU at least;
        a spot has one pixel above it and its others at or above it.
    nmin : `int`
        Pixels at or above the threshold a spot holds at least.
    boxFind : `int`
        Spots whose centroids are closer than boxFind + 1 pixels are one, the one with the highest peak.
        No spot is kept closer to the edges of the frame than three times the larger of boxFind and boxCent.
    boxCent : `int`
        Radius (pix) of the window of the centroid.
    fwhmx, fwhmy : `float`
        Full widths at half maximum (pix) of the Gaussian weights of the window, for x and for y.
    maxIt : `int`
        Iterations of the windowed centroid at most; it stops when both axes move by less than 1e-4 pixel.
    aperture : `float`
        Radius (pix) the flux is summed in; no flux when it is not positive.
    gain : `float`
        e-/ADU.
    nThreads : `int`
        Threads used at most. The result does not depend on it.

    Returns
    -------
    spots : `numpy.ndarray`
        Record array, one row per spot, in the order their first pixel comes in the frame: x, y (windowed
        centroid, pix, a pixel being centred on its indices), x2, y2, xy (second moments of the pixels above
        the threshold, pix^2), peak (highest pixel, ADU above the level of its column), flux (sum within the
        aperture, ADU above the column levels), fluxerr (its statistical error: the noise of the pixels
        summed and the photon noise of the flux) and background (level subtracted from the column of the
        centroid, ADU: the median of the column, rounded to an integer).
    levels : `dict`
        background (median of the column levels), noise, threshold (as applied, above the column levels); ADU.

    Raises
    ------
    TypeError, ValueError
        The frame is not of integers within 2^30 of zero, has no pixel, or has more than 65535 rows.
    RuntimeError
        The frame holds too many pixels above the threshold, or too many spots, to be a frame of spots.
    MemoryError
        The memory the frame needs could not be had.
    """
    cdef np.ndarray[int, ndim=2, mode="c"] work = np.array(integerFrame(image), dtype=np.intc, order='C')
    cdef spotParams params
    cdef frameLevels levels
    cdef spot *found
    cdef int nSpots[1]
    cdef char *buffer

    params.findSigma = findSigma
    params.nmin = nmin
    params.boxFind = boxFind
    params.boxCent = boxCent
    params.fwhmx = fwhmx
    params.fwhmy = fwhmy
    params.maxIt = maxIt
    params.aperture = aperture
    params.gain = gain
    params.nThreads = nThreads

    with nogil:
        found = c_findSpots(<int *>work.data, work.shape[1], work.shape[0], &params, &levels, nSpots)

    if nSpots[0] == SPOTS_BAD_IMAGE:
        raise ValueError(f'a frame of {work.shape[0]} rows by {work.shape[1]} columns cannot be reduced: '
                         f'1 to 65535 rows, 1 column at least')
    if nSpots[0] == SPOTS_FLOODED:
        raise RuntimeError(f'too many pixels above the threshold of {levels.threshold:.1f} ADU for a frame of spots')
    if nSpots[0] == SPOTS_TOO_MANY:
        raise RuntimeError('too many spots for a frame of spots')
    if nSpots[0] == SPOTS_NO_MEMORY:
        raise MemoryError('no memory to find the spots of the frame')

    buffer = <char *>found
    spots = np.frombuffer(buffer[:sizeof(spot) * nSpots[0]], dtype=SPOT_DTYPE).copy()
    free(found)

    return spots, dict(background=levels.background, noise=levels.noise, threshold=levels.threshold)


def columnLevels(image, *, int nThreads=8):
    """Level of each column of a frame, and its noise.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
        Frame, of any integer type, with pixels within 2^30 of zero.
    nThreads : `int`
        Threads used at most.

    Returns
    -------
    levels : `numpy.ndarray`, shape (nx,)
        Median of each column, as numpy.median(image, axis=0).
    noise : `float`
        Root mean square, over the columns, of the median minus the 15.87th percentile of the column.

    Raises
    ------
    TypeError, ValueError
        The frame is not of integers within 2^30 of zero, has no pixel, or has more than 65535 rows.
    MemoryError
        The memory the frame needs could not be had.
    """
    cdef np.ndarray[int, ndim=2, mode="c"] work = np.ascontiguousarray(integerFrame(image), dtype=np.intc)
    cdef np.ndarray[double, ndim=1, mode="c"] levels = np.empty(work.shape[1], dtype=np.float64)
    cdef double noise[1]
    cdef int status

    with nogil:
        status = c_columnLevels(<int *>work.data, work.shape[1], work.shape[0], nThreads, <double *>levels.data, noise)

    if status == -1:
        raise ValueError(f'a frame of {work.shape[0]} rows by {work.shape[1]} columns cannot be reduced: '
                         f'1 to 65535 rows, 1 column at least')
    if status < 0:
        raise MemoryError('no memory to reduce the columns of the frame')

    return levels, noise[0]


def apertureFlux(image, x, y, double radius):
    """Sum of a background-free image within radius pixels of each position.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
        Of any integer type, with pixels within 2^30 of zero.
    x, y : `numpy.ndarray`, shape (n,)
        Positions (pix), a pixel being centred on its indices.
    radius : `float`
        Pixels; a pixel is summed when its centre is within it.

    Returns
    -------
    `numpy.ndarray`, shape (n,)
        NaN where the position is not on the image or the aperture leaves it.
    """
    cdef np.ndarray[int, ndim=2, mode="c"] work = np.ascontiguousarray(integerFrame(image), dtype=np.intc)
    cdef np.ndarray[double, ndim=1, mode="c"] xPos = np.ascontiguousarray(x, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] yPos = np.ascontiguousarray(y, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] flux = np.empty(xPos.shape[0], dtype=np.float64)
    cdef int nPixels[1]
    cdef int k

    if yPos.shape[0] != xPos.shape[0]:
        raise ValueError(f'{xPos.shape[0]} x positions for {yPos.shape[0]} y positions')

    for k in range(xPos.shape[0]):
        flux[k] = c_apertureSum(<int *>work.data, work.shape[1], work.shape[0], xPos[k], yPos[k], radius, nPixels)

    return flux
