"""Spots of an MCS frame from the spotFinder library, as the actor keeps them."""
import numpy as np

import mcsActor.spotFinder.spotFinder as spotFinder

GAIN = 2.24
"""e-/ADU of the Canon camera."""

SPOT_COLUMNS = ('spotId', 'x', 'y', 'x2', 'y2', 'xy', 'background', 'peak', 'flux', 'fluxerr')
"""Columns of the array of spots, in order."""


def measureSpots(image, centParms):
    """Spots of a frame, one row each, with the columns `SPOT_COLUMNS`.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
        Frame as read from the camera; left untouched.
    centParms : `dict`
        findSigma, nmin, boxFind, boxCent, fwhmx, fwhmy, maxIt, aperture: see `spotFinder.findSpots`.

    Returns
    -------
    points : `numpy.ndarray`, shape (n, 10)
        spotId counts from 1. x, y, x2, y2, xy are in pixels; background is the level subtracted from the
        column of the spot; peak, flux and fluxerr are ADU above the column levels.
    levels : `dict`
        background, noise, threshold of the frame (ADU; the threshold is above the background).

    Raises
    ------
    ValueError, RuntimeError
        See `spotFinder.findSpots`.
    """
    spots, levels = spotFinder.findSpots(image, findSigma=centParms['findSigma'], nmin=centParms['nmin'],
                                         boxFind=centParms['boxFind'], boxCent=centParms['boxCent'],
                                         fwhmx=centParms['fwhmx'], fwhmy=centParms['fwhmy'], maxIt=centParms['maxIt'],
                                         aperture=centParms['aperture'], gain=GAIN)

    points = np.column_stack([np.arange(len(spots)) + 1, spots['x'], spots['y'], spots['x2'], spots['y2'], spots['xy'],
                              spots['background'], spots['peak'], spots['flux'], spots['fluxerr']])
    return points, levels
