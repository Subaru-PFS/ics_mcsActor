"""Spots of an MCS frame from the spotFinder library, as the actor keeps them."""
import numpy as np

import mcsActor.spotFinder.spotFinder as spotFinder

SPOT_COLUMNS = ('spotId', 'x', 'y', 'x2', 'y2', 'xy', 'background', 'peak', 'flux', 'fluxerr')
"""Columns of the array of spots, in order."""


def measureSpots(image, centParms, gain):
    """Spots of a frame, one row each, with the columns `SPOT_COLUMNS`.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
        Frame as read from the camera; left untouched.
    centParms : `dict`
        findSigma, nmin, boxFind, boxCent, fwhmx, fwhmy, maxIt, aperture: see `spotFinder.findSpots`.
    gain : `float`
        e-/ADU of the detector, for the photon noise of the flux errors.

    Returns
    -------
    points : `numpy.ndarray`, shape (n, 10)
        spotId counts from 1. x, y, x2, y2, xy are in pixels; background is the level subtracted at the
        spot; peak, flux and fluxerr are ADU above the background.
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
                                         aperture=centParms['aperture'], gain=gain)

    points = np.column_stack([np.arange(len(spots)) + 1, spots['x'], spots['y'], spots['x2'], spots['y2'], spots['xy'],
                              spots['background'], spots['peak'], spots['flux'], spots['fluxerr']])
    return points, levels
