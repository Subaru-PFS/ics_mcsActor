"""Point replayed cobra_match rows at the targets fps actually commanded.

Offline there is no image to centroid, so the recorded positions of some other visit are
replayed and the convergence loop reads back something unrelated to what it asked for.
This rewrites those positions as a measurement of the current cobra_target instead, which
is enough for the loop to close and for fiberStatus to describe this visit.

The error model is a bulk that converges geometrically, a tail that does not, and a
fraction the camera never matches.  Tail and blind membership are drawn per cobra and are
stable across the iterations of a visit, so a cobra that misses keeps missing.
"""

import numpy as np

FIRST_SIGMA_MM = 0.100
"""Scatter of the first move, before any measurement feeds back."""

CONVERGENCE_RATIO = 0.35
"""Factor the scatter shrinks by each iteration."""

REPEATABILITY_MM = 0.004
"""Scatter the loop cannot beat, whatever the iteration."""

TAIL_FRACTION = 0.05
"""Cobras whose scatter never shrinks; they end beyond any sensible tolerance."""

TAIL_SIGMA_MM = 0.080
"""Scatter for a tail cobra, at every iteration."""

BLIND_FRACTION = 0.015
"""Cobras the camera never matches, reported with spot_id -1 and no position."""


def perCobraDraw(pfsVisitId, cobraId):
    """Uniform draw in [0, 1) per cobra, identical across the iterations of a visit.

    Parameters
    ----------
    pfsVisitId : `int`
    cobraId : `numpy.ndarray` of `int`

    Returns
    -------
    `numpy.ndarray` of `float`
    """
    seeds = (np.int64(pfsVisitId) << np.int64(16)) + np.asarray(cobraId, dtype=np.int64)
    return np.array([np.random.default_rng(int(s)).random() for s in seeds])


def scatterFor(iteration, isTail):
    """Positional scatter in millimetres, per cobra.

    Parameters
    ----------
    iteration : `int`
        Zero-based convergence iteration.
    isTail : `numpy.ndarray` of `bool`

    Returns
    -------
    `numpy.ndarray` of `float`
    """
    converging = np.maximum(FIRST_SIGMA_MM * CONVERGENCE_RATIO ** iteration, REPEATABILITY_MM)
    return np.where(isTail, TAIL_SIGMA_MM, converging)


def measureAgainstTargets(df, mcsFrameId, cursor):
    """Replace replayed positions with a measurement of this frame's commanded targets.

    Parameters
    ----------
    df : `pandas.DataFrame`
        Replayed cobra_match rows, already carrying this frame's ids.
    mcsFrameId : `int`
        ``pfs_visit_id * 100 + iteration``.
    cursor : `psycopg2.extensions.cursor`
        Open cursor on the database holding cobra_target.

    Returns
    -------
    `pandas.DataFrame`
        ``df`` with pfi_center_x_mm, pfi_center_y_mm and spot_id rewritten.  Returned
        unchanged when no target has been commanded for the frame, which happens on the
        first exposure of a visit: it measures where the cobras start.
    """
    pfsVisitId, iteration = int(mcsFrameId) // 100, int(mcsFrameId) % 100

    cursor.execute('SELECT cobra_id, pfi_target_x_mm, pfi_target_y_mm FROM cobra_target '
                   'WHERE pfs_visit_id = %s AND iteration = %s', (pfsVisitId, iteration))
    rows = cursor.fetchall()
    if not rows:
        return df

    target = {int(cobraId): (x, y) for cobraId, x, y in rows}
    cobraId = df.cobra_id.to_numpy().astype(int)
    known = np.array([cid in target for cid in cobraId])
    xy = np.array([target.get(cid, (np.nan, np.nan)) for cid in cobraId], dtype=float)

    # cobra_match.spot_id is a foreign key into mcs_data for the same frame, so a spot
    # the replayed centroids do not contain cannot be referenced.
    cursor.execute('SELECT spot_id FROM mcs_data WHERE mcs_frame_id = %s', (mcsFrameId,))
    spots = {int(spotId) for (spotId,) in cursor.fetchall()}
    spotId = df.spot_id.to_numpy().astype(int)
    hasSpot = np.array([sid in spots for sid in spotId])

    draw = perCobraDraw(pfsVisitId, cobraId)
    isTail = draw < TAIL_FRACTION
    isBlind = (draw >= TAIL_FRACTION) & (draw < TAIL_FRACTION + BLIND_FRACTION)
    # A cobra with no commanded target was not moved, so it cannot be measured against
    # one; report it unmatched rather than inventing a position.
    isBlind |= ~known | ~np.isfinite(xy[:, 0]) | ~hasSpot

    rng = np.random.default_rng(int(mcsFrameId))
    sigma = scatterFor(iteration, isTail)
    x = xy[:, 0] + sigma * rng.standard_normal(len(cobraId))
    y = xy[:, 1] + sigma * rng.standard_normal(len(cobraId))
    x[isBlind] = np.nan
    y[isBlind] = np.nan

    df = df.copy()
    df.loc[:, 'pfi_center_x_mm'] = x
    df.loc[:, 'pfi_center_y_mm'] = y
    df.loc[:, 'spot_id'] = np.where(isBlind, -1, spotId)

    return df
