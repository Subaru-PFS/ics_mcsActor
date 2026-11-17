#!/usr/bin/env python
"""Measure the spots of raw MCS frames with spotFinder and with the windowed centroid, and compare.

    python compareWithWindowed.py /data/raw/2026-11-16/mcs/PFSC15503500.fits [more frames] [--threads 8]

For each frame: the levels each method works with, the time it takes, the spots each one finds, and
the difference of their positions by brightness of the spot. Both extensions must be built
(python setup.py build_ext -i in spotFinder and in windowedCentroid).
"""
import argparse
import contextlib
import ctypes
import os
import time

import fitsio
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import sigmaclip

import mcsActor.spotFinder.spotFinder as spotFinder
import mcsActor.windowedCentroid.centroid as centroid

PEAK_BANDS = (0, 1000, 3000, 10000, 30000, np.inf)
"""Limits (ADU above the background) of the bands of peak the differences are given in."""

WINDOWED_FIELDS = ('x', 'y', 'x2', 'y2', 'peak', 'xy', 'back')
"""Fields of a spot of the windowed centroid, in the order of its C structure."""


@contextlib.contextmanager
def quietStdout():
    """Send what compiled code prints on the standard output to /dev/null."""
    libc = ctypes.CDLL(None)
    print(end='', flush=True)
    libc.fflush(None)
    saved = os.dup(1)
    with open(os.devnull, 'w') as sink:
        os.dup2(sink.fileno(), 1)
        try:
            yield
        finally:
            libc.fflush(None)
            os.dup2(saved, 1)
            os.close(saved)


def readFrame(path):
    """Image of a raw MCS file and its exposure time (s), NaN when the header has none."""
    with fitsio.FITS(path) as hdus:
        header = hdus[0].read_header()
        image = hdus[1].read() if len(hdus) > 1 and hdus[1].has_data() else hdus[0].read()
    return image, float(header['EXPTIME']) if 'EXPTIME' in header else np.nan


def timed(function, repeat):
    """Result of function() and the shortest and median time (s) of repeat calls."""
    seconds = []
    for _ in range(repeat):
        start = time.perf_counter()
        result = function()
        seconds.append(time.perf_counter() - start)
    return result, min(seconds), float(np.median(seconds))


def windowedThresholds(image, boresight, threshSigma, findSigma):
    """Threshold and background of the windowed centroid: clipped mean and deviation within 500 pix of the boresight.

    Parameters
    ----------
    image : `numpy.ndarray`, shape (ny, nx)
    boresight : `tuple` of `float`
        x, y (pix).
    threshSigma : `float`
        Clipping, in standard deviations.
    findSigma : `float`
        Threshold above the mean, in standard deviations.

    Returns
    -------
    threshold, background : `float`
        ADU, as in the image.
    """
    x0, y0 = int(round(boresight[0])), int(round(boresight[1]))
    cut = image[y0 - 500:y0 + 501, x0 - 500:x0 + 501]
    rows, cols = np.mgrid[y0 - 500:y0 + 501, x0 - 500:x0 + 501]
    disc = np.hypot(cols - boresight[0], rows - boresight[1]) < 500
    kept, _, _ = sigmaclip(cut[disc], threshSigma, threshSigma)
    return kept.mean() + findSigma * kept.std(), kept.mean()


def windowedSpots(image, threshold, params):
    """Spots of the windowed centroid, as a record array with `WINDOWED_FIELDS`."""
    with quietStdout():
        raw = centroid.centroid_only(image.astype('<i4'), params['fwhmx'], params['fwhmy'], int(threshold),
                                     int(threshold), params['boxFind'], params['boxCent'], params['nmin'],
                                     params['maxIt'], 0)
    values = np.frombuffer(raw, dtype='<f8').reshape(-1, len(WINDOWED_FIELDS))
    spots = np.zeros(len(values), dtype=[(name, '<f8') for name in WINDOWED_FIELDS])
    for k, name in enumerate(WINDOWED_FIELDS):
        spots[name] = values[:, k]
    return spots


def matchSpots(a, b, radius):
    """Indices of the spots of a and of b that are each other's nearest, closer than radius (pix)."""
    if not len(a) or not len(b):
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    distance, nearest = cKDTree(np.column_stack([b['x'], b['y']])).query(np.column_stack([a['x'], a['y']]))
    _, back = cKDTree(np.column_stack([a['x'], a['y']])).query(np.column_stack([b['x'], b['y']]))
    mutual = (distance < radius) & (back[nearest] == np.arange(len(a)))
    return np.flatnonzero(mutual), nearest[mutual]


def listSpots(spots, peakOffset, limit=8):
    """A few spots as text: x, y and peak above the background."""
    shown = [f'({s["x"]:.1f}, {s["y"]:.1f}; {s["peak"] - peakOffset:.0f})' for s in spots[:limit]]
    return ' '.join(shown) + (f' and {len(spots) - limit} more' if len(spots) > limit else '')


def compare(path, args):
    """Measure one frame with both methods and print the comparison."""
    params = dict(findSigma=args.findSigma, nmin=args.nmin, boxFind=args.boxFind, boxCent=args.boxCent,
                  fwhmx=args.fwhmx, fwhmy=args.fwhmy, maxIt=args.maxIt, aperture=args.aperture, gain=args.gain)
    scale = args.micronsPerPixel
    image, exptime = readFrame(path)
    print(f'\n==== {path}: {image.shape[0]} x {image.shape[1]} {image.dtype}, exposure {exptime:g} s')

    (threshold, background), tThreshold, _ = timed(
        lambda: windowedThresholds(image, args.boresight, args.threshSigma, args.findSigma), args.repeat)
    old, tOld, tOldMedian = timed(lambda: windowedSpots(image, threshold, params), args.repeat)
    (new, levels), tNew, tNewMedian = timed(lambda: spotFinder.findSpots(image, nThreads=args.threads, **params),
                                            args.repeat)
    (_, _), tOne, _ = timed(lambda: spotFinder.findSpots(image, nThreads=1, **params), args.repeat)

    print(f'windowed:   background {background:8.1f}  threshold {threshold - background:6.1f} ADU above it  '
          f'{len(old):5d} spots  {1000 * tOld:5.0f} ms (median {1000 * tOldMedian:.0f}), '
          f'threshold apart: {1000 * tThreshold:.0f} ms here, once per visit in the actor')
    print(f'spotFinder: background {levels["background"]:8.1f}  threshold {levels["threshold"]:6.1f} ADU above it  '
          f'{len(new):5d} spots  {1000 * tNew:5.0f} ms (median {1000 * tNewMedian:.0f}) on {args.threads} threads, '
          f'{1000 * tOne:.0f} ms on one;  noise {levels["noise"]:.2f} ADU')

    iNew, iOld = matchSpots(new, old, args.matchRadius)
    onlyNew = np.setdiff1d(np.arange(len(new)), iNew)
    onlyOld = np.setdiff1d(np.arange(len(old)), iOld)
    print(f'in both: {len(iNew)};  only spotFinder: {len(onlyNew)}  {listSpots(new[onlyNew], 0.0)}')
    print(f'              only windowed: {len(onlyOld)}  {listSpots(old[onlyOld], background)}')
    if not len(iNew):
        return

    dx = scale * (new['x'][iNew] - old['x'][iOld])
    dy = scale * (new['y'][iNew] - old['y'][iOld])
    peak = new['peak'][iNew]
    print(f'position, spotFinder minus windowed, um at {scale:g} um per pixel, by peak above the background:')
    print(f'{"peak (ADU)":>16s} {"spots":>6s} {"median dx":>10s} {"median dy":>10s} {"rms dx":>8s} {"rms dy":>8s} '
          f'{"largest":>8s}')
    for low, high in zip(PEAK_BANDS[:-1], PEAK_BANDS[1:]):
        band = (peak >= low) & (peak < high)
        if band.any():
            label = f'{low:g} - {high:g}' if np.isfinite(high) else f'above {low:g}'
            print(f'{label:>16s} {band.sum():6d} {np.median(dx[band]):10.2f} {np.median(dy[band]):10.2f} '
                  f'{np.sqrt(np.mean(dx[band] ** 2)):8.2f} {np.sqrt(np.mean(dy[band] ** 2)):8.2f} '
                  f'{np.hypot(dx[band], dy[band]).max():8.2f}')

    flux = new['flux']
    print(f'peak: spotFinder + background minus windowed, median {np.median(peak + levels["background"] - old["peak"][iOld]):+.1f} ADU;  '
          f'moments x2, y2: spotFinder {np.median(new["x2"]):.2f}, {np.median(new["y2"]):.2f}  '
          f'windowed {np.median(old["x2"]):.2f}, {np.median(old["y2"]):.2f} pix^2')
    print(f'flux within {args.aperture:g} pix: 5th, 50th, 95th percentile {np.nanpercentile(flux, [5, 50, 95]).round(0).tolist()} ADU;  '
          f'flux / fluxerr median {np.nanmedian(flux / new["fluxerr"]):.0f};  without flux: {int(np.isnan(flux).sum())}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('paths', nargs='+', help='raw MCS files')
    parser.add_argument('--threads', type=int, default=8, help='threads of spotFinder')
    parser.add_argument('--repeat', type=int, default=5, help='calls each method is timed over')
    parser.add_argument('--matchRadius', type=float, default=2.0, help='pixels within which two spots are the same')
    parser.add_argument('--micronsPerPixel', type=float, default=86.9, help='scale the differences are given at')
    parser.add_argument('--boresight', type=float, nargs=2, default=(4438.7, 2838.3), metavar=('X', 'Y'),
                        help='centre of the region the windowed threshold is taken from, pixels')
    for name, default in (('findSigma', 25.0), ('threshSigma', 4.0), ('fwhmx', 4.0), ('fwhmy', 9.0), ('aperture', 10.0),
                          ('gain', 2.24)):
        parser.add_argument(f'--{name}', type=float, default=default, help='centroid parameter')
    for name, default in (('nmin', 10), ('boxFind', 14), ('boxCent', 14), ('maxIt', 20)):
        parser.add_argument(f'--{name}', type=int, default=default, help='centroid parameter')
    args = parser.parse_args()

    for path in args.paths:
        compare(path, args)


if __name__ == '__main__':
    main()
