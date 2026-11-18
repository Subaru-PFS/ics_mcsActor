# Saturation of MCS frames

What the mcsActor records about saturated pixels and spots of the Canon MCS camera, where it is
recorded, and how the thresholds were calibrated.

## Why pixels saturate below 65535

`canonexp` sums n = ceil(T / 0.8 s) + 1 reads into an unsigned 16-bit coadd (7 reads at 4.8 s). Each
read starts from an empty pixel, carries a bias of 256.5 ADU and clips about 2.8k ADU above it, so a
pixel saturates when its most exposed *read* clips, well below 65535. The camera free-runs with a
rolling readout, 143 µs per row, not synchronised with the shutter: how the exposure is shared between
the reads of a row depends on the row and on the frame. The coadd value at which a pixel has clipped
is therefore a function of the row, computed for the worst phase a frame can have
(`saturation.rowThresholds`). The clip level drifts, about −1.5 ADU per day in 2026 for no known
reason, and is configuration (`clipLevel`); it should be re-measured about monthly.

The longer the exposure, the larger the share of a row's light that its most exposed read holds;
exposure time does not reduce saturation, only the fiducial LED or BIA brightness does.

## Saturated pixels: the MASK HDU

Every PFSC file written by the actor has an HDU `MASK` (unsigned 8 bit, RICE compressed) after `SPOTS`.

| card | meaning |
|---|---|
| `MP_SAT` | bit number of the SAT plane: a pixel `p` is saturated when `p & (1 << MP_SAT)` |
| `SATKNOWN` | whether the saturation of the frame was measured; when False, no bit is set |
| `SATWHY` | why it was not measured (a frame that is not a Canon coadd, a background that is not the bias of the exposure time, a flooded frame, a configuration error) |
| `CLIPLVL` | clip level of one read above its bias (ADU) the frame was measured with |
| `SATFACT` | threshold as a fraction of `CLIPLVL`; 1 puts it at the clip level |

A pixel is in the plane when it reaches the threshold of its row and is not isolated (a hot pixel or a
particle hit), or reaches 8/9 of it and touches a saturated pixel (grown four times).

## Saturated spots: three flags

A spot is judged by the number of saturated pixels within `spotRadius` (10 px, the photometry
aperture) of its centroid. Three thresholds on that count set the bits of
`pfs.utils.datamodel.mcs.MCSSpotFlags`; the flags nest, each implying the ones above it.

| bit | value | flag | from (saturated pixels) | meaning | use |
|---|---|---|---|---|---|
| 0 | 1 | `NO_SAT_CHECK` | — | the saturation of the frame was not measured; the other bits are unset and say nothing | treat as unknown, not as clean |
| 1 | 2 | `BAD_PHOTOMETRY` | 2 | flux low by more than ~0.5% | do not use the flux; `fluxerr` carries the bias |
| 2 | 4 | `BAD_CENTROID` | 5 | centroid biased by more than ~1 µm | usable: the mean bias is constant at a fixed flux; the frame-to-frame scatter is 1.5 times that of an unsaturated spot (~3 µm added) |
| 3 | 8 | `SATURATED` | 10 | frame-to-frame scatter of the centroid at least twice that of an unsaturated spot | the position is not stable |

A spot below `BAD_PHOTOMETRY` (one saturated pixel) loses 0.2% of its flux in median, which is within
the per-spot noise of the calibration.

The flux error of every spot includes the saturation, in quadrature with the statistical error:

    fluxerr = hypot(statistical error, fluxErrorPerPixel * nSatPix * flux)

with `fluxErrorPerPixel` = 0.008, the flux lost per saturated pixel over one to four pixels.

### Where the flags are

- **opdb**: `mcs_data.flags`, one value per spot. The `spot_id = -1` row has 0. Rows written before the
  flags also hold 0: a frame was flagged only if its `mcs_exposure.version_actor` is the release that
  introduced them or a later one.
- **PFSC file**: the `SPOTS` table has `flags` and `nSatPix` (the count itself, -1 when unknown), so
  another threshold can be applied without re-running anything. Its header defines the bits and gives
  the thresholds the frame was flagged with:

| card | meaning |
|---|---|
| `FLAGVER` | version of the definitions below |
| `FLAG<n>` | name of the flag worth 2**n |
| `SATRAD` | radius (pix) within which saturated pixels are counted |
| `SATPHOT`, `SATCENT`, `SATSPOT` | counts from which `BAD_PHOTOMETRY`, `BAD_CENTROID` and `SATURATED` are set |
| `SATFERR` | fraction of the flux added to `fluxerr` per saturated pixel |

- **Keywords**, per frame: `saturatedSpots`, and after matching `saturatedFiducials` and
  `saturatedCobras`, each `frameId,nBadPhotometry,nBadCentroid,nSaturated,nSpots` (counts -1 when
  unknown). The actor warns when more than `maxFiducialFraction` of the matched fiducials are
  `SATURATED`.

The flags are reported only: no spot is rejected from the transform or from the cobra matching.

## Calibration

The thresholds were set on the LED duty scan of 2026-10-31 (visits 152723-152728, frames
15272300-15272804; duty 2% to 10%, five frames at each of 3.9, 5, 7.9 and 10%), with the production
rule (`clipLevel` 2815 ADU, `thresholdFactor` 1.0; 2835 ADU gives the same tables). For each of the 71
fiducials, the flux deficit is its flux against the flux it has at a duty where it is linear, scaled
by the LED level; the centroid shift is its position against the same reference, after removing the
frame's motion with the never-clipped fiducials. 1562 fiducial-frames; 86 µm per MCS pixel at the PFI.

**Flux.** Flux lost by saturated pixels within 10 px (the truth itself scatters by 0.6% per spot):

| saturated pixels | fiducial-frames | median | 84% | 95% |
|---|---|---|---|---|
| 0 | 990 | 0.0% | 0.6% | 1.3% |
| 1 | 34 | 0.2% | 1.6% | 2.6% |
| 2 | 46 | 0.75% | 1.9% | 3.1% |
| 3 | 63 | 2.0% | 3.8% | 4.8% |
| 4 | 36 | 2.6% | 4.2% | 7.0% |
| 5 | 39 | 3.6% | 8.0% | 10.3% |

**Centroid bias.** Mean shift along y, the direction of the secondary image the camera puts about 5 px
towards -y; along x it is a fifth to a half of it, with the opposite sign:

| saturated pixels | fiducial-frames | median flux lost | mean y shift (pix) |
|---|---|---|---|
| 0 | 990 | 0.0% | -0.018 (reference) |
| 1-2 | 80 | 0.5% | 0.000 |
| 3-4 | 99 | 2.3% | -0.029 |
| 5-6 | 81 | 5.3% | -0.094 |
| 7-8 | 75 | 8.4% | -0.168 |
| 9-12 | 98 | 11.6% | -0.218 |
| 13-24 | 82 | 12-16% | -0.22 to -0.35 |
| > 24 | 57 | 22-28% | -0.52 to -0.68 |

Below 5 pixels the bias is at most 0.011 px (1 µm) relative to unsaturated spots.

**Centroid stability.** Frame-to-frame rms of a fiducial's position over the five frames of one duty,
after an affine fit to the never-clipped fiducials of each frame, by the spot's mean count over the
frames (median over the spots of each class; x behaves alike from a floor of 0.039 px):

| saturated pixels | spot-duties | rms y (pix) | relative to unsaturated |
|---|---|---|---|
| 0 | 155 | 0.027 | 1 |
| 1-2 | 32 | 0.025 | 0.9 |
| 3-4 | 17 | 0.033 | 1.2 |
| 5-7 | 24 | 0.041 | 1.5 |
| 8-9 | 13 | 0.048 | 1.8 |
| 10-23 | 26 | about 0.06 | 2.2-2.35 |
| 24 and more | 12 | 0.089-0.12 | 3.3-4.5 |

The number of saturated pixels of a spot itself changes by two to four from frame to frame: the share
of the light its most exposed read holds depends on the frame. Above 10 pixels this doubles the scatter
of the centroid, about 5 µm added in quadrature.

### Limits

**The count does not weigh the same on every row.** The threshold of a row is the coadd value at which
its pixels have clipped for the worst phase a frame can have, and the light later rows collect is shared
more evenly between reads, so on high rows the mask also holds pixels that did not clip in that frame.
Flux lost per masked pixel, spots with 1 to 12 masked pixels:

| rows | flux lost per masked pixel | y shift per masked pixel (pix) |
|---|---|---|
| 0-1500 | 1.25% | -0.025 |
| 1500-3000 | 0.80% | -0.012 |
| 3000-4500 | 0.56% | -0.006 |
| 4500-5778 | 0.13% | -0.001 |

The flags are therefore conservative on high rows. On low rows a single masked pixel already costs
about 1.2% of the flux, so the 0.5% of `BAD_PHOTOMETRY` holds there only in median; unflagged spots lose
at most 1.2-1.6% (95th percentile) in every band of rows. Weighting the count by row needs the phase of
the frame, which is not recorded.

One night, one focus, fiducials only. The number of saturated pixels depends on the size of the
spots: a sharper focus puts the same flux on fewer pixels. Cobras were not in the calibration; their
spots are of a similar size. The thresholds should be checked again on the first nights of a run, with
a duty scan: per spot and frame, the count of saturated pixels within `spotRadius`, against the flux
deficit, the centroid shift and the frame-to-frame scatter, as above.

## Configuration

`saturation` section of `config/actors/mcs.yaml` in pfs_instdata:

| setting | value | |
|---|---|---|
| `clipLevel` | 2815.0 | ADU above the bias, one read; measured 2026-11-13 |
| `thresholdFactor` | 1.0 | threshold at the clip level |
| `maxFiducialFraction` | 0.5 | warn when more of the matched fiducials are `SATURATED` |
| `spotRadius` | 10 | pix |
| `badPhotometryPixels` | 2 | |
| `badCentroidPixels` | 5 | |
| `saturatedPixels` | 10 | |
| `fluxErrorPerPixel` | 0.008 | |

Without the section, or with a setting missing or out of range, every frame is `NO_SAT_CHECK` and the
reason is in `SATWHY` and in a warning.
