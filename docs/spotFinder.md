# spotFinder: what the measured spots mean

`cMethod: spotFinder` (pfs_instdata `config/actors/mcs.yaml`, or `mcs switchCMethod cMethod=...`) measures
the spots of a frame with the spotFinder library (`python/mcsActor/spotFinder`); `cMethod: win` uses the
windowed centroid of `python/mcsActor/windowedCentroid`. `mcs_exposure.measurement_algorithm` records
which one measured a frame, as the literal strings `'spotFinder'` and `'windowed'`, quotes included.

## The steps

1. A smooth background is subtracted from the frame (a copy: the frame written to the PFSC file is the
   raw one): the medians of blocks of about 256 pixels a side (every fourth row of each), interpolated
   bilinearly between block centres, by their deviation from the median of all blocks, rounded to an
   ADU. A frame whose background is flat to half an ADU is left as it is; ramps and blobs of scattered
   light are removed.
2. The median of every column of what is left is subtracted: the column levels and the pattern of every
   eighth column.
3. The noise is the root mean square, over the columns, of the median minus the 15.87th percentile of
   what is left of each column: the noise of the pixels, plus whatever structure the smooth background
   does not follow (the front-lit cobra hardware of a frame taken with the lights on, a blob of
   scattered light narrower than a block), which raises the threshold over it. The detection threshold is
   `findSigma` times the noise, 1 ADU at least.
4. A spot is a 4-connected region of at least `nmin` pixels above the threshold; regions closer than
   `boxFind` + 1 pixels are merged, the brighter peak kept.
5. Its position is the Gaussian-windowed centroid (`fwhmx`, `fwhmy`, within `boxCent`), iterated until
   both axes move by less than 1e-4 pixel or `maxIt` iterations.
6. Its flux is the sum within `aperture` pixels of the centroid, above the background.

## The columns of mcs_data

Positions agree with the windowed centroid to 0.2 µm rms for spots above 3000 ADU, and are about 2
times less noisy for faint ones. Measuring a frame takes about 0.15 s on 8 threads. The other columns change
meaning:

| column | `spotFinder` | `windowed` |
|---|---|---|
| `mcs_center_x_pix`, `mcs_center_y_pix` | windowed centroid, converged on both axes | windowed centroid, stopped when either axis converges |
| `mcs_second_moment_*_pix` | second moments of the pixels above the threshold, background subtracted | moments including the background; 30-50% larger |
| `peakvalue` | highest pixel, ADU above the background (about 1795 ADU lower than a raw value at 4.8 s) | raw ADU |
| `bgvalue` | background subtracted at the centroid: level of its column and deviation of the smooth background | one value for the frame |
| `flux` | sum within `aperture` pixels, above the background | NaN unless `doPhot` |
| `fluxerr` | statistical error (noise of the pixels summed, photon noise with the detector gain), plus the saturation term of `docs/saturation.md` | NaN unless `doPhot` |
| `flags` | `MCSSpotFlags`, see `docs/saturation.md` | `MCSSpotFlags`, see `docs/saturation.md` |

Plots or analyses of `peakvalue`, `bgvalue` or the moments across the switch must select on
`measurement_algorithm`.

## Known limits

- **Compact scattered light**: a blob about one block wide (sigma 300 pixels) is followed only roughly by
  the bilinear background, and raises the noise; at 1000 ADU it moves spots on its flanks by up to 1.4
  µm, at 3000 ADU it loses 1% of them. Broader light, and ramps of any direction, are removed: on a frame
  with 2475 spots, ramps of 300 and 1000 ADU and blobs of 400 to 1000 ADU (sigma 600 and 1500 pixels)
  leave every spot, within 0.1 µm rms, with at most one spurious detection, and the noise within 7%.
- **Close pairs**: two spots closer than about 20 pixels along y, or 15 along x, merge into one, and the
  survivor is biased by up to 0.5 pixel; the windowed centroid loses the same pairs.
