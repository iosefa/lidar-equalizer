# lidar-equalizer
*A lightweight C++17 tool and PDAL extension for overlap-aware LiDAR density normalization.*

## Overview
`lidar-equalizer` balances point densities across overlapping LiDAR flight lines using per-cell proportional thinning.
 It operates directly on LAS/LAZ files via PDAL, producing uniform, analysis-ready datasets for canopy metrics or terrain modeling.

## Build
```bash
brew install pdal cmake
git clone https://github.com/<your-username>/lidar-equalizer.git
cd lidar-equalizer
mkdir build && cd build
cmake ..
make -j
```

## Usage
```
./lidar_equalizer input.laz output_equalized.laz [cell_size] [target_density] [seed] \
    [--class-scope=all|nonground|ground] \
    [--flag-overlap-only|--equalize-overlap-only|--join-overlap-by-scan-angle]
```

### Example
```
./lidar_equalizer tile_2018.laz tile_2018_eq.laz 1.0 6.0 42 --class-scope=nonground
```

`--class-scope` defaults to `all`. Use `nonground` to thin only non-ground points (Classification ≠ 2) or `ground` to thin only the ground returns while passing the rest through unchanged.

### Overlap flagging only
```
./lidar_equalizer tile_2018.laz tile_2018_overlap.laz --flag-overlap-only --min-points-per-psid=5 --overlap-dilate=1 --swath-key=psid-channel
```

### Overlap-aware modes
- `--flag-overlap-only` – run the standalone overlap flagger (read → flag → write) without thinning.
- `--equalize-overlap-only` – thin only the points inside the detected overlap mask; non-overlap points are forwarded unchanged.
- `--join-overlap-by-scan-angle` – keep the point with the smallest absolute scan angle within each overlap cell (ties broken by return number, intensity, GPS time).

Shared options:
- `--min-points-per-psid=<n>` – require at least `n` points from a swath in a cell before it contributes to overlap (default `1`).
- `--overlap-dilate=<r>` – dilate the overlap mask by `r` grid cells to smooth small gaps (default `0`).
- `--overwrite-overlap` – clear existing overlap bits instead of merging with the input flags.
- `--swath-key=psid|psid-channel|psid-file|channel` – choose how swaths are identified when balancing quotas and flagging overlap (default `psid`).
- `--cell-size=<meters>` – override the positional cell size argument with a named flag.
- `--self-test-overlap` – developer utility that runs the synthetic overlap regression checks without touching disk.

## How it works
1. Grids the LiDAR domain (default 1 m²).
2. Counts total and per-PSID point density.
3. Allocates proportional quotas for each swath in overlap cells.
4. Randomly subsamples to achieve target uniform density.

## Why?
Merging multiple LiDAR flight lines often doubles point density in overlap zones. This tool ensures balanced coverage—avoiding artefacts in canopy metrics or structural models.

## License
MIT © 2025 Iosefa Percival
