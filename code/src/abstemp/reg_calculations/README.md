# abstemp.reg_calculations

Spatial aggregation of gridded SST fields onto the ~11 000 mintmat ocean regions
and Longhurst province mapping.

## Overview

The mintmat dataset divides the global ocean into ~11 000 regions (~0.5° resolution).
This sub-package maps both OSTIA (full-resolution 0.05°) and CMIP6 (model-dependent)
SST fields onto those regions and computes per-region SST statistics.

## Modules

### Package `__init__`

High-level helpers for the 2°×2° intermediate grid:

```python
from abstemp.reg_calculations import (
    nearest,        # map mintmat region centroids → 2° grid indices
    nreg_arr,       # 2D array: 2°-grid cell → region index
    arr_to_regvec,  # 2D spatial array → 1D per-region vector
    regvec_to_arr,  # 1D per-region DataArray → 2D spatial array
    latlon_to_reg,  # (lat, lon) → region index
    reg_area,       # region areas (km²) on the 2° grid
    add_maxmon,     # write SST statistics into mintmat_2001-2009.nc
    add_longhurst,  # write Longhurst province codes into mintmat_2001-2009.nc
)
```

### `sst_ostia`

Maps OSTIA SST fields onto mintmat regions using a haversine BallTree.

```python
from abstemp.reg_calculations import sst_ostia

ij  = sst_ostia.nearest_ostia()              # pixel → region mapping
reg = sst_ostia.to_reg(arr, ij)              # average pixels into regions
reg = sst_ostia.calc_reg_maxmon(2001, 2009)  # warmest-month stats per region
```

Key functions:

| Function | Description |
|----------|-------------|
| `nearest_ostia` | BallTree lookup: OSTIA pixel → nearest region centroid |
| `to_reg` | Pixel-average SST into regions |
| `calc_reg_maxmon` | Per-year maximum month + magnitude per region |
| `regrid` | Bucket-average from 0.05° to 2° resolution |
| `max_month` / `min_month` / `clim_max_month` | Per-year max/min/climatological monthly SST arrays |

### `sst_cmip6`

Maps CMIP6 SST fields (any model) onto mintmat regions.

```python
from abstemp.reg_calculations import sst_cmip6

ds  = cmip6.open_dataset("ec_earth3_cc")
reg = sst_cmip6.to_reg(ds)   # max SST per region across all time steps
```

### `longhurst`

Maps Longhurst province codes (basins, regions, biomes) onto mintmat regions.

```python
from abstemp.reg_calculations import longhurst

ij  = longhurst.nearest_longhurst()  # Longhurst pixel → region centroid
reg = longhurst.lh_to_reg(lh.basins) # basin code per region
```

> **Note:** `lh = abstemp.open_longhurst()` runs at import time as a
> module-level side effect (see ISSUES.md S3).

## Known issues

- `calc_reg_maxmon` in `sst_ostia` and `regrid` in `sst_ostia` each contain
  unreachable blocks of dead code after an early `return`.
  See ISSUES.md D2 / D3.
- The intent of `svec.reindex(range(svec.index.max() + 1))` in
  `sst_cmip6.to_reg()` is undocumented.  See ISSUES.md L8.
