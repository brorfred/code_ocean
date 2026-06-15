# abstemp.seagrid

Utilities for downloading and opening ocean model grids.

## Modules

### `cmip6`

Downloads and opens CMIP6 SST files from the Copernicus Climate Data Store.

```python
from abstemp.seagrid import cmip6

# Open a single model (downloaded files must be in src/abstemp/data/cmip6/)
ds = cmip6.open_dataset(model="cesm2", experiment="ssp5_8_5")

# Download all models for both SSP experiments (requires ~/.cdsapirc)
cmip6.retrieve_all_files()

# Available models
print(list(cmip6.cmip6_sst_models.keys()))
# ['access_cm2', 'cams_csm1_0', 'canesm5_canoe', 'cesm2', ...]

# Available experiments
print(cmip6.experiments)  # ['ssp2_4_5', 'ssp5_8_5']
```

`open_dataset` normalises all models to a common schema:

| Coordinate/Variable | Description |
|---------------------|-------------|
| `sst` | Sea-surface temperature (°C) |
| `lat`, `lon` | 1-D coordinates for regular grids |
| `lats`, `lons` | 2-D coordinates for curvilinear grids |

Longitudes are centred on 0° (−180 to 180) by `center_on_gmt`.

### `glorys`

Downloads and opens daily GLORYS12v1 ocean reanalysis fields from the Mercator
THREDDS server via OPeNDAP.

```python
from abstemp.seagrid import glorys

ds = glorys.open_dataset(dtm="2020-06-15", data_var="temp")
da = glorys.open_dataarray(dtm="2020-06-15", data_var="uvel")
```

Available `data_var` keys: `uvel`, `vvel`, `wvel`, `temp`, `salt`.

> `open_dataarray` had a `NameError` on the return statement (`data_array`
> instead of `data_var`) — fixed.

## Notes

- `seagrid` has no `__init__.py` (implicit namespace package).  See ISSUES.md S1.
- `_mercator.py` is an incomplete stub with missing imports and should not be
  used.  See ISSUES.md O1.
