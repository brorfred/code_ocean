# abstemp.data

Bundled datasets and download utilities for the **abstemp** analysis.

## Bundled files

| File | Description |
|------|-------------|
| `mintmat_2001-2009.nc` | Dijkstra minimum-travel-time matrix between ~11 000 ocean regions, with per-region SST statistics (max/min/climatological) and Longhurst classifications. |
| `ostia_maxmonsst_1985-1990.nc` | Per-pixel max/min monthly SST statistics for OSTIA reanalysis 1985–1990 (variables: `maxarr`, `maxmon`, `minarr`, `minmon`, `clim_maxarr`, `clim_maxmon`, `clim_minarr`, `clim_minmon`). |
| `ostia_maxmonsst_2019-2023.nc` | Same statistics for OSTIA NRT 2019–2023. |
| `maxmonsst_cmip6/{model}_{experiment}_maxmonsst.nc` | Per-pixel max/min monthly SST statistics for CMIP6 models 2095–2100 (same variables as the OSTIA maxmonsst files).  Downloaded files include CNRM-CM6-1-HR SSP5-8.5 and SSP2-4.5 by default. |
| `Longhurst_Regions_2007.nc` | Longhurst (2007) biogeochemical province boundaries on a 0.5° grid. |
| `abstemp_reg_degvel.parquet` | Per-region degree-velocity table for 1985, 2019, and 2095. |
| `all_cmip6_hists_ssp585.csv` | Area-weighted SST histograms for all CMIP6 models, SSP5-8.5. |
| `all_cmip6_hists_ssp245.csv` | Area-weighted SST histograms for all CMIP6 models, SSP2-4.5. |
| `all_ostia_hists.csv` | Area-weighted OSTIA SST histograms for 1985–1990 and 2019–2023. |
| `growth_rates.csv` | Observed phytoplankton growth rates vs. temperature (literature compilation). |
| `growth_rates_binned_1_deg.csv` | Same data binned at 1 °C resolution. |
| `cmip6/` | Per-model CMIP6 NetCDF files downloaded via `download.retrieve_all_cmip6_files()`. |

Raw monthly SST files (`ostia_sst_1985-1990.nc`, `ostia_sst_2019-2023.nc`) are only needed to regenerate the maxmonsst files; see `generate_ostia_maxmonsst_files()` and `generate_cmip6_maxmonsst_files()`.

## Public API

```python
from abstemp import data

# Pre-computed max-month statistics (lat/lon grid, no time dimension)
ds85  = data.open_warmest_1985()   # xarray.Dataset, OSTIA reanalysis 1985–1990
ds19  = data.open_warmest_2019()   # xarray.Dataset, OSTIA NRT 2019–2023
ds95  = data.open_warmest_2095()   # xarray.Dataset, CNRM-CM6-1-HR SSP5-8.5 2095–2100 (default)
lh    = data.open_longhurst()      # xarray.Dataset, Longhurst provinces

# Raw monthly SST (time/lat/lon) — only needed to regenerate maxmonsst files
raw85 = data.open_ostia_1985()
raw19 = data.open_ostia_2019()
raw95 = data.open_cmip6_2095()    # model="ecearth" by default
```

## Downloading data

```python
from abstemp.data import download

# Download all CMIP6 SST files (requires a CDS API key in ~/.cdsapirc)
download.retrieve_all_cmip6_files()

# Download the mintmat file from Zenodo
download.mintmat()

# Download pre-computed warmest-month SST files from the project server
download.maxmonsst_fields()

# Download Longhurst province mask
download.longhurst_regions()
```

## Regenerating derived files

```python
from abstemp.data import download, all

# Download raw monthly SST files (needed as input)
download.ostia_sst_fields()

# Recompute ostia_maxmonsst_*.nc and abstemp_reg_degvel.parquet
all()
```

Or step by step:

```python
from abstemp import data

data.generate_ostia_maxmonsst_files()   # ostia_maxmonsst_1985-1990.nc, ostia_maxmonsst_2019-2023.nc
```
