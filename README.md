# Absolute Temperature — Reproducible Code Capsule

This capsule reproduces all manuscript figures for the study on how rising
sea-surface temperatures affect the biological accessibility of cooler ocean
regions via ocean connectivity.

---

## Reproducing the figures

Click **Reproducible Run** at the top of this capsule.  All six manuscript
figures will be written to the **Results** panel (usually within 10–15
minutes).  No configuration is required.

The run script (`code/run`) executes the following steps automatically:

1. Wires the output directory so figures land in `/results`.
2. Calls `abstemp.figure_scripts.generate_figures()`, which runs each figure
   generator in sequence and prints a pass/fail summary.

---

## What each figure shows

| Result file | Figure | Description |
|-------------|--------|-------------|
| `warmest_sst_map.png` | Fig. 2 | Three-panel global map of the warmest monthly SST for 1985–1990, 2019–2023, and 2095–2100 (CNRM-CM6-1-HR SSP5-8.5). Colours show SST thresholds at 25, 30–37 °C. |
| `sst_hist_1985_2019_2095.pdf` | Fig. 3 | Area-weighted histograms of peak monthly SST for three climate periods. The red envelope shows the full CMIP6 SSP5-8.5 model spread; the red line is the ensemble mean. |
| `growth_models.pdf` | Fig. 4 | Phytoplankton growth rate vs. temperature for two parametric models (Norberg–Eppley and Blackford), overlaid on the observed growth-rate compilation. |
| `regdegvel_maps.png` | Fig. 5 | Three-panel global map of regional temperature velocity (years to travel 1 °C warmer) for 1985–1990, 2019–2023, and 2095–2100. |
| `ostia_sst.png` | Fig. 6 | Three-panel map of the warmest month, coldest month, and seasonal SST range from the OSTIA 2019–2023 product. |
| `regdegvel_hist.pdf` | Fig. 7 | Area-weighted histogram of the time (years) required to travel 1–4 °C warmer, computed for 2019 SST values across ~11 000 ocean regions. |

---

## Scientific background

The analysis centres on *absolute temperature velocity*: given the current
(or projected) sea-surface temperature field, how many years does it take for
ocean currents to carry a water parcel from a given region to the nearest
region that is already ≥Δ°C warmer?  This is distinct from the classic
*climate velocity* metric (which tracks how fast isotherms move across the
surface) because it uses the connectivity encoded in ocean particle-tracking
rather than spatial SST gradients.

**Ocean regions** (~11 000 globally) are derived from the ECCO Lagrangian
particle-tracking model and stored in `mintmat_2001-2009.nc`.  The file
contains the Dijkstra-shortest-path minimum travel time between every pair of
regions (`dijkmintmat`, units: days), which is the core input to the
degree-velocity calculation.

**Temperature velocity** for each region is computed by
`abstemp.tempvel.movedegree()`: find the closest (by travel time) region
whose peak monthly SST exceeds the source region's SST by Δ°C, then divide
the travel time by 365 to convert to years.

**Three climate periods** are compared:

| Period | SST source |
|--------|-----------|
| 1985–1990 | OSTIA reanalysis (`ostia_maxmonsst_1985-1990.nc`) |
| 2019–2023 | OSTIA NRT product (`ostia_maxmonsst_2019-2023.nc`) |
| 2095–2100 | CNRM-CM6-1-HR SSP5-8.5 (`maxmonsst_cmip6/cnrm_cm6_1_hr_ssp5_8_5_maxmonsst.nc`) |

For the SST histogram comparison (Fig. 3), 24 CMIP6 models under both
SSP5-8.5 and SSP2-4.5 are included; the pre-computed histograms are stored as
CSV files in the repository.

**Phytoplankton growth models** (Fig. 4) are fully parametric and require no
input data — they implement the Norberg–Eppley and Blackford equations and
overlay the observed growth-rate compilation from `growth_rates.csv`.

---

## Data assets

All input data are in the **Data** panel of this capsule (mounted read-only at
`/data` during the run).

| File | Size | Used for |
|------|------|---------|
| `ostia_maxmonsst_1985-1990.nc` | 1.2 GB | Fig. 2, Fig. 5 |
| `ostia_maxmonsst_2019-2023.nc` | 1.2 GB | Fig. 2, Fig. 5, Fig. 6 |
| `maxmonsst_cmip6/cnrm_cm6_1_hr_ssp5_8_5_maxmonsst.nc` | 81 MB | Fig. 2, Fig. 5 |
| `maxmonsst_cmip6/cnrm_cm6_1_hr_ssp2_4_5_maxmonsst.nc` | 81 MB | (SSP2-4.5 variant) |
| `mintmat_2001-2009.nc` | ~1 GB | Fig. 5 (region grid mapping) |
| `Longhurst_Regions_2007.nc` | 14 MB | Longhurst province overlays |

The following smaller files are bundled directly in the repository
(`code/src/abstemp/data/`) and require no separate download:

| File | Used for |
|------|---------|
| `abstemp_reg_degvel.parquet` | Fig. 5, Fig. 7 — pre-computed per-region degree-velocity for 1985, 2019, 2095 |
| `all_cmip6_hists_ssp585.csv` | Fig. 3 — area-weighted SST histograms, 24 CMIP6 models, SSP5-8.5 |
| `all_cmip6_hists_ssp245.csv` | Fig. 3 — same for SSP2-4.5 |
| `all_ostia_hists.csv` | Fig. 3 — OSTIA area-weighted histograms for 1985–1990 and 2019–2023 |
| `growth_rates.csv` | Fig. 4 — observed phytoplankton growth rates vs. temperature |

---

## Capsule layout

```
code/                        Python source (installed as editable package)
│
├── run                      Entrypoint executed by Reproducible Run
│
└── src/abstemp/
    ├── __init__.py          Top-level API; vector_figs flag
    ├── figure_scripts/      One module per manuscript figure
    │   ├── __init__.py      generate_figures() — runs all figures, reports pass/fail
    │   ├── maxmonsst.py     Fig. 2 — warmest-month SST maps
    │   ├── global_hists.py  Fig. 3 — SST histograms
    │   ├── growth_models.py Fig. 4 — phytoplankton growth models
    │   ├── temp_velocities.py Fig. 5, 7 — temperature-velocity maps and histogram
    │   └── methods.py       Fig. 6 — OSTIA SST diagnostic maps
    ├── tempvel/             Degree-velocity computation (movedegree, regvel_*)
    ├── reg_calculations/    Region ↔ grid mapping (regvec_to_arr, arr_to_regvec)
    ├── data/                Bundled small datasets + open_warmest_* accessors
    └── seagrid/             SST download interfaces (OSTIA, CMIP6, Copernicus)

data/                        Large NetCDF input files (mounted at /data)
environment/
└── Dockerfile               Reproducible environment (pixi + conda-forge)
```

---

## Environment

Dependencies are managed with [pixi](https://pixi.sh) and pinned via
`code/pixi.lock`.  The full dependency list is in `code/pyproject.toml`.  Key
packages: `xarray`, `numpy`, `scipy`, `matplotlib`, `cartopy`, `projmap`,
`pyresample`, `pandas`, `scikit-learn`, `dask`, `netCDF4`.

The environment is built once into the container image (stored outside
`/code/` so it is not affected when Code Ocean mounts the code directory at
runtime).  No internet access is required during the run.

---

## Running locally with Docker or Podman

If you want to reproduce figures on your own machine rather than on Code Ocean:

```bash
# Clone the repository
git clone <repo-url>
cd code_ocean

# Build the image (~5–10 min on first run; subsequent builds use cache)
podman build -t abstemp -f environment/Dockerfile .

# Run (mount the data directory and an output directory)
mkdir -p figs
podman run --rm \
  -v "$(pwd)/data:/data:ro" \
  -v "$(pwd)/figs:/results" \
  abstemp
```

Replace `podman` with `docker` if using Docker Engine.  Figures are written to
`figs/` on the host.  No display is required (the `Agg` backend is used).

---

## Exploring the analysis interactively

If you want to step through the analysis rather than just running the full
figure pipeline, the Python package can be used interactively after installing
the pixi environment:

```bash
cd code
pixi install          # one-time setup
pixi run python       # launch Python inside the environment
```

```python
import abstemp
import abstemp.figure_scripts as figs

# Load the pre-computed degree-velocity table
df = abstemp.read_regdegvel()
print(df.columns.tolist())

# Reproduce individual figures
figs.global_histograms()     # Fig. 3
figs.growth_model_plot()     # Fig. 4
figs.tempvel_histogram()     # Fig. 7

# Explore the connectivity matrix
ds = abstemp.open_mintmat_ds()
print(ds)                    # dijkmintmat: travel times in days between ~11 000 regions
```

Large NetCDF files must be accessible via the `ABSTEMP_DATA_DIR` environment
variable (set to `/data` in the container) or placed in
`code/src/abstemp/data/` for local use.
