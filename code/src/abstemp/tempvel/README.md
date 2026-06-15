# abstemp.tempvel

Compute and visualise **temperature velocity** — the time (in years) for each ocean
region to reach a location that is currently 1 °C warmer.

## Concept

Temperature velocity is calculated using the Dijkstra minimum-travel-time matrix from
`mintmat_2001-2009.nc`.  For each region *r*, `movedegree(degree, sst)` finds the
nearest region (by travel time) whose SST is at least `degree` °C higher than *r*'s
current SST.  The travel time to reach that region is the degree-velocity.

## Public API

```python
from abstemp import tempvel

# Build the full per-region degree-velocity table (all three time periods)
df = tempvel.generate_regdegvel_df()

# Compute degree-velocity for a specific SST array
vel, bas, frm = tempvel.movedegree(degree=1, sst=df["2019_max_sst"])

# Add degree-velocity columns for a single time period
df = tempvel.regvel_1985(df)   # adds 1985_* columns
df = tempvel.regvel_2019(df)   # adds 2019_* columns
df = tempvel.regvel_2095(df)   # adds 2095_* columns

# Plot area-weighted histogram of degree-velocity
tempvel.regdegvel_hist(df)
```

## Function reference

| Function | Description |
|----------|-------------|
| `generate_regdegvel_df()` | Build complete degree-velocity DataFrame for 1985/2019/2095 |
| `movedegree(degree, sst)` | Core algorithm: return vel/bas/frm arrays for given SST |
| `regvel_1985(df)` | Add 1985–1990 columns to *df* |
| `regvel_2019(df)` | Add 2019–2023 columns to *df* |
| `regvel_2095(df)` | Add 2095–2100 columns (CNRM-CM6-1-HR SSP5-8.5 by default) to *df* |
| `regdegvel_hist(df)` | Area-weighted histogram of travel time (years) |

## Output columns (per time period `{year}`)

| Column | Description |
|--------|-------------|
| `{year}_max_sst` | Maximum monthly SST (°C) for that period |
| `{year}_regdegvel` | Travel time (days) to reach a 1°C-warmer region |
| `{year}_regdegbas` | Longhurst basin of the destination region |
| `{year}_regdegfrm` | Region index of the destination |

## Known issues

- `tempvel/climvel.py` is unmaintained legacy code and should not be used.
