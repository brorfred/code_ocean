# abstemp.figure_scripts

Publication-ready figure generators for the abstemp manuscript.

## Usage

All figure functions save output to the `figs/` directory (created automatically
where needed).  Run from the project root so relative paths resolve correctly.

```python
import abstemp.figure_scripts as figs

figs.warmest_month_maps()   # Three-panel warmest-month SST maps
figs.global_histograms()    # SST histogram comparison (OSTIA + CMIP6)
figs.growth_model_plot()    # Blackford + Norberg-Eppley two-panel figure
figs.checkerboard()         # 2°×2° grid diagnostic map
figs.sst_maps()             # Max/min/range SST maps
```

Temperature-velocity figures are called from `temp_velocities` directly:

```python
from abstemp.figure_scripts.temp_velocities import tempvel_maps, tempvel_histogram
import abstemp
df = abstemp.read_regdegvel()
tempvel_maps(df)       # three-panel comparison map
tempvel_histogram()    # area-weighted velocity histogram
```

## Modules

### `figpref`

Matplotlib style presets for manuscript and presentation output.

```python
from abstemp.figure_scripts import figpref

figpref.manuscript()      # clean print-ready style
figpref.presentation()    # high-contrast dark style for slides
figpref.current()         # apply whichever style is set in figpref.set_current
```

Module-level constants `lon` and `lat` provide the 2°×2° grid coordinate arrays
used by map figures.

### `global_hists`

```python
from abstemp.figure_scripts.global_hists import main_histograms
main_histograms()              # SSP5-8.5 only
main_histograms(include_c24=True)  # also overlay SSP2-4.5
```

### `growth_models`

Growth-rate models (Blackford, Norberg-Eppley) with derivatives and noise envelopes.

```python
from abstemp.figure_scripts.growth_models import plot, growth_observations
plot()               # two-panel model comparison → figs/growth_models.pdf
growth_observations()  # observed growth rates + Blackford curve
```

### `maxmonsst`

```python
from abstemp.figure_scripts.maxmonsst import warmest_months, all_model_maps
warmest_months()    # three-panel map 1985/2019/2095 → figs/warmest_sst_map.*
all_model_maps()    # one PNG per model × experiment → figs/model_maps/
```

### `methods`

```python
from abstemp.figure_scripts.methods import checkerboard, sst_maps
checkerboard()   # 2°×2° grid diagnostic → figs/checkerboard.*
sst_maps()       # max/min/range SST panels → figs/ostia_sst.*
```

### `temp_velocities`

```python
from abstemp.figure_scripts.temp_velocities import (
    tempvel_map, tempvel_maps, tempvel_histogram
)
tempvel_map(df, year=2019)   # single-year velocity map
tempvel_maps(df)             # three-panel comparison map
tempvel_histogram()          # area-weighted velocity histogram
```

## Known issues

- `all_model_maps()` in `maxmonsst.py` creates an `mp` object that is
  never passed into the loop.  See ISSUES.md D5.
- `blackford_noise()` in `growth_models.py` computes values that are
  never returned.  See ISSUES.md D6.
