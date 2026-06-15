"""Figure-generation entry points for the abstemp manuscript.

Each public function in this module produces and saves one or more
publication-ready figures to the ``figs/`` directory.
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd

import cmap
import projmap

import abstemp
from abstemp import (
    data,
    warmest_month,
    global_sst_histograms,
    tempvel,
)

from .growth_models import plot as growth_model_plot
from .maxmonsst import warmest_months as warmest_month_maps
from .methods import checkerboard, sst_maps
from .global_hists import main_histograms as global_histograms
from .temp_velocities import tempvel_maps, tempvel_histogram

def generate_figures():
    """Run all manuscript figure generators, reporting pass/fail for each."""
    figures = [
        ("Figure 2 — warmest month maps",    warmest_month_maps),
        ("Figure 3 — global histograms",      global_histograms),
        ("Figure 4 — growth model",           growth_model_plot),
        ("Figure 5 — temperature velocities", tempvel_maps),
        ("Figure 6 — SST maps",               sst_maps),
        ("Figure 7 — tempvel histogram",      tempvel_histogram),
    ]
    results = {}
    for label, fn in figures:
        print(f"\n{label}")
        try:
            fn()
            results[label] = "OK"
        except Exception as exc:
            results[label] = f"FAILED — {exc}"
            plt.close("all")

    print("\n── Summary ──────────────────────────────")
    for label, status in results.items():
        print(f"  {'✓' if status == 'OK' else '✗'} {label}: {status}")
