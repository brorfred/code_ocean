"""Interactive world map of per-region degree-velocity data using Plotly."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

import abstemp

_MONTH_NAMES = [
    "", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]

_LONGHURST_NAMES = {
    1.0: "Coastal",
    2.0: "Trades",
    3.0: "Westerlies",
    4.0: "Polar",
    5.0: "Other",
}

# JS injected after Plotly renders — highlights the base dot on hover
_POST_SCRIPT = """
(function() {
    var baseLats = %s;
    var baseLons = %s;
    var baseNregs = %s;

    gd.on('plotly_hover', function(data) {
        var pt = data.points[0];
        if (pt.curveNumber !== 0) return;
        var cd = pt.customdata;
        var baseIdx = cd[0];
        if (baseIdx == null || isNaN(baseIdx)) return;
        var lat = baseLats[baseIdx];
        var lon = baseLons[baseIdx];
        if (lat == null || lon == null) return;
        Plotly.restyle(gd, {lat: [[lat]], lon: [[lon]], visible: true}, [1]);
    });

    gd.on('plotly_unhover', function() {
        Plotly.restyle(gd, {visible: false}, [1]);
    });
})();
"""


def _build_dataframe() -> pd.DataFrame:
    ds = abstemp.open_mintmat_ds()
    df = abstemp.read_regdegvel()

    df["clm_max_sst"] = ds.clm_max_sst.values
    df["clm_mon_sst"] = ds.clm_mon_sst.values

    return df.dropna(subset=["reglat", "reglon"])


def _fmt(val, decimals=2, suffix=""):
    if val is None or (isinstance(val, float) and not np.isfinite(val)):
        return "—"
    return f"{val:.{decimals}f}{suffix}"


def _fmt_month(val):
    if val is None or (isinstance(val, float) and not np.isfinite(val)):
        return "—"
    m = int(round(val))
    if 1 <= m <= 12:
        return _MONTH_NAMES[m]
    return f"{val:.1f}"


def _build_hover(df: pd.DataFrame, full_df: pd.DataFrame) -> list[str]:
    texts = []
    for idx, row in df.iterrows():
        base_idx = row["regdegbas"]
        vel_idx = row["regdegvel"]

        base_nreg = "—"
        base_lat = "—"
        base_lon = "—"
        if np.isfinite(base_idx):
            b = full_df.loc[int(base_idx)]
            base_nreg = str(int(base_idx))
            base_lat = f"{b['reglat']:.2f}"
            base_lon = f"{b['reglon']:.2f}"

        vel_nreg = "—"
        if np.isfinite(vel_idx):
            vel_nreg = str(int(vel_idx))

        longh = row["longhurst_basins"]
        longh_str = _LONGHURST_NAMES.get(longh, "—") if np.isfinite(longh) else "—"

        parts = [
            f"<b>nreg {int(idx)}  lat {row['reglat']:.2f}  lon {row['reglon']:.2f}</b>",
            f"Absolute max SST (°C): {_fmt(row['abs_max_sst'])}",
            f"Month of absolute max SST: {_fmt_month(row['abs_mon_sst'])}",
            f"Climatological max SST (°C): {_fmt(row['clm_max_sst'])}",
            f"Month of climatological max SST: {_fmt_month(row['clm_mon_sst'])}",
            f"Longhurst basin: {longh_str}",
            f"Degree-velocity → dot {vel_nreg}",
            f"<b>Base dot: nreg {base_nreg}  lat {base_lat}  lon {base_lon}</b>",
        ]
        texts.append("<br>".join(parts))
    return texts


def world_map(output_html: str | None = "degvel_map.html") -> go.Figure:
    """Create an interactive world map of per-region degree-velocity data.

    Each dot is coloured by its ``regdegvel`` value.  Hovering shows full SST
    statistics and, if ``regdegbas`` is defined for that dot, highlights the
    corresponding base dot on the map with a red ring.

    Parameters
    ----------
    output_html : str or None
        Path for a self-contained HTML export.  Defaults to
        ``"degvel_map.html"``.  Pass ``None`` to skip saving.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    df = _build_dataframe()

    # Build full-index lookup arrays for the JS callback
    full_df = abstemp.read_regdegvel()
    max_idx = full_df.index.max()
    lat_arr = [None] * (max_idx + 1)
    lon_arr = [None] * (max_idx + 1)
    nreg_arr = [None] * (max_idx + 1)
    for idx, row in full_df.iterrows():
        if np.isfinite(row["reglat"]) and np.isfinite(row["reglon"]):
            lat_arr[idx] = row["reglat"]
            lon_arr[idx] = row["reglon"]
            nreg_arr[idx] = int(idx)

    hover = _build_hover(df, full_df)

    # customdata[0] = regdegbas index (for JS lookup of base dot location)
    custom = df["regdegbas"].where(df["regdegbas"].notna(), other=float("nan")).values
    custom = [[int(v) if np.isfinite(v) else None] for v in custom]

    main_trace = go.Scattergeo(
        lat=df["reglat"],
        lon=df["reglon"],
        mode="markers",
        marker=dict(
            color=df["regdegvel"],
            colorscale="Viridis",
            size=3,
            opacity=0.85,
            colorbar=dict(title="Degree-velocity (dot index)", thickness=15),
            showscale=True,
        ),
        text=hover,
        hovertemplate="%{text}<extra></extra>",
        customdata=custom,
        name="Degree-velocity",
    )

    # Invisible trace used as a highlight ring for the base dot
    base_indicator = go.Scattergeo(
        lat=[None],
        lon=[None],
        mode="markers",
        marker=dict(
            symbol="circle-open",
            size=14,
            color="red",
            line=dict(width=3, color="red"),
        ),
        hoverinfo="skip",
        showlegend=False,
        visible=False,
        name="Base dot",
    )

    fig = go.Figure([main_trace, base_indicator])
    fig.update_layout(
        title="Per-region degree-velocity (regdegvel) — hover to highlight base dot",
        geo=dict(
            projection_type="natural earth",
            showland=True,
            landcolor="rgb(230, 230, 220)",
            showocean=True,
            oceancolor="rgb(200, 220, 240)",
            showcoastlines=True,
            coastlinecolor="rgb(100, 100, 100)",
            showframe=False,
        ),
        margin=dict(l=0, r=0, t=60, b=0),
    )

    if output_html is not None:
        import json
        post_script = _POST_SCRIPT % (
            json.dumps(lat_arr),
            json.dumps(lon_arr),
            json.dumps(nreg_arr),
        )
        fig.write_html(
            output_html,
            include_plotlyjs="cdn",
            post_script=post_script,
        )
        print(f"Saved to {output_html}")

    return fig


if __name__ == "__main__":
    fig = world_map()
    fig.show()
