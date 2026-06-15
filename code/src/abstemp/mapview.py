"""Interactive world map of per-region SST statistics using Plotly."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

import abstemp


_SST_VARS = {
    "abs_max_sst": "Absolute max SST (°C)",
    "abs_mon_sst": "Month of absolute max SST",
    "clm_max_sst": "Climatological max SST (°C)",
    "clm_mon_sst": "Month of climatological max SST",
}


def _build_dataframe(ds) -> pd.DataFrame:
    """Extract region coordinates and SST variables into a flat DataFrame.

    Parameters
    ----------
    ds : xarray.Dataset
        The mintmat dataset with ``reglat``, ``reglon``, and the SST
        variables listed in ``_SST_VARS``.

    Returns
    -------
    pandas.DataFrame
        One row per region with columns ``lat``, ``lon``, ``nreg``, and one
        column per SST variable.  Region 0 (NaN coordinates) is dropped.
    """
    df = pd.DataFrame(
        {
            "lat": ds.reglat.values,
            "lon": ds.reglon.values,
            "nreg": ds.nreg.values,
        }
    )
    for var in _SST_VARS:
        df[var] = ds[var].values

    # Region 0 has NaN coords — drop it
    df = df.dropna(subset=["lat", "lon"])
    return df


def _make_hover(df: pd.DataFrame) -> list[str]:
    """Build hover-text strings for each region.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame returned by :func:`_build_dataframe`.

    Returns
    -------
    list of str
        One HTML hover string per row.
    """
    texts = []
    for _, row in df.iterrows():
        parts = [f"<b>ij {int(row.nreg)}  lat {row.lat:.2f}  lon {row.lon:.2f}</b>"]
        for var, label in _SST_VARS.items():
            val = row[var]
            if np.isfinite(val):
                parts.append(f"{label}: {val:.2f}")
            else:
                parts.append(f"{label}: —")
        texts.append("<br>".join(parts))
    return texts


def world_map(variable: str = "abs_max_sst", output_html: str | None = None) -> go.Figure:
    """Create an interactive world map of a per-region SST variable.

    Renders ~11 000 region centroids as a scatter plot on a natural-earth
    projection.  Hovering over a point shows all four SST variables for that
    region.  A dropdown menu lets the user switch between variables without
    reloading the page.

    Parameters
    ----------
    variable : str, optional
        Name of the SST variable to display on load.  Must be one of
        ``abs_max_sst``, ``abs_mon_sst``, ``clm_max_sst``, ``clm_mon_sst``.
        Default is ``"abs_max_sst"``.
    output_html : str or None, optional
        If given, the figure is saved as a self-contained HTML file at this
        path in addition to being returned.

    Returns
    -------
    plotly.graph_objects.Figure
        The interactive figure.  Call ``fig.show()`` to open it in a browser,
        or ``fig.write_html(path)`` to save it.
    """
    if variable not in _SST_VARS:
        raise ValueError(f"variable must be one of {list(_SST_VARS)}, got {variable!r}")

    ds = abstemp.open_mintmat_ds()
    df = _build_dataframe(ds)
    hover = _make_hover(df)

    # One trace per variable; only the selected one is visible at load time
    traces = []
    for var, label in _SST_VARS.items():
        traces.append(
            go.Scattergeo(
                lat=df["lat"],
                lon=df["lon"],
                mode="markers",
                marker=dict(
                    color=df[var],
                    colorscale="RdBu_r",
                    size=3,
                    opacity=0.85,
                    colorbar=dict(title=label, thickness=15),
                    showscale=True,
                ),
                text=hover,
                hovertemplate="%{text}<extra></extra>",
                name=label,
                visible=(var == variable),
            )
        )

    # Dropdown buttons: each button makes exactly one trace visible
    buttons = []
    for i, (var, label) in enumerate(_SST_VARS.items()):
        visibility = [v == var for v in _SST_VARS]
        buttons.append(
            dict(
                label=label,
                method="update",
                args=[{"visible": visibility}],
            )
        )

    fig = go.Figure(traces)
    fig.update_layout(
        title="Per-region SST statistics (2001–2008)",
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
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=0.01,
                xanchor="left",
                y=1.12,
                yanchor="top",
            )
        ],
        margin=dict(l=0, r=0, t=60, b=0),
    )

    if output_html is not None:
        fig.write_html(output_html, include_plotlyjs="cdn")
        print(f"Saved to {output_html}")

    return fig


if __name__ == "__main__":
    fig = world_map(output_html="sst_map.html")
    fig.show()
