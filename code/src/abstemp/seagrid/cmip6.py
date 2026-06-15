import glob
import pathlib
import zipfile

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import projmap
import cdsapi

DATADIR = pathlib.Path(__file__).absolute().parent
DATADIR = DATADIR.parent / "data" / "cmip6"
DATADIR.mkdir(parents=True, exist_ok=True)


#mlist = ["cams_csm1_0", 'ec_earth3_cc','cesm2', 'hadgem3_gc31_ll',
#          "canesm5_canoe", 'ipsl_cm6a_lr', 'kiost_esm', 'miroc6',
#          'miroc_es2l', 'noresm2_mm',]

#old_models = {#"cams_csm1_0",
#            'ec_earth3_cc':"EC-Earth3-CC",
#            'cesm2':"CESM2",
#            'hadgem3_gc31_ll':"HadGEM3-GC31-LL",
#            'canesm5_canoe':"CanESM5-CanOE",
#            'ipsl_cm6a_lr':"IPSL-CM6A-LR",
#            'kiost_esm':"KIOST-ESM",
#            'miroc6':"MIROC6",
#            'miroc_es2l':"MIROC-ES2L",
#            'noresm2_mm':"NorESM2-MM",}
#fail 'cams_csm1_0',

cmip6_sst_models = {
    "access_cm2":       "ACCESS-CM2",
    "cams_csm1_0":      "CAMS-CSM1-0",
    "canesm5_canoe":    "CanESM5-CanOE",
    "cesm2":            "CESM2",
    #"cesm2_waccm":      "CESM2-WACCM",
    "ciesm":            "CIESM",
    "cmcc_cm2_sr5":      "CMCC-CM2-SR5",
    "cmcc_esm2":        "CMCC-ESM2",
    "cnrm_cm6_1":       "CNRM-CM6-1",
    "cnrm_cm6_1_hr":    "CNRM-CM6-1-HR",
    "cnrm_esm2_1":      "CNRM-ESM2-1",
    'ec_earth3_cc':     "EC-Earth3-CC",
    #'ec_earth3':        "EC-Earth3",
    "fgoals_f3_l":      "FGOALS-f3-L",
    "fgoals_g3":        "FGOALS-g3",
    "fio_esm_2_0":      "FIO-ESM-2-0",
    #"giss_e2_1_g":      "GISS-E2-1-G",
    "hadgem3_gc31_ll":  "HadGEM3-GC31-LL",
    "ipsl_cm6a_lr":     "IPSL-CM6A-LR",
    "kiost_esm":        "KIOST-ESM",
    "mcm_ua_1_0":       "MCM-UA-1-0",
    "miroc_es2l":       "MIROC-ES2L",
    "miroc6":           "MIROC6",
    #"mpi_esm1_2_hr":    "MPI-ESM1-2-HR",
    #"mpi_esm1-2_lr":    "MPI-ESM1-2-LR",
    #"mri_esm2_0":       "MRI-ESM2-0",
    "nesm3":            "NESM3",
    "noresm2_lm":       "NorESM2-LM",
    "noresm2_mm":       "NorESM2-MM",
    #"taiesm1":          "TaiESM1",
    "ukesm1_0_ll":      "UKESM1-0-LL",
}

experiments = ["ssp2_4_5", "ssp5_8_5"]

def insensitive_glob(pattern: str, datadir: "pathlib.Path | None" = None) -> list[str]:
    """Return files in *datadir* whose names contain *pattern* (case-insensitive).

    Parameters
    ----------
    pattern : str
        Substring to search for in file names.
    datadir : pathlib.Path or None, optional
        Directory to search.  Defaults to the package CMIP6 data directory.

    Returns
    -------
    list of str
        Matching file paths with hyphens replaced by underscores.
    """
    datadir = datadir or DATADIR
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    #return glob.glob(''.join(map(either, f"*{pattern}*")))
    #print(''.join(map(either, f"*{pattern}*")))
    #return pathlib.Path(datadir).glob(''.join(map(either, f"*{pattern}*")))
    flist = glob.glob(f"{datadir}/*")
    flist = [s.replace("-", "_") for s in flist]
    return [s for s in flist if pattern.lower() in s.lower()]

def download(model: str = "ec_earth3_cc", experiment: str = "ssp5_8_5") -> None:
    """Download a CMIP6 SST file from the CDS API and store it locally.

    Retrieves monthly sea-surface temperature for years 2090–2100 from the
    Copernicus Climate Data Store, extracts the NetCDF file from the
    downloaded ZIP archive, and moves it into the package CMIP6 data
    directory.

    Parameters
    ----------
    model : str, optional
        Model key as it appears in :data:`cmip6_sst_models`.
        Default is ``"ec_earth3_cc"``.
    experiment : str, optional
        SSP experiment identifier (e.g. ``"ssp5_8_5"`` or ``"ssp2_4_5"``).
        Default is ``"ssp5_8_5"``.

    Returns
    -------
    None
    """
    dataset = "projections-cmip6"
    request = {
        'temporal_resolution': 'monthly',
        'experiment': experiment,
        'variable': 'sea_surface_temperature',
        'model': model,
        'year': ['2090', '2091', '2092', '2093', '2094',
                 '2095', '2096', '2097', '2098', '2099', '2100'],
        'month': ['01', '02', '03', '04', '05', '06',
                  '07', '08', '09', '10', '11', '12'],
        "format": "netcdf"}
    zipfname = pathlib.Path(f"{model}.zip")
    client = cdsapi.Client()
    client.retrieve(dataset, request, str(zipfname))
    zf = zipfile.ZipFile(zipfname)
    fobj = [fobj
            for fobj in zf.filelist
                if fobj.filename.endswith(".nc")
                ][0]
    ncfname = pathlib.Path(zf.extract(fobj))
    ncfname.rename(DATADIR / fobj.filename)
    zipfname.unlink()


def center_on_gmt(ds):
    """Roll dataset so longitudes span -180 to 180 with 0° at center.

    Handles both 1-D regular grids (lon coordinate) and 2-D curvilinear
    grids (lons data variable).  For 1-D grids the array is also sorted so
    longitude increases monotonically from -180 to 180, and the 2-D lons/lats
    meshgrid coordinates are rebuilt to stay consistent.
    """
    lon_key = "lon" if "lon" in ds.coords and ds["lon"].ndim == 1 else "lons" if "lons" in ds else None
    if lon_key is None:
        return ds

    if ds[lon_key].ndim == 1:
        # Regular 1-D grid: return early if already centred.  Use max > 180
        # (not min < 0) because some grids start slightly negative (e.g. -0.9°)
        # but still span mostly 0–360 and need centering.
        if ds[lon_key].max() <= 180:
            return ds
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180)
        ds = ds.sortby("lon")
        # Rebuild the 2-D lons/lats meshgrid so it matches the sorted lon.
        if "lons" in ds and "lats" in ds:
            lons_2d, lats_2d = np.meshgrid(ds.lon.values, ds.lat.values)
            ds = ds.assign_coords(
                lons=(ds.lons.dims, lons_2d),
                lats=(ds.lats.dims, lats_2d),
            )
    else:
        # 2-D curvilinear grid: centre lons if needed, then always mask seam.
        # Tripolar grids always have a seam (at the dateline or the tripolar
        # cap) where adjacent cells jump ~360° in longitude; pcolor renders
        # that column as a globe-spanning quad that paints SST over land.
        if ds[lon_key].max() > 180:
            lons_norm = ((ds.lons + 180) % 360) - 180
            i_dim = ds.lons.dims[-1]
            j_mid = ds.lons.shape[0] // 2
            roll_by = int(np.argmin(lons_norm[j_mid, :].values))
            ds = ds.roll({i_dim: -roll_by}, roll_coords=True)
            ds = ds.assign_coords(lons=((ds.lons + 180) % 360) - 180)

        lon_vals = ds.lons.values
        seam = np.abs(np.diff(lon_vals, axis=-1)) > 180   # (ny, nx-1)
        seam_mask = np.zeros(lon_vals.shape, bool)
        seam_mask[:, :-1] |= seam
        seam_mask[:, 1:] |= seam
        # pcolor computes cell corners as 4-point averages of adjacent centres;
        # expand by 1 to prevent corrupted corners when the seam shifts columns
        # between adjacent rows.
        expanded = seam_mask.copy()
        expanded[1:, :] |= seam_mask[:-1, :]
        expanded[:-1, :] |= seam_mask[1:, :]
        expanded[:, 1:] |= seam_mask[:, :-1]
        expanded[:, :-1] |= seam_mask[:, 1:]
        seam_mask = expanded
        if seam_mask.any():
            spatial_dims = set(ds.lons.dims)
            for var in list(ds.data_vars):
                var_dims = set(ds[var].dims)
                # Skip auxiliary arrays (e.g. bounds with a nvertex dimension).
                if spatial_dims.issubset(var_dims) and not (var_dims - spatial_dims - {"time"}):
                    ds[var] = ds[var].where(~seam_mask)

    return ds


def open_dataset(
    model: str = "cesm2",
    experiment: str = "ssp5_8_5",
    datadir: "pathlib.Path | None" = None,
) -> xr.Dataset:
    """Open a CMIP6 SST file and return it as a normalised xarray Dataset.

    Renames model-specific coordinate/variable names to the common names
    ``lat``, ``lon``, and ``sst``.  2-D curvilinear grids are renamed to
    ``lats``/``lons`` and a regular ``lat``/``lon`` meshgrid is added.
    The dataset is centred on the Greenwich meridian via
    :func:`center_on_gmt`.

    Parameters
    ----------
    model : str, optional
        Model key as it appears in :data:`cmip6_sst_models`.
        Default is ``"cesm2"``.
    experiment : str, optional
        SSP experiment identifier.  Default is ``"ssp5_8_5"``.
    datadir : pathlib.Path or None, optional
        Directory containing CMIP6 NetCDF files.  Defaults to the
        package CMIP6 data directory.

    Returns
    -------
    xarray.Dataset
        Dataset with variable ``sst`` and coordinates ``lat``, ``lon``
        (and optionally ``lats``, ``lons`` for curvilinear grids).
    """
    datadir = datadir or DATADIR
    mod = cmip6_sst_models[model]
    experiment = "".join(experiment.split("_")).lower()
    fname = f"tos_Omon_{mod}_{experiment}*"

    filename = glob.glob(f"{datadir}/{fname}")[0]
    ds = xr.open_dataset(filename) #pathlib.Path(datadir) / fname)
    if "latitude" in ds:
        ds = ds.rename({"latitude":"lat", "longitude":"lon"})
    elif "nav_lat" in ds:
        ds = ds.rename({"nav_lat":"lat", "nav_lon":"lon"})
    if "tos" in ds:
        ds = ds.rename({"tos":"sst"})
    if ds.lat.ndim == 2:
        ds = ds.rename({"lat":"lats", "lon":"lons"})
    if not "lats" in ds:
        lons,lats = np.meshgrid(ds.lon, ds.lat)
        ds = ds.assign_coords(lons=(("y","x"),lons), lats=(("y","x"),lats))
    return center_on_gmt(ds)

def retrieve_all_files() -> None:
    """Download all CMIP6 SST files for every model and both SSP experiments.

    Iterates over :data:`cmip6_sst_models` and :data:`experiments`,
    calling :func:`download` for each combination.

    Returns
    -------
    None
    """
    for model in cmip6_sst_models:
        for experiment in experiments:
            download(model=model, experiment=experiment)
            print(f"{model} -- {experiment}")
