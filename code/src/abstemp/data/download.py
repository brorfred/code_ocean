import pathlib
import os

import requests
from tqdm import tqdm

from abstemp.seagrid import cmip6

REPO_IP = "157.180.67.142"

def retrieve_all_cmip6_files() -> None:
    """Download all CMIP6 SST files for every registered model and experiment.

    Iterates over :data:`abstemp.seagrid.cmip6.cmip6_sst_models` and
    :data:`abstemp.seagrid.cmip6.experiments`, calling
    :func:`abstemp.seagrid.cmip6.download` for each combination.

    Returns
    -------
    None
    """
    for model in cmip6.cmip6_sst_models:
        for experiment in cmip6.experiments:
            cmip6.download(model=model, experiment=experiment)
            print(f"{model} -- {experiment}")

def http_retrieve(
    url: str,
    params: dict = {"download": 1},
    datadir: pathlib.Path = pathlib.Path(__file__).parent,
    filename: str = "download.nc",
    token: str | None = None,
) -> None:
    """Download a file over HTTP with a progress bar.

    Parameters
    ----------
    url : str
        Full URL of the file to download.
    params : dict, optional
        Query parameters appended to the request.  Default is
        ``{"download": 1}`` (Zenodo download trigger).
    datadir : pathlib.Path, optional
        Destination directory.  Defaults to the package data directory.
    filename : str, optional
        Name for the saved file.  Default is ``"download.nc"``.
    token : str or None, optional
        Zenodo personal access token for restricted records.

    Returns
    -------
    None
    """
    headers = {"Authorization": f"Bearer {token}"} if token else {}
    response = requests.get(url, params=params, headers=headers, stream=True, timeout=60)
    response.raise_for_status()
    if "text/html" in response.headers.get("Content-Type", ""):
        raise ValueError(
            "Received an HTML page instead of a file — the URL may be "
            "restricted. Pass a Zenodo personal access token via the `token` argument."
        )
    dest = pathlib.Path(datadir) / filename

    total = int(response.headers.get("content-length", 0))
    with open(dest, "wb") as f, tqdm(
        desc=os.path.basename(dest),
        total=total,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            bar.update(len(chunk))

    print(f"Saved to {dest}")


def mintmat(
    url: str = "https://zenodo.org/records/20497223/files/mintmat_2001-2009.nc",
    datadir: pathlib.Path = pathlib.Path(__file__).parent,
    token: str | None = None,
) -> None:
    """Download the mintmat connectivity dataset from Zenodo.

    Parameters
    ----------
    url : str, optional
        Zenodo download URL for the mintmat file.
    datadir : pathlib.Path, optional
        Destination directory.  Defaults to the package data directory.
    token : str or None, optional
        Zenodo personal access token for restricted records.

    Returns
    -------
    None
    """
    http_retrieve(
        url=url,
        params={"download": 1},
        datadir=datadir,
        filename="mintmat_2001-2009.nc",
        token=token,
    )


def maxmonsst_fields() -> None:
    """Download all pre-computed warmest-month SST NetCDF files from the project server.

    Downloads:

    - ``ostia_maxmonsst_1985-1990.nc`` — OSTIA reanalysis max-month stats 1985–1990
    - ``ostia_maxmonsst_2019-2023.nc`` — OSTIA NRT max-month stats 2019–2023
    - ``maxmonsst_cmip6/cnrm_cm6_1_hr_ssp5_8_5_maxmonsst.nc`` — CNRM-CM6-1-HR SSP5-8.5 max-month stats 2095–2100
    - ``maxmonsst_cmip6/cnrm_cm6_1_hr_ssp2_4_5_maxmonsst.nc`` — CNRM-CM6-1-HR SSP2-4.5 max-month stats 2095–2100

    Returns
    -------
    None
    """
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/ostia_maxmonsst_1985-1990.nc",
        filename="ostia_maxmonsst_1985-1990.nc",
    )
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/ostia_maxmonsst_2019-2023.nc",
        filename="ostia_maxmonsst_2019-2023.nc",
    )
    datadir = pathlib.Path(__file__).parent / "maxmonsst_cmip6"
    datadir.mkdir(exist_ok=True, parents=True)
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/cnrm_cm6_1_hr_ssp5_8_5_maxmonsst.nc",
        filename=str(datadir / "cnrm_cm6_1_hr_ssp5_8_5_maxmonsst.nc"),
    )
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/cnrm_cm6_1_hr_ssp2_4_5_maxmonsst.nc",
        filename=str(datadir / "cnrm_cm6_1_hr_ssp2_4_5_maxmonsst.nc"),
    )

def ostia_sst_fields() -> None:
    """Download the raw monthly OSTIA SST NetCDF files from the project server.

    Only needed to regenerate the max-month statistics via
    :func:`abstemp.data.generate_maxmonsst_files`.  Downloads:

    - ``ostia_sst_1985-1990.nc`` — OSTIA reanalysis monthly SST 1985–1990
    - ``ostia_sst_2019-2023.nc`` — OSTIA NRT monthly SST 2019–2023

    Returns
    -------
    None
    """
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/ostia_sst_1985-1990.nc",
        filename="ostia_sst_1985-1990.nc",
    )
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/ostia_sst_2019-2023.nc",
        filename="ostia_sst_2019-2023.nc",
    )

def longhurst_regions() -> None:
    """Download the Longhurst province mask (on OSTIA grid) from the project server.

    Saves ``Longhurst_Regions_2007.nc`` to the package data directory.

    Returns
    -------
    None
    """
    http_retrieve(
        url=f"http://{REPO_IP}/abstemp/Longhurst_Regions_2007.nc",
        filename="Longhurst_Regions_2007.nc",
    )
