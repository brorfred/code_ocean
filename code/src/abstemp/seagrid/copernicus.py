"""Download interface for Copernicus Marine Data Store.

Uses the ``copernicusmarine`` library to retrieve original files by date.
"""

import pandas as pd


class DownloadError(Exception):
    """Raised when a Copernicus Marine download fails."""

    pass

def get(dtm, settings, **kw):
    """Download a file from the Copernicus Marine Data Store.

    Parameters
    ----------
    dtm : str or datetime-like
        Date of the file to download.
    settings : dict-like
        Configuration with keys ``'dataset_id'``, ``'datadir'``,
        ``'username'``, ``'password'``, ``'file_filter'``.
    **kw
        Overrides for any *settings* key.

    Raises
    ------
    DownloadError
        If no file was found on the server.
    """
    import copernicusmarine
    sdict = {}
    for key in ["dataset_id", "datadir", "username", "password", "file_filter"]:
        sdict[key] = kw.get(key, settings[key])
    dtm = pd.to_datetime(dtm)
    response = copernicusmarine.get(
        dataset_id=sdict["dataset_id"],
        filter=f"*{dtm.strftime(sdict['file_filter'])}*",
        output_directory = sdict["datadir"],
        username=sdict["username"],
        password=sdict["password"],
        #force_download=True,
        no_directories=True,
        #overwrite_output_data=True,
        disable_progress_bar=True,
        #service="original-files",
        #dataset_version="default"
        )
    if response.number_of_files_to_download == 0:
        raise DownloadError("No file was found on server.")
