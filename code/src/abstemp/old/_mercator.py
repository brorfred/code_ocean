
import pandas as pd


def retrieve(dtm="2020-01-01", data_var="uvel", force=False):
    """Download a single daily field from the Mercator THREDDS server and save it as Zarr.

    Opens the remote dataset via OPeNDAP, selects the time step closest to `dtm`,
    and writes the result to a local Zarr store.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object for the desired time step.
    data_var : str
        Variable key (one of 'uvel', 'vvel', 'wvel', 'temp', 'salt').
        force : bool
            If True, re-download the file even if it already exists locally.
    """
    dtm = pd.to_datetime(dtm)
    fname = filename(dtm=dtm, data_var=data_var)
    ddir  = datadir(dtm=dtm)
    if (ddir / fname).exists() and not force:
       return
    chunks = {'time_counter':1, 'deptht':-1, 'y':1000, 'x':1000}
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(
        max_retries=requests.adapters.Retry(total=5, backoff_factor=1),
        pool_connections=8,
        pool_maxsize=16
    )
    session.mount("http://", adapter)
    path = f"/thredds/dodsC/{PRODUCT}-daily-grid"
    url = HOST + path + vardict[data_var][0]
    store = PydapDataStore.open(url, session=session)
    ds = xr.open_dataset(store, chunks=chunks).astype('float32')

    dtm_index = pd.to_datetime(ds.time_counter).normalize()
    iloc = np.argmin(np.abs(dtm_index-dtm))
    da = ds[vardict[data_var][1]].isel(time_counter=iloc).to_dataset()

    write_job = da.to_zarr(str(ddir / fname), mode="w", compute=False)
    with ProgressBar():
        write_job.compute()
