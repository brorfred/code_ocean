import numpy as np
import pandas as pd
import xarray as xr


def npz_to_ds(infile="data/mintmat_2001-2009.npz", outfile="data/mintmat_2001-2009.nc"):
    fh = np.load(infile)
    mintmat = fh["mintmat"].astype(np.float32)
    dijkmintmat = fh["dijkmintmat"].astype(np.float32)

    ds = xr.Dataset(
        data_vars={
            "mintmat": (("to_reg", "from_reg"), mintmat),
            "dijkmintmat": (("to_reg", "from_reg"), dijkmintmat),
            "regmat": (("lat", "lon"), fh["regmat"]),
        },
        coords={
            "reglat": (("to_reg",), fh["reglat"]),
            "reglon": (("to_reg",), fh["reglon"]),
        },
    )
    return ds


def to_dataframe(ds):
    df = pd.DataFrame(
        {
            "id": range(len(ds.reglat)),
            "lat": ds.reglat,
            "lon": ds.reglon,
        }
    )
    return df
