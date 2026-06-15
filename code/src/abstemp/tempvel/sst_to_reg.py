
from abstemp.seagrid import ostia

import numpy as np
import xarray as xr
import pandas as pd
from sklearn.neighbors import BallTree

from abstemp import tempvel


def sst_ds_to_df(sst_ds):
    lons, lats = np.meshgrid(sst_ds.lon, sst_ds.lat)
    sst_df = pd.DataFrame({
        "lat":lats.flat,
        "lon":lons.flat,
        "sst":sst_ds.sst.data.flat
    }).dropna()
    return sst_df


def find_nearest_patch(reg_df, sst):
    tree = BallTree(np.deg2rad(reg_df.iloc[1:,[1,2]]), leaf_size=2, metric="haversine")
    sst_df = sst if isinstance(sst, pd.DataFrame) else sst_ds_to_df(sst)
    dist, ind = tree.query(np.deg2rad(sst_df[["lat","lon"]]), k=1)
    return ind + 1

def sst_to_patch(reg_df, sst_ds, reg=None):
    sst_df = sst_ds_to_df(sst_ds)
    reg = find_nearest_patch(reg_df, sst_df) if reg is None else reg
    sst_df["reg"] = reg
    regsst = sst_df.groupby("reg").mean()
    return regsst

"""
from njord import ostia
from abstemp.tempvel import sst_to_reg
from abstemp import tempvel

ds = xr.open_dataset("data/mintmat_2001-2009.nc")
df = tempvel.to_dataframe(ds)
ostiareg = pd.read_parquet("ostia_reg_1x_1x.parquet")

sstlist = []
dtmvec = pd.date_range("2019-01-01", "2023-12-31")
for dtm in dtmvec:
    sst = ostia.open_dataset(dtm)
    regsst = sst_to_reg.sst_to_patch(df, sst, ostiareg["reg"])
    sstlist.append(regsst["sst"].values)
    print(dtm)
sstarr = np.array(sstlist)
ds = xr.Dataset(
    data_vars={"sst":(("date", "reg"), sstarr)},
    coords={"date":dtmvec, "reg":df.id[1:]}
    )
"""
