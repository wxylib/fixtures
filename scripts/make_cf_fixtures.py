"""CF Zarr oracles (wxylib stage S5b): a store laid out as xarray writes
an OISST-style month (time x zlev x lat x lon, packed int16 with
scale_factor and _FillValue, float32 1D coordinates, latitude ascending,
longitude 0 to 360), and what GDAL's warper writes for each day from a
north-up GeoTIFF of the same values.

Needs xarray and zarr (3.x). Run from the repository root:

    python scripts/make_cf_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.transform import from_bounds, from_origin
from rasterio.warp import Resampling, reproject

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "cf")

NT, STEP = 3, 1.0
TARGET = dict(crs="EPSG:3031", bounds=(-4e6, -4e6, 4e6, 4e6), size=(100, 100))
METHODS = {"near": Resampling.nearest, "bilinear": Resampling.bilinear}


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    lat = np.arange(-90 + STEP / 2, 90, STEP, dtype="float32")
    lon = np.arange(STEP / 2, 360, STEP, dtype="float32")
    LON, LAT = np.meshgrid(lon.astype("float64"), lat.astype("float64"))
    data = np.empty((NT, 1, lat.size, lon.size), "int16")
    for t in range(NT):
        v = 2500 * np.cos(np.radians(LAT)) + 300 * np.sin(np.radians(LON) * (t + 1))
        v = np.round(v).astype("int16")
        v[LAT < -78] = -999
        data[t, 0] = v
    ds = xr.Dataset(
        {"sst": (("time", "zlev", "lat", "lon"), data,
                 {"scale_factor": 0.01, "add_offset": 0.0, "_FillValue": np.int16(-999),
                  "units": "Celsius"})},
        coords={
            "time": pd.date_range("1981-09-01", periods=NT),
            "zlev": ("zlev", np.array([0], "float32")),
            "lat": ("lat", lat, {"units": "degrees_north", "standard_name": "latitude"}),
            "lon": ("lon", lon, {"units": "degrees_east", "standard_name": "longitude"}),
        },
    )
    ds.to_zarr(os.path.join(OUT, "oisst-like.zarr"), zarr_format=3, mode="w",
               encoding={"sst": {"chunks": (2, 1, lat.size, lon.size)}})

    w, h = TARGET["size"]
    dst_transform = from_bounds(*TARGET["bounds"], w, h)
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__,
           "xarray": xr.__version__, "store": "oisst-like.zarr", "array": "/sst",
           "crs": TARGET["crs"], "geotransform": list(dst_transform.to_gdal()),
           "width": w, "height": h, "dtype": "int16", "nodata": -999,
           "days": NT, "results": []}
    tif = os.path.join(OUT, "day.tif")
    for m, alg in METHODS.items():
        name = "oisst-like-to-3031-%s.bin" % m
        res = np.empty((NT, h, w), "int16")
        for t in range(NT):
            with rasterio.open(tif, "w", driver="GTiff", width=lon.size, height=lat.size,
                               count=1, dtype="int16", crs="EPSG:4326", nodata=-999,
                               transform=from_origin(0, 90, STEP, STEP)) as d:
                d.write(data[t, 0][::-1], 1)
            dst = np.zeros((h, w), "int16")
            with rasterio.open(tif) as src:
                reproject(rasterio.band(src, 1), dst, dst_transform=dst_transform,
                          dst_crs=TARGET["crs"], resampling=alg, dst_nodata=-999,
                          num_threads=1)
            res[t] = dst
        res.astype("<i2").tofile(os.path.join(OUT, name))
        out["results"].append({"resampling": m, "file": name})
    os.remove(tif)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
