"""CF Zarr oracles (wxylib stage S5b): a store laid out as xarray writes
an OISST-style month (time x zlev x lat x lon, packed int16 with
scale_factor and _FillValue, float32 1D coordinates, latitude ascending,
longitude 0 to 360), and what GDAL's warper writes for each day from a
north-up GeoTIFF of the same values.

The warps use GDAL's own Python bindings (osgeo.gdal). Needs osgeo.gdal,
numpy, pandas, xarray and zarr (3.x), for example from a conda-forge env
(micromamba create -n gdal -c conda-forge gdal numpy xarray zarr). Run from
the repository root:

    python scripts/make_cf_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import pandas as pd
import xarray as xr
from osgeo import gdal

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "cf")

NT, STEP = 3, 1.0
TARGET = dict(crs="EPSG:3031", bounds=(-4e6, -4e6, 4e6, 4e6), size=(100, 100))
METHODS = ["near", "bilinear"]


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
    xmin, ymin, xmax, ymax = TARGET["bounds"]
    dst_gt = [xmin, (xmax - xmin) / w, 0.0, ymax, 0.0, (ymin - ymax) / h]
    out = {"gdal": gdal.__version__,
           "xarray": xr.__version__, "store": "oisst-like.zarr", "array": "/sst",
           "crs": TARGET["crs"], "geotransform": dst_gt,
           "width": w, "height": h, "dtype": "int16", "nodata": -999,
           "days": NT, "results": []}
    tif = "/vsimem/day.tif"
    for m in METHODS:
        name = "oisst-like-to-3031-%s.bin" % m
        res = np.empty((NT, h, w), "int16")
        for t in range(NT):
            # the day as a north-up GeoTIFF (in memory), 0 to 360 east, 90 to -90 north
            src = gdal.GetDriverByName("GTiff").Create(tif, lon.size, lat.size, 1, gdal.GDT_Int16)
            src.SetGeoTransform([0.0, STEP, 0.0, 90.0, 0.0, -STEP])
            src.SetProjection("EPSG:4326")
            src.GetRasterBand(1).SetNoDataValue(-999)
            src.GetRasterBand(1).WriteArray(data[t, 0][::-1])
            # approximate transformer at 0.125 pixels (gdalwarp's default),
            # one thread, target starts as nodata
            dst = gdal.Warp("", src, format="MEM", dstSRS=TARGET["crs"],
                            outputBounds=TARGET["bounds"], width=w, height=h,
                            resampleAlg=m, dstNodata=-999, errorThreshold=0.125,
                            multithread=False, warpOptions=["INIT_DEST=NO_DATA"])
            res[t] = dst.ReadAsArray()
            src = dst = None
        res.astype("<i2").tofile(os.path.join(OUT, name))
        out["results"].append({"resampling": m, "file": name})
    gdal.Unlink(tif)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
