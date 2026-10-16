"""Warp oracles: what GDAL's warper (through rasterio.warp.reproject, with
the approximate transformer at its default 0.125 pixel error threshold)
writes for a few source fixtures, target grids and resampling methods.

Target grids are near 1:1 with the source resolution, which is the bar
rwarp is validated at (downsampled interpolated warps differ by design
until overview selection and kernel scaling land, stage S5).

Run from the repository root after make_cog_fixtures.py:

    python scripts/make_warp_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COG = os.path.join(ROOT, "cog")
OUT = os.path.join(ROOT, "warp")

CASES = [
    # global lon/lat to south polar stereographic, about 55 km cells
    dict(name="world-to-3031", source="world-int16-deflate.tif",
         crs="EPSG:3031", bounds=(-4.0e6, -4.0e6, 4.0e6, 4.0e6), size=(145, 145)),
    # polar stereographic back to lon/lat, about 11 km in latitude
    dict(name="polar-to-4326", source="polar-f32-zstd-bigtiff.tif",
         crs="EPSG:4326", bounds=(-180.0, -90.0, 180.0, -62.0), size=(1800, 140)),
    # UTM 55S, three bands, to lon/lat at about 10 m
    dict(name="utm-to-4326", source="utm-u16-lzw-pred2-3band.tif",
         crs="EPSG:4326", bounds=(147.000, -43.440, 147.036, -43.424), size=(290, 180)),
]
METHODS = {"near": Resampling.nearest, "bilinear": Resampling.bilinear, "cubic": Resampling.cubic}


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__,
           "error_threshold": 0.125, "cases": []}
    for c in CASES:
        w, h = c["size"]
        dst_transform = from_bounds(*c["bounds"], w, h)
        with rasterio.open(os.path.join(COG, c["source"])) as src:
            count, dtype, nodata = src.count, src.dtypes[0], src.nodata
            rec = {
                "name": c["name"],
                "source": "cog/" + c["source"],
                "crs": c["crs"],
                "geotransform": list(dst_transform.to_gdal()),
                "width": w,
                "height": h,
                "count": count,
                "dtype": dtype,
                "nodata": nodata,
                "results": [],
            }
            for m, alg in METHODS.items():
                dst = np.full((count, h, w), nodata if nodata is not None else 0, dtype=dtype)
                for b in range(count):
                    reproject(
                        source=rasterio.band(src, b + 1),
                        destination=dst[b],
                        dst_transform=dst_transform,
                        dst_crs=c["crs"],
                        resampling=alg,
                        src_nodata=nodata,
                        dst_nodata=nodata,
                        num_threads=1,
                    )
                name = "%s-%s.bin" % (c["name"], m)
                dst.astype(dst.dtype.newbyteorder("<")).tofile(os.path.join(OUT, name))
                rec["results"].append({"resampling": m, "file": name})
        out["cases"].append(rec)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
