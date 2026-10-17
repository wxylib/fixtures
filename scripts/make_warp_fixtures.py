"""Warp oracles: what GDAL's warper (through rasterio.warp.reproject, with
the approximate transformer at its default 0.125 pixel error threshold)
writes for a few source fixtures, target grids and resampling methods.

Some target grids are near 1:1 with the source resolution; others are
coarser, where GDAL stretches its interpolation kernels (wxylib's scaled
warp, stage S5), or need the source window split into chunks.

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
    # a 0..360 source (made below) to south polar stereographic: the 0/360
    # seam runs through the target, where the approximate transformer must
    # see the longitude wrap's jump
    dict(name="world360-to-3031", source="world360-int16.tif", made=True,
         crs="EPSG:3031", bounds=(-4.0e6, -4.0e6, 4.0e6, 4.0e6), size=(400, 400)),
    # the world at 1 degree with the ocean as nodata (int16), to south polar
    # stereographic: kernels that skip invalid cells, and target values
    # stepped off nodata
    dict(name="worldnd-to-3031", source="worldnd-int16.tif", made=True,
         crs="EPSG:3031", bounds=(-4.0e6, -4.0e6, 4.0e6, 4.0e6), size=(200, 200)),
    # the polar fixture to the whole globe: most of the target misses the
    # source, so GDAL splits the target into chunks with their own windows
    dict(name="polar-to-4326-global", source="polar-f32-zstd-bigtiff.tif",
         crs="EPSG:4326", bounds=(-180.0, -90.0, 180.0, 90.0), size=(720, 360)),
    # UTM 55S, three bands, to lon/lat at about 10 m
    dict(name="utm-to-4326", source="utm-u16-lzw-pred2-3band.tif",
         crs="EPSG:4326", bounds=(147.000, -43.440, 147.036, -43.424), size=(290, 180)),
    # the same at about 50 m: uint16 kernels stretched over 5 by 4 cells
    dict(name="utm-to-4326-coarse", source="utm-u16-lzw-pred2-3band.tif",
         crs="EPSG:4326", bounds=(147.000, -43.440, 147.036, -43.424), size=(58, 45)),
    # Australian Albers, three uint8 bands, to EPSG:3857 at about 3 cells to 1
    dict(name="rgb-to-3857-coarse", source="rgb-u8-deflate-pixel.tif",
         crs="EPSG:3857", bounds=(14695200.0, 2300.0, 14697100.0, 3900.0), size=(63, 53)),
]
METHODS = {"near": Resampling.nearest, "bilinear": Resampling.bilinear, "cubic": Resampling.cubic}


def make_world360(path):
    """The world fixture at 1 degree, columns rolled to start at lon 0."""
    with rasterio.open(os.path.join(COG, "world-int16-deflate.tif")) as src:
        a = src.read(1, out_shape=(180, 360), resampling=Resampling.average)
    a = np.roll(a, -180, axis=1)
    profile = dict(driver="GTiff", width=360, height=180, count=1, dtype="int16",
                   crs="EPSG:4326", transform=from_bounds(0, -90, 360, 90, 360, 180),
                   tiled=True, blockxsize=256, blockysize=256, compress="deflate")
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(a, 1)


def make_worldnd(path):
    """The world fixture at 1 degree, values below 0 set to nodata."""
    with rasterio.open(os.path.join(COG, "world-int16-deflate.tif")) as src:
        a = src.read(1, out_shape=(180, 360), resampling=Resampling.average)
    a[a < 0] = -32768
    profile = dict(driver="GTiff", width=360, height=180, count=1, dtype="int16",
                   crs="EPSG:4326", transform=from_bounds(-180, -90, 180, 90, 360, 180),
                   tiled=True, blockxsize=256, blockysize=256, compress="deflate",
                   nodata=-32768)
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(a, 1)


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    make_world360(os.path.join(OUT, "world360-int16.tif"))
    make_worldnd(os.path.join(OUT, "worldnd-int16.tif"))
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__,
           "error_threshold": 0.125, "cases": []}
    for c in CASES:
        w, h = c["size"]
        dst_transform = from_bounds(*c["bounds"], w, h)
        folder, prefix = (OUT, "warp/") if c.get("made") else (COG, "cog/")
        with rasterio.open(os.path.join(folder, c["source"])) as src:
            count, dtype, nodata = src.count, src.dtypes[0], src.nodata
            rec = {
                "name": c["name"],
                "source": prefix + c["source"],
                "crs": c["crs"],
                "geotransform": list(dst_transform.to_gdal()),
                "width": w,
                "height": h,
                "count": count,
                "dtype": dtype,
                "nodata": nodata,
                # GDAL wraps longitude for a source that starts at 0 (the
                # whole western hemisphere would be empty otherwise), so
                # wxylib compares this case with its own wrap on
                "wrap_globe": c["source"].startswith("world360"),
                "results": [],
            }
            for m, alg in METHODS.items():
                if m not in c.get("methods", METHODS):
                    continue
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
