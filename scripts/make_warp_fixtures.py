"""Warp oracles: what GDAL's warper (through osgeo.gdal's gdal.Warp, with
the approximate transformer at its default 0.125 pixel error threshold)
writes for a few source fixtures, target grids and resampling methods.

Some target grids are near 1:1 with the source resolution; others are
coarser, where GDAL stretches its interpolation kernels (wxylib's scaled
warp, stage S5), or need the source window split into chunks.

Run from the repository root after make_cog_fixtures.py, with a Python
that has GDAL's bindings and numpy (e.g. a conda-forge env):

    python scripts/make_warp_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import osgeo
from osgeo import gdal, gdal_array, osr

gdal.UseExceptions()

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
    # gdalwarp's default overview choice (-ovr AUTO): the coarsest
    # overview still finer than the target, then the warp from it
    dict(name="utm-to-4326-ovr", source="utm-u16-lzw-pred2-3band.tif", overview="auto",
         crs="EPSG:4326", bounds=(147.000, -43.440, 147.036, -43.424), size=(58, 45)),
    dict(name="world-to-4326-ovr", source="world-int16-deflate.tif", overview="auto",
         crs="EPSG:4326", bounds=(-180.0, -90.0, 180.0, 90.0), size=(100, 50)),
    dict(name="polar-to-4326-ovr", source="polar-f32-zstd-bigtiff.tif", overview="auto",
         crs="EPSG:4326", bounds=(-180.0, -90.0, 180.0, -60.0), size=(90, 15)),
    # lon/lat to web mercator around (0, 0): GDAL's special points for a
    # geographic source pass lon 0, lat +-89.9999 through the target
    # geotransform without reprojection, so here they land in the target
    # and switch source windows to grid sampling
    dict(name="world-to-3857-origin", source="world-int16-deflate.tif",
         crs="EPSG:3857", bounds=(-2.0e6, -2.0e6, 2.0e6, 2.0e6), size=(37, 37)),
    # Australian Albers, three uint8 bands, to EPSG:3857 at about 3 cells to 1
    dict(name="rgb-to-3857-coarse", source="rgb-u8-deflate-pixel.tif",
         crs="EPSG:3857", bounds=(14695200.0, 2300.0, 14697100.0, 3900.0), size=(63, 53)),
    # the world to south polar stereographic under a 1 MB warp memory limit
    # (gdalwarp -wm 1): GDAL's planner halves the target until each chunk's
    # working memory fits, so the warp runs as several chunks
    dict(name="world-to-3031-wm1", source="world-int16-deflate.tif", wm=1,
         methods=("near", "cubic"),
         crs="EPSG:3031", bounds=(-4.0e6, -4.0e6, 4.0e6, 4.0e6), size=(1200, 1200)),
]
# lanczos, average and mode are rwarp's kernels on each GDAL chunk in
# wxylib (hypertidy/rwarp#3, wxylib/core#19)
METHODS = ("near", "bilinear", "cubic", "lanczos", "average", "mode")


def geotransform(bounds, w, h):
    xmin, ymin, xmax, ymax = bounds
    return [xmin, (xmax - xmin) / w, 0.0, ymax, 0.0, -(ymax - ymin) / h]


def srs_wkt(crs):
    srs = osr.SpatialReference()
    srs.SetFromUserInput(crs)
    return srs.ExportToWkt()


def world_1deg():
    """The world fixture averaged to 1 degree (360 by 180)."""
    ds = gdal.Open(os.path.join(COG, "world-int16-deflate.tif"))
    return ds.GetRasterBand(1).ReadAsArray(
        buf_xsize=360, buf_ysize=180, resample_alg=gdal.GRIORA_Average)


def write_1deg(path, a, xmin, nodata=None):
    ds = gdal.GetDriverByName("GTiff").Create(
        path, 360, 180, 1, gdal.GDT_Int16,
        ["TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "COMPRESS=DEFLATE"])
    ds.SetGeoTransform([xmin, 1.0, 0.0, 90.0, 0.0, -1.0])
    ds.SetProjection(srs_wkt("EPSG:4326"))
    band = ds.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    band.WriteArray(a)
    ds = None


def make_world360(path):
    """The world fixture at 1 degree, columns rolled to start at lon 0."""
    write_1deg(path, np.roll(world_1deg(), -180, axis=1), 0.0)


def make_worldnd(path):
    """The world fixture at 1 degree, values below 0 set to nodata."""
    a = world_1deg()
    a[a < 0] = -32768
    write_1deg(path, a, -180.0, nodata=-32768)


class OverviewLog:
    """Collects the overview level gdal.Warp reports choosing (its WARP
    debug message "Selecting overview level N for ...")."""

    def __init__(self):
        self.level = None

    def __call__(self, cls, num, msg):
        if cls == gdal.CE_Debug and msg.startswith("WARP: Selecting overview level"):
            self.level = int(msg.split()[4])


def warp(src, c, w, h, gt, method, overview):
    """One band-by-band warp to a MEM dataset of the target grid filled
    with nodata (or 0) first: GDAL skips target chunks that miss the
    source, so INIT_DEST alone leaves those at 0. Returns
    the (count, h, w) array and the overview level gdal.Warp used."""
    count = src.RasterCount
    band1 = src.GetRasterBand(1)
    nodata = band1.GetNoDataValue()
    dst = gdal.GetDriverByName("MEM").Create("", w, h, count, band1.DataType)
    dst.SetGeoTransform(gt)
    dst.SetProjection(srs_wkt(c["crs"]))
    log = OverviewLog()
    gdal.PushErrorHandler(log)
    gdal.SetConfigOption("CPL_DEBUG", "WARP")
    try:
        for b in range(1, count + 1):
            band = dst.GetRasterBand(b)
            if nodata is not None:
                band.SetNoDataValue(nodata)
            band.Fill(nodata if nodata is not None else 0)
            opts = gdal.WarpOptions(
                srcBands=[b], dstBands=[b],
                resampleAlg=method,
                srcNodata=nodata, dstNodata=nodata,
                errorThreshold=0.125,
                warpMemoryLimit=c.get("wm", 64),
                multithread=False,
                overviewLevel="AUTO" if overview else "NONE",
                warpOptions=["INIT_DEST=NO_DATA"])
            gdal.Warp(dst, src, options=opts)
    finally:
        gdal.SetConfigOption("CPL_DEBUG", None)
        gdal.PopErrorHandler()
    return dst.ReadAsArray().reshape(count, h, w), log.level


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    make_world360(os.path.join(OUT, "world360-int16.tif"))
    make_worldnd(os.path.join(OUT, "worldnd-int16.tif"))
    out = {"gdal": gdal.VersionInfo("RELEASE_NAME"), "osgeo": osgeo.__version__,
           "error_threshold": 0.125, "cases": []}
    for c in CASES:
        w, h = c["size"]
        gt = geotransform(c["bounds"], w, h)
        folder, prefix = (OUT, "warp/") if c.get("made") else (COG, "cog/")
        src = gdal.Open(os.path.join(folder, c["source"]))
        band1 = src.GetRasterBand(1)
        nodata = band1.GetNoDataValue()
        rec = {
            "name": c["name"],
            "source": prefix + c["source"],
            "crs": c["crs"],
            "geotransform": gt,
            "width": w,
            "height": h,
            "count": src.RasterCount,
            "dtype": np.dtype(gdal_array.GDALTypeCodeToNumericTypeCode(band1.DataType)).name,
            "nodata": nodata,
            # GDAL wraps longitude for a source that starts at 0 (the
            # whole western hemisphere would be empty otherwise), so
            # wxylib compares this case with its own wrap on
            "wrap_globe": c["source"].startswith("world360"),
            # gdalwarp's -ovr: "auto" cases warp from the overview its
            # default choice picks (0 is the first overview)
            "overview": c.get("overview", "none"),
            "overview_level": None,
            # gdalwarp's -wm, in MB (GDAL's default is 64)
            "warp_memory_mb": c.get("wm", 64),
            "results": [],
        }
        for m in METHODS:
            if m not in c.get("methods", METHODS):
                continue
            dst, level = warp(src, c, w, h, gt, m, c.get("overview") == "auto")
            rec["overview_level"] = level
            name = "%s-%s.bin" % (c["name"], m)
            dst.astype(dst.dtype.newbyteorder("<")).tofile(os.path.join(OUT, name))
            rec["results"].append({"resampling": m, "file": name})
        out["cases"].append(rec)
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
