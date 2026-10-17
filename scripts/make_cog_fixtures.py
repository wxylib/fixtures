# Writes small COG fixtures with real (synthetic) pixel values, and records
# what GDAL (through its own Python bindings, osgeo.gdal) reads for a few
# windows, so wxylib readers are checked against GDAL. Run from the
# repository root, with a Python that has osgeo.gdal and numpy (for example
# a conda-forge env: micromamba create -n gdal -c conda-forge gdal numpy):
#
#   python scripts/make_cog_fixtures.py
#
# Output: cog/*.tif, cog/oracle/*.bin (little-endian, C order, band-major)
# and cog/manifest.json. Content is ASCII only.

import json, os
import numpy as np
from osgeo import gdal, gdal_array, osr

gdal.UseExceptions()

OUT = os.path.join(os.path.dirname(__file__), "..", "cog")
ORACLE = os.path.join(OUT, "oracle")
os.makedirs(ORACLE, exist_ok=True)

def field(h, w, bands, dtype, seed):
    y, x = np.mgrid[0:h, 0:w].astype(np.float64)
    out = []
    for b in range(bands):
        v = 1000 * np.sin(x / (37 + 11 * b) + seed) * np.cos(y / (23 + 7 * b)) + 7 * x - 3 * y + 100 * b
        out.append(v)
    a = np.stack(out)
    if np.issubdtype(np.dtype(dtype), np.integer):
        info = np.iinfo(dtype)
        a = np.clip(np.round(a), info.min, info.max)
    return a.astype(dtype)

def write(path, data, bounds, crs, nodata, driver, options):
    """Write all bands in one call. GTiff is created in place; the COG driver
    can only CreateCopy, so it copies from an in-memory dataset."""
    bands, h, w = data.shape
    etype = gdal_array.NumericTypeCodeToGDALTypeCode(data.dtype)
    if driver == "COG":
        ds = gdal.GetDriverByName("MEM").Create("", w, h, bands, etype)
    else:
        ds = gdal.GetDriverByName(driver).Create(path, w, h, bands, etype, options=options)
    xmin, ymin, xmax, ymax = bounds
    ds.SetGeoTransform([xmin, (xmax - xmin) / w, 0.0, ymax, 0.0, (ymin - ymax) / h])
    srs = osr.SpatialReference()
    srs.SetFromUserInput(crs)
    ds.SetProjection(srs.ExportToWkt())
    if nodata is not None:
        for b in range(bands):
            ds.GetRasterBand(b + 1).SetNoDataValue(nodata)
    ds.WriteArray(data)
    if driver == "COG":
        gdal.GetDriverByName(driver).CreateCopy(path, ds, options=options)
    ds = None

def overview_factors(ds):
    b = ds.GetRasterBand(1)
    return [int(round(ds.RasterXSize / b.GetOverview(i).XSize)) for i in range(b.GetOverviewCount())]

cases = [
    dict(file="world-int16-deflate.tif", w=720, h=360, bands=1, dtype="int16",
         bounds=(-180, -90, 180, 90), crs="EPSG:4326",
         driver="COG", options=["BLOCKSIZE=128", "COMPRESS=DEFLATE"]),
    dict(file="polar-f32-zstd-bigtiff.tif", w=500, h=400, bands=1, dtype="float32",
         bounds=(-3e6, -2.4e6, 3e6, 2.4e6), crs="EPSG:3031", nodata=-9999.0,
         driver="COG", options=["BLOCKSIZE=128", "COMPRESS=ZSTD", "BIGTIFF=YES"]),
    dict(file="utm-u16-lzw-pred2-3band.tif", w=300, h=200, bands=3, dtype="uint16",
         bounds=(499980, 5190220, 502980, 5192220), crs="EPSG:32755",
         driver="COG", options=["BLOCKSIZE=64", "COMPRESS=LZW", "PREDICTOR=2"]),
    dict(file="planar-u8-deflate-2band.tif", w=200, h=150, bands=2, dtype="uint8",
         bounds=(0, 0, 200, 150), crs="EPSG:3857",
         driver="GTiff", options=["TILED=YES", "BLOCKXSIZE=64", "BLOCKYSIZE=64",
                                  "COMPRESS=DEFLATE", "INTERLEAVE=BAND"]),
    dict(file="rgb-u8-deflate-pixel.tif", w=230, h=170, bands=3, dtype="uint8",
         bounds=(1000, 2000, 3300, 3700), crs="EPSG:3577",
         driver="GTiff", options=["TILED=YES", "BLOCKXSIZE=64", "BLOCKYSIZE=64",
                                  "COMPRESS=DEFLATE", "INTERLEAVE=PIXEL"]),
]

windows = [(0, 0, 50, 40), (97, 61, 130, 77), (None, None, None, None)]

manifest = {"gdal": gdal.__version__, "cases": []}
for i, c in enumerate(cases):
    path = os.path.join(OUT, c["file"])
    data = field(c["h"], c["w"], c["bands"], c["dtype"], seed=i)
    if "nodata" in c:
        data[:, :10, :10] = c["nodata"]
    write(path, data, c["bounds"], c["crs"], c.get("nodata"), c["driver"], c["options"])
    rec = {"file": c["file"], "windows": []}
    ds = gdal.Open(path)
    b1 = ds.GetRasterBand(1)
    srs = ds.GetSpatialRef()
    srs.AutoIdentifyEPSG()
    bx, by = b1.GetBlockSize()
    factors = overview_factors(ds)
    rec.update(width=ds.RasterXSize, height=ds.RasterYSize, count=ds.RasterCount,
               dtype=np.dtype(gdal_array.GDALTypeCodeToNumericTypeCode(b1.DataType)).name,
               block=[by, bx], overviews=factors, geotransform=list(ds.GetGeoTransform()),
               epsg=int(srs.GetAuthorityCode(None)), nodata=b1.GetNoDataValue(),
               interleave=ds.GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE"))
    for k, (col, row, w, h) in enumerate(windows):
        if col is None:  # last window: the bottom-right corner, crossing tile edges
            w, h = 70, 45
            col, row = ds.RasterXSize - w, ds.RasterYSize - h
        # keep every window inside the image
        w, h = min(w, ds.RasterXSize - col), min(h, ds.RasterYSize - row)
        arr = ds.ReadAsArray(col, row, w, h)
        name = f"{os.path.splitext(c['file'])[0]}-w{k}.bin"
        arr.astype(arr.dtype.newbyteorder("<")).tofile(os.path.join(ORACLE, name))
        rec["windows"].append({"col": col, "row": row, "width": w, "height": h, "file": "oracle/" + name})
    for lvl, f in enumerate(factors, start=1):
        ov = gdal.OpenEx(path, open_options=[f"OVERVIEW_LEVEL={lvl - 1}"])
        rec.setdefault("levels", []).append({"factor": f, "width": ov.RasterXSize, "height": ov.RasterYSize,
                                             "geotransform": list(ov.GetGeoTransform())})
        ov = None
    ds = None
    manifest["cases"].append(rec)

with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=1)
print(json.dumps([(r["file"], r["width"], r["height"], r["block"], r["overviews"]) for r in manifest["cases"]]))
