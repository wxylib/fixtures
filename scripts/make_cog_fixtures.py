# Writes small COG fixtures with real (synthetic) pixel values, and records
# what GDAL (through rasterio) reads for a few windows, so wxylib readers are
# checked against GDAL. Run from the repository root:
#
#   python3 -m pip install rasterio numpy
#   python3 scripts/make_cog_fixtures.py
#
# Output: cog/*.tif, cog/oracle/*.bin (little-endian, C order, band-major)
# and cog/manifest.json. Content is ASCII only.

import json, os
import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.windows import Window

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

cases = [
    dict(file="world-int16-deflate.tif", w=720, h=360, bands=1, dtype="int16",
         bounds=(-180, -90, 180, 90), crs="EPSG:4326",
         profile=dict(driver="COG", blocksize=128, compress="DEFLATE")),
    dict(file="polar-f32-zstd-bigtiff.tif", w=500, h=400, bands=1, dtype="float32",
         bounds=(-3e6, -2.4e6, 3e6, 2.4e6), crs="EPSG:3031", nodata=-9999.0,
         profile=dict(driver="COG", blocksize=128, compress="ZSTD", BIGTIFF="YES")),
    dict(file="utm-u16-lzw-pred2-3band.tif", w=300, h=200, bands=3, dtype="uint16",
         bounds=(499980, 5190220, 502980, 5192220), crs="EPSG:32755",
         profile=dict(driver="COG", blocksize=64, compress="LZW", predictor=2)),
    dict(file="planar-u8-deflate-2band.tif", w=200, h=150, bands=2, dtype="uint8",
         bounds=(0, 0, 200, 150), crs="EPSG:3857",
         profile=dict(driver="GTiff", tiled=True, blockxsize=64, blockysize=64,
                      compress="DEFLATE", interleave="band")),
]

windows = [(0, 0, 50, 40), (97, 61, 130, 77), (None, None, None, None)]

manifest = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__, "cases": []}
for i, c in enumerate(cases):
    path = os.path.join(OUT, c["file"])
    data = field(c["h"], c["w"], c["bands"], c["dtype"], seed=i)
    if "nodata" in c:
        data[:, :10, :10] = c["nodata"]
    with rasterio.open(path, "w", width=c["w"], height=c["h"], count=c["bands"], dtype=c["dtype"],
                       crs=c["crs"], transform=from_bounds(*c["bounds"], c["w"], c["h"]),
                       nodata=c.get("nodata"), **c["profile"]) as ds:
        ds.write(data)
    rec = {"file": c["file"], "windows": []}
    with rasterio.open(path) as ds:
        t = ds.transform
        rec.update(width=ds.width, height=ds.height, count=ds.count, dtype=ds.dtypes[0],
                   block=list(ds.block_shapes[0]), overviews=ds.overviews(1),
                   geotransform=[t.c, t.a, t.b, t.f, t.d, t.e], epsg=ds.crs.to_epsg(),
                   nodata=ds.nodata, interleave=ds.interleaving.value if ds.interleaving else None)
        for k, (col, row, w, h) in enumerate(windows):
            if col is None:  # last window: the bottom-right corner, crossing tile edges
                w, h = 70, 45
                col, row = ds.width - w, ds.height - h
            arr = ds.read(window=Window(col, row, w, h))
            name = f"{os.path.splitext(c['file'])[0]}-w{k}.bin"
            arr.astype(arr.dtype.newbyteorder("<")).tofile(os.path.join(ORACLE, name))
            rec["windows"].append({"col": col, "row": row, "width": w, "height": h, "file": "oracle/" + name})
        for lvl, f in enumerate(ds.overviews(1), start=1):
            with rasterio.open(path, overview_level=lvl - 1) as ov:
                t = ov.transform
                rec.setdefault("levels", []).append({"factor": f, "width": ov.width, "height": ov.height,
                                                     "geotransform": [t.c, t.a, t.b, t.f, t.d, t.e]})
    manifest["cases"].append(rec)

with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=1)
print(json.dumps([(r["file"], r["width"], r["height"], r["block"], r["overviews"]) for r in manifest["cases"]]))
