"""Point-read oracles (wxylib stage S4c): random points over each COG
fixture, given in EPSG:4326 and in the file's own CRS, and the values GDAL
reads at them as gdallocationinfo does: reproject the point to the file's
CRS, apply the inverse geotransform, take the cell the point falls in
(floor), and read every band there. Points off the file have no value.

Each case writes little-endian float64 x and y (`*-x.bin`, `*-y.bin`),
the values as float64 band by band (`*-values.bin`, shaped bands x n) and
a uint8 mask, 1 where the point is on the file (`*-found.bin`).

Run from the repository root:

    python scripts/make_point_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import rasterio
from rasterio.warp import transform, transform_bounds

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "points")
N = 400


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    rng = np.random.default_rng(20261017)
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__, "cases": []}
    for f in sorted(os.listdir(os.path.join(ROOT, "cog"))):
        if not f.endswith(".tif"):
            continue
        path = os.path.join(ROOT, "cog", f)
        with rasterio.open(path) as ds:
            data = ds.read()
            native = ds.crs.to_string()
            b = ds.bounds
            for crs in dict.fromkeys(["EPSG:4326", native]):
                if crs == native:
                    box = (b.left, b.bottom, b.right, b.top)
                else:
                    box = transform_bounds(ds.crs, crs, *b, densify_pts=21)
                # 10 percent beyond the box on each side, so some points miss
                w, h = box[2] - box[0], box[3] - box[1]
                x = rng.uniform(box[0] - 0.1 * w, box[2] + 0.1 * w, N)
                y = rng.uniform(box[1] - 0.1 * h, box[3] + 0.1 * h, N)
                if crs == "EPSG:4326":
                    x = np.clip(x, -180, 180)
                    y = np.clip(y, -90, 90)
                if crs == native:
                    nx, ny = x, y
                else:
                    nx, ny = (np.array(v) for v in transform(crs, ds.crs, x, y))
                inv = ~ds.transform
                col, row = inv * (nx, ny)
                col, row = np.floor(col), np.floor(row)
                found = (np.isfinite(col) & np.isfinite(row) & (col >= 0) & (row >= 0)
                         & (col < ds.width) & (row < ds.height))
                values = np.zeros((ds.count, N))
                ci, ri = col[found].astype(int), row[found].astype(int)
                values[:, found] = data[:, ri, ci]
                name = "%s-%s" % (f[:-4], "4326" if crs == "EPSG:4326" else "native")
                x.astype("<f8").tofile(os.path.join(OUT, name + "-x.bin"))
                y.astype("<f8").tofile(os.path.join(OUT, name + "-y.bin"))
                values.astype("<f8").tofile(os.path.join(OUT, name + "-values.bin"))
                found.astype("u1").tofile(os.path.join(OUT, name + "-found.bin"))
                out["cases"].append({"name": name, "source": "cog/" + f, "crs": crs, "n": N,
                                     "bands": ds.count, "found": int(found.sum())})
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
