"""Point-read oracles (wxylib stage S4c): random points over each COG
fixture, given in EPSG:4326 and in the file's own CRS, and the values GDAL
reads at them as gdallocationinfo does: reproject the point to the file's
CRS, apply the inverse geotransform, take the cell the point falls in
(floor), and read every band there. Points off the file have no value.

Each case writes little-endian float64 x and y (`*-x.bin`, `*-y.bin`),
the values as float64 band by band (`*-values.bin`, shaped bands x n) and
a uint8 mask, 1 where the point is on the file (`*-found.bin`).

Interpolated reads come from GDAL's own GDALRasterInterpolateAtGeolocation
(GDAL 3.11+), called through ctypes in the GDAL that rasterio bundles:
`*-bilinear.bin`, `*-cubic.bin` and `*-cubicspline.bin`, float64 bands x
n, NaN where GDAL reports no value (off the file, or nodata in the
kernel). GDAL's own (pixel, line) for each point, from
GDALDatasetGeolocationToPixelLine, is in `*-pixel.bin` and `*-line.bin`.

Run from the repository root:

    python scripts/make_point_fixtures.py
"""

import ctypes
import glob
import json
import os
import shutil

import numpy as np
import rasterio
from rasterio.warp import transform, transform_bounds

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "points")
N = 400
INTERP = {"bilinear": 1, "cubic": 2, "cubicspline": 3}


def gdal_lib():
    libs = os.path.join(os.path.dirname(rasterio.__file__), os.pardir, "rasterio.libs")
    g = ctypes.CDLL(glob.glob(os.path.join(libs, "libgdal*.so*"))[0])
    g.GDALAllRegister()
    vp, d = ctypes.c_void_p, ctypes.c_double
    g.GDALOpen.restype = vp
    g.GDALOpen.argtypes = [ctypes.c_char_p, ctypes.c_int]
    g.GDALClose.argtypes = [vp]
    g.GDALGetRasterBand.restype = vp
    g.GDALGetRasterBand.argtypes = [vp, ctypes.c_int]
    g.OSRNewSpatialReference.restype = vp
    g.OSRNewSpatialReference.argtypes = [ctypes.c_char_p]
    g.OSRSetFromUserInput.argtypes = [vp, ctypes.c_char_p]
    g.OSRSetAxisMappingStrategy.argtypes = [vp, ctypes.c_int]
    p = g.GDALDatasetGeolocationToPixelLine
    p.restype = ctypes.c_int
    p.argtypes = [vp, d, d, vp, ctypes.POINTER(d), ctypes.POINTER(d), vp]
    f = g.GDALRasterInterpolateAtGeolocation
    f.restype = ctypes.c_int
    f.argtypes = [vp, d, d, vp, ctypes.c_int, ctypes.POINTER(d), ctypes.POINTER(d), vp]
    return g


def srs_for(g, crs, native):
    if crs == native:
        return None
    srs = g.OSRNewSpatialReference(None)
    g.OSRSetFromUserInput(srs, crs.encode())
    g.OSRSetAxisMappingStrategy(srs, 0)  # OAMS_TRADITIONAL_GIS_ORDER
    return srs


def pixel_line(g, path, crs, native, x, y):
    """GDAL's (pixel, line) for each point, NaN where it has none."""
    ds = g.GDALOpen(path.encode(), 0)
    srs = srs_for(g, crs, native)
    out = np.full((2, len(x)), np.nan)
    px, ln = ctypes.c_double(), ctypes.c_double()
    for i in range(len(x)):
        if g.GDALDatasetGeolocationToPixelLine(
                ds, x[i], y[i], srs, ctypes.byref(px), ctypes.byref(ln), None) == 0:
            out[:, i] = px.value, ln.value
    g.GDALClose(ds)
    return out


def interpolate(g, path, crs, native, x, y, bands, alg):
    """GDAL's interpolated read at each point, NaN where it has none."""
    ds = g.GDALOpen(path.encode(), 0)
    srs = srs_for(g, crs, native)
    out = np.full((bands, len(x)), np.nan)
    re, im = ctypes.c_double(), ctypes.c_double()
    for b in range(bands):
        band = g.GDALGetRasterBand(ds, b + 1)
        for i in range(len(x)):
            err = g.GDALRasterInterpolateAtGeolocation(
                band, x[i], y[i], srs, alg, ctypes.byref(re), ctypes.byref(im), None)
            if err == 0:
                out[b, i] = re.value
    g.GDALClose(ds)
    return out


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    rng = np.random.default_rng(20261017)
    g = gdal_lib()
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__,
           "interpolated": list(INTERP), "cases": []}
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
                pl = pixel_line(g, path, crs, native, x, y)
                pl[0].astype("<f8").tofile(os.path.join(OUT, name + "-pixel.bin"))
                pl[1].astype("<f8").tofile(os.path.join(OUT, name + "-line.bin"))
                for m, alg in INTERP.items():
                    v = interpolate(g, path, crs, native, x, y, ds.count, alg)
                    v.astype("<f8").tofile(os.path.join(OUT, "%s-%s.bin" % (name, m)))
                out["cases"].append({"name": name, "source": "cog/" + f, "crs": crs, "n": N,
                                     "bands": ds.count, "found": int(found.sum())})
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
