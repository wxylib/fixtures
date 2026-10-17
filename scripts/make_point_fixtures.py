"""Point-read oracles (wxylib stage S4c): random points over each COG
fixture, given in EPSG:4326 and in the file's own CRS, and the values GDAL
reads at them as gdallocationinfo does: reproject the point to the file's
CRS, apply the inverse geotransform, take the cell the point falls in
(floor), and read every band there. Points off the file have no value.

Each case writes little-endian float64 x and y (`*-x.bin`, `*-y.bin`),
the values as float64 band by band (`*-values.bin`, shaped bands x n) and
a uint8 mask, 1 where the point is on the file (`*-found.bin`).

Interpolated reads come from GDAL's own Band.InterpolateAtGeolocation
(GDALRasterInterpolateAtGeolocation, GDAL 3.11+): `*-bilinear.bin`,
`*-cubic.bin` and `*-cubicspline.bin`, float64 bands x n, NaN where GDAL
reports no value (off the file, or nodata in the kernel). GDAL's own
(pixel, line) for each point, as GDALDatasetGeolocationToPixelLine gives
it, is in `*-pixel.bin` and `*-line.bin`.

Uses GDAL's own Python bindings (osgeo.gdal, osr) and numpy, from a
conda-forge env (GDAL 3.11+). Run from the repository root:

    python scripts/make_point_fixtures.py
"""

import json
import os
import shutil

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "points")
N = 400
INTERP = {"bilinear": gdal.GRIORA_Bilinear, "cubic": gdal.GRIORA_Cubic,
          "cubicspline": gdal.GRIORA_CubicSpline}


def srs_for(crs):
    srs = osr.SpatialReference()
    srs.SetFromUserInput(crs)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def crs_string(srs):
    """EPSG:code where the file's CRS has one, else WKT."""
    if srs.GetAuthorityName(None) and srs.GetAuthorityCode(None):
        return "%s:%s" % (srs.GetAuthorityName(None), srs.GetAuthorityCode(None))
    return srs.ExportToWkt()


def pixel_line(ds, srs, x, y):
    """GDAL's (pixel, line) for each point, NaN where it has none, as
    GDALDatasetGeolocationToPixelLine does it: a GenImgProj transformer
    from the dataset to `srs` (None: the dataset's own CRS; lon, lat order
    for a geographic DST_SRS), inverted."""
    opts = [] if srs is None else ["DST_SRS=" + srs.ExportToWkt(["FORMAT=WKT2"])]
    tr = gdal.Transformer(ds, None, opts)
    out = np.full((2, len(x)), np.nan)
    for i in range(len(x)):
        ok, (px, ln, _) = tr.TransformPoint(1, float(x[i]), float(y[i]))
        if ok:
            out[:, i] = px, ln
    return out


def interpolate(ds, srs, x, y, alg):
    """GDAL's interpolated read at each point, NaN where it has none
    (InterpolateAtGeolocation returns None)."""
    out = np.full((ds.RasterCount, len(x)), np.nan)
    # GDAL 3.12 warns on every call that its own transformer ignores
    # DST_SRS_AXIS_MAPPING_STRATEGY; harmless, so keep it quiet
    with gdal.quiet_errors():
        for b in range(ds.RasterCount):
            band = ds.GetRasterBand(b + 1)
            for i in range(len(x)):
                v = band.InterpolateAtGeolocation(float(x[i]), float(y[i]), srs, alg)
                if v is not None:
                    out[b, i] = v
    return out


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    rng = np.random.default_rng(20261017)
    out = {"gdal": gdal.__version__, "interpolated": list(INTERP), "cases": []}
    for f in sorted(os.listdir(os.path.join(ROOT, "cog"))):
        if not f.endswith(".tif"):
            continue
        path = os.path.join(ROOT, "cog", f)
        ds = gdal.Open(path)
        nb, width, height = ds.RasterCount, ds.RasterXSize, ds.RasterYSize
        data = ds.ReadAsArray().reshape(nb, height, width)
        dsrs = ds.GetSpatialRef()
        dsrs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        native = crs_string(dsrs)
        gt = ds.GetGeoTransform()
        b = (gt[0], gt[3] + gt[5] * height, gt[0] + gt[1] * width, gt[3])
        for crs in dict.fromkeys(["EPSG:4326", native]):
            srs = None if crs == native else srs_for(crs)
            if srs is None:
                box = b
            else:
                box = osr.CoordinateTransformation(dsrs, srs).TransformBounds(*b, 21)
            # 10 percent beyond the box on each side, so some points miss
            w, h = box[2] - box[0], box[3] - box[1]
            x = rng.uniform(box[0] - 0.1 * w, box[2] + 0.1 * w, N)
            y = rng.uniform(box[1] - 0.1 * h, box[3] + 0.1 * h, N)
            if crs == "EPSG:4326":
                x = np.clip(x, -180, 180)
                y = np.clip(y, -90, 90)
            if srs is None:
                nx, ny = x, y
            else:
                pts = osr.CoordinateTransformation(srs, dsrs).TransformPoints(
                    np.column_stack([x, y]))
                nx, ny = np.array(pts)[:, 0], np.array(pts)[:, 1]
            inv = gdal.InvGeoTransform(gt)
            col = np.floor(inv[0] + nx * inv[1] + ny * inv[2])
            row = np.floor(inv[3] + nx * inv[4] + ny * inv[5])
            found = (np.isfinite(col) & np.isfinite(row) & (col >= 0) & (row >= 0)
                     & (col < width) & (row < height))
            values = np.zeros((nb, N))
            ci, ri = col[found].astype(int), row[found].astype(int)
            values[:, found] = data[:, ri, ci]
            name = "%s-%s" % (f[:-4], "4326" if crs == "EPSG:4326" else "native")
            x.astype("<f8").tofile(os.path.join(OUT, name + "-x.bin"))
            y.astype("<f8").tofile(os.path.join(OUT, name + "-y.bin"))
            values.astype("<f8").tofile(os.path.join(OUT, name + "-values.bin"))
            found.astype("u1").tofile(os.path.join(OUT, name + "-found.bin"))
            pl = pixel_line(ds, srs, x, y)
            pl[0].astype("<f8").tofile(os.path.join(OUT, name + "-pixel.bin"))
            pl[1].astype("<f8").tofile(os.path.join(OUT, name + "-line.bin"))
            for m, alg in INTERP.items():
                v = interpolate(ds, srs, x, y, alg)
                v.astype("<f8").tofile(os.path.join(OUT, "%s-%s.bin" % (name, m)))
            out["cases"].append({"name": name, "source": "cog/" + f, "crs": crs, "n": N,
                                 "bands": nb, "found": int(found.sum())})
        ds = None
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
