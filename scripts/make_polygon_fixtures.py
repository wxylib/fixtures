"""Polygon-to-cell oracles (wxylib stage S4c): random polygons over each
COG fixture, given in EPSG:4326 and in the file's own CRS, and the cells
GDAL's rasterizer burns for each one on the file's grid, alone, with and
without ALL_TOUCHED (gdal.RasterizeLayer on a MEM raster, one polygon
in an OGR MEM layer, rings as given). Vertices in another CRS are
reprojected one by one, as gdal_rasterize does, before rasterizing.

Per case: `*-coords.bin` (little-endian float64 x, y pairs, every ring of
every polygon in order), `*-rings.json` (per polygon, its rings' point
counts; the first ring is the exterior) and `*-cells.json` (per polygon,
the flat cell indices row * width + col as [start, length] runs, for
"centre" and "all_touched").

Uses GDAL's own Python bindings (osgeo.gdal, ogr, osr) and numpy, from a
conda-forge env. Run from the repository root:

    python scripts/make_polygon_fixtures.py
"""

import json
import os
import shutil

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "polygons")
N = 12


def star(rng, cx, cy, r, k):
    """A closed star-shaped ring of k vertices around (cx, cy)."""
    # jittered angles, no gap wider than 81 degrees for k >= 8, so the
    # ring keeps a hole of radius 0.15 r inside it (valid polygons)
    a = (np.arange(k) + rng.uniform(0, 0.8, k)) * 2 * np.pi / k
    rad = r * rng.uniform(0.4, 1.0, k)
    x, y = cx + rad * np.cos(a), cy + rad * np.sin(a)
    ring = list(zip(x, y))
    return ring + [ring[0]]


def polygons(rng, box):
    w, h = box[2] - box[0], box[3] - box[1]
    out = []
    for i in range(N):
        cx = rng.uniform(box[0] - 0.05 * w, box[2] + 0.05 * w)
        cy = rng.uniform(box[1] - 0.05 * h, box[3] + 0.05 * h)
        r = rng.uniform(0.03, 0.25) * min(w, h)
        ext = star(rng, cx, cy, r, int(rng.integers(8, 30)))
        rings = [ext]
        if i % 3 == 0:
            # a hole, wound the other way
            rings.append(star(rng, cx, cy, 0.15 * r, 8)[::-1])
        if i % 4 == 1:
            # a rectangle with edges on cell boundaries is in the native
            # cases below; here a clockwise exterior
            rings[0] = rings[0][::-1]
        out.append(rings)
    return out


def runs(idx):
    """Flat indices as [start, length] runs."""
    out = []
    for i in idx:
        i = int(i)
        if out and out[-1][0] + out[-1][1] == i:
            out[-1][1] += 1
        else:
            out.append([i, 1])
    return out


def crs_string(srs):
    """EPSG:code where the file's CRS has one, else WKT."""
    if srs.GetAuthorityName(None) and srs.GetAuthorityCode(None):
        return "%s:%s" % (srs.GetAuthorityName(None), srs.GetAuthorityCode(None))
    return srs.ExportToWkt()


def apply_gt(gt, c, r):
    """(col, row) to (x, y), summed in the order affine's `t * (c, r)` does."""
    return (c * gt[1] + r * gt[2] + gt[0], c * gt[4] + r * gt[5] + gt[3])


def burn(rings, gt, width, height, all_touched):
    """The cells GDAL's rasterizer burns for one polygon (rings in the
    grid's CRS), as a flat uint8 array."""
    geom = ogr.Geometry(ogr.wkbPolygon)
    for r in rings:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in r:
            ring.AddPoint_2D(float(x), float(y))
        ring.CloseRings()
        geom.AddGeometryDirectly(ring)
    vds = ogr.GetDriverByName("MEM").CreateDataSource("")
    lyr = vds.CreateLayer("p", geom_type=ogr.wkbPolygon)
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometryDirectly(geom)
    lyr.CreateFeature(feat)
    rds = gdal.GetDriverByName("MEM").Create("", width, height, 1, gdal.GDT_Byte)
    rds.SetGeoTransform(gt)
    opts = ["MERGE_ALG=REPLACE"] + (["ALL_TOUCHED=TRUE"] if all_touched else [])
    gdal.RasterizeLayer(rds, [1], lyr, burn_values=[1], options=opts)
    return rds.GetRasterBand(1).ReadAsArray().ravel()


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    rng = np.random.default_rng(20261017)
    out = {"gdal": gdal.__version__, "cases": []}
    for f in sorted(os.listdir(os.path.join(ROOT, "cog"))):
        if not f.endswith(".tif"):
            continue
        path = os.path.join(ROOT, "cog", f)
        ds = gdal.Open(path)
        width, height = ds.RasterXSize, ds.RasterYSize
        dsrs = ds.GetSpatialRef()
        dsrs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        native = crs_string(dsrs)
        gt = ds.GetGeoTransform()
        ds = None
        b = (gt[0], gt[3] + gt[5] * height, gt[0] + gt[1] * width, gt[3])
        for crs in dict.fromkeys(["EPSG:4326", native]):
            if crs == native:
                box, ct = b, None
            else:
                srs = osr.SpatialReference()
                srs.SetFromUserInput(crs)
                srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                box = osr.CoordinateTransformation(dsrs, srs).TransformBounds(*b, 21)
                ct = osr.CoordinateTransformation(srs, dsrs)
            polys = polygons(rng, box)
            if crs == "EPSG:4326":
                polys = [[[(min(max(x, -180.0), 180.0), min(max(y, -89.5), 89.5))
                           for x, y in r] for r in p] for p in polys]
            if crs == native:
                # one polygon on cell edges exactly, and one through
                # cell centres: the knife edges of the rule
                e = [apply_gt(gt, c, r) for c, r in [(10, 5), (40, 5), (40, 30), (10, 30)]]
                polys.append([e + [e[0]]])
                m = [apply_gt(gt, c + 0.5, r + 0.5) for c, r in [(12, 7), (33, 7), (33, 21), (12, 21)]]
                polys.append([m + [m[0]]])
            coords, rings_out, cells = [], [], {"centre": [], "all_touched": []}
            for p in polys:
                rings_out.append([len(r) for r in p])
                for r in p:
                    coords.extend(r)
                if ct is None:
                    geo = p
                else:
                    geo = [[v[:2] for v in ct.TransformPoints(r)] for r in p]
                for mode, touched in [("centre", False), ("all_touched", True)]:
                    m = burn(geo, gt, width, height, touched)
                    cells[mode].append(runs(np.flatnonzero(m)))
            name = "%s-%s" % (f[:-4], "4326" if crs == "EPSG:4326" else "native")
            np.array(coords, dtype="<f8").tofile(os.path.join(OUT, name + "-coords.bin"))
            with open(os.path.join(OUT, name + "-rings.json"), "w") as fh:
                json.dump(rings_out, fh)
            with open(os.path.join(OUT, name + "-cells.json"), "w") as fh:
                json.dump(cells, fh)
            out["cases"].append({"name": name, "source": "cog/" + f, "crs": crs,
                                 "polygons": len(polys),
                                 "cells": sum(n for c in cells["centre"] for _, n in c)})
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
