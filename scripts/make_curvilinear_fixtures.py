"""Curvilinear warp oracles (wxylib stage S5b): a source grid placed by 2D
longitude and latitude arrays (cell centres, as CF coordinates are), and
what GDAL's warper writes for it through its geolocation-array
transformer.

The source is a GeoTIFF with no geotransform, wrapped in a VRT whose
GEOLOCATION metadata points at a two-band longitude/latitude GeoTIFF with
GEOREFERENCING_CONVENTION=PIXEL_CENTER (what GDAL's netCDF driver sets
for CF coordinates). osgeo.gdal's gdal.Warp then warps it with the
geolocation transformer, approximated at 0.125 pixels.

For wxylib, each case's arrays are also written raw (little-endian
float64 longitude and latitude, the source values in their own type).
Cases with "cell": "edge" place the grid by its cell corners instead:
arrays one row and column larger than the grid, read by GDAL with
GEOREFERENCING_CONVENTION=TOP_LEFT_CORNER.

Run from the repository root, with a Python that has GDAL's bindings and
numpy (e.g. a conda-forge env):

    python scripts/make_curvilinear_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import osgeo
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "curvilinear")

METHODS = ("near", "bilinear", "cubic")


def srs_wkt(crs):
    srs = osr.SpatialReference()
    srs.SetFromUserInput(crs)
    return srs.ExportToWkt()


def geotransform(bounds, w, h):
    xmin, ymin, xmax, ymax = bounds
    return [xmin, (xmax - xmin) / w, 0.0, ymax, 0.0, -(ymax - ymin) / h]


def write_tif(path, bands, dtype):
    """A plain GeoTIFF (no georeferencing) of 2D arrays, one per band."""
    ny, nx = bands[0].shape
    d = gdal.GetDriverByName("GTiff").Create(path, nx, ny, len(bands), dtype)
    for k, a in enumerate(bands):
        d.GetRasterBand(k + 1).WriteArray(a)
    d = None


def swath(ny, nx, shift=0.0):
    """A curved, rotated swath over south-eastern Australia and the
    Southern Ocean: rows bend, columns shear, spacing varies. With
    `shift` -0.5 and one more row and column, the cell corners."""
    j, i = np.mgrid[0:ny, 0:nx].astype("float64") + shift
    lon = 140.0 + 0.11 * i + 0.03 * j + 0.0004 * (j - ny / 2) ** 2
    lat = -40.0 - 0.09 * j + 0.02 * i + 1.5 * np.sin(i / nx * np.pi)
    return lon, lat


def polar(ny, nx, shift=0.0):
    """A swath running over the south pole region, to polar stereographic:
    rows are circles of latitude, columns meridians, slightly skewed."""
    j, i = np.mgrid[0:ny, 0:nx].astype("float64") + shift
    lon = -180.0 + 360.0 * (i + 0.5) / nx + 0.2 * j
    lon = (lon + 180.0) % 360.0 - 180.0
    lat = -60.0 - 25.0 * (j + 0.5) / ny
    return lon, lat


def field(lon, lat):
    v = 1000.0 * np.sin(np.radians(lon) * 3) * np.cos(np.radians(lat) * 2) + 50.0 * lat
    return v.astype("float32")


CASES = [
    dict(name="swath-to-4326", grid=swath, shape=(120, 160), crs="EPSG:4326",
         bounds=(139.0, -53.0, 160.0, -37.0), size=(210, 160)),
    # coarse: kernels stretched over about 3 source cells
    dict(name="swath-to-4326-coarse", grid=swath, shape=(120, 160), crs="EPSG:4326",
         bounds=(139.0, -53.0, 160.0, -37.0), size=(70, 53)),
    dict(name="swath-to-3031", grid=swath, shape=(120, 160), crs="EPSG:3031",
         bounds=(0.8e6, -6.0e6, 4.0e6, -2.8e6), size=(128, 128)),
    dict(name="polar-to-3031", grid=polar, shape=(50, 180), crs="EPSG:3031",
         bounds=(-3.5e6, -3.5e6, 3.5e6, 3.5e6), size=(140, 140)),
    # cell edges: (ny + 1, nx + 1) corner arrays, TOP_LEFT_CORNER
    dict(name="swath-edges-to-4326", grid=swath, shape=(120, 160), crs="EPSG:4326",
         bounds=(139.0, -53.0, 160.0, -37.0), size=(210, 160), cell="edge"),
    dict(name="swath-edges-to-3031", grid=swath, shape=(120, 160), crs="EPSG:3031",
         bounds=(0.8e6, -6.0e6, 4.0e6, -2.8e6), size=(128, 128), cell="edge"),
]


def write_geoloc_vrt(path, src_tif, lonlat_tif, ny, nx, convention):
    wkt = srs_wkt("EPSG:4326")
    md = {
        "X_DATASET": lonlat_tif, "X_BAND": "1",
        "Y_DATASET": lonlat_tif, "Y_BAND": "2",
        "PIXEL_OFFSET": "0", "LINE_OFFSET": "0",
        "PIXEL_STEP": "1", "LINE_STEP": "1",
        "SRS": wkt,
        "GEOREFERENCING_CONVENTION": convention,
    }
    items = "\n".join(
        '    <MDI key="%s">%s</MDI>' % (k, v.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;"))
        for k, v in md.items())
    with open(path, "w") as f:
        f.write("""<VRTDataset rasterXSize="%d" rasterYSize="%d">
  <Metadata domain="GEOLOCATION">
%s
  </Metadata>
  <VRTRasterBand dataType="Float32" band="1">
    <SimpleSource>
      <SourceFilename relativeToVRT="1">%s</SourceFilename>
      <SourceBand>1</SourceBand>
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>
""" % (nx, ny, items, os.path.basename(src_tif)))


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    out = {"gdal": gdal.VersionInfo("RELEASE_NAME"), "osgeo": osgeo.__version__,
           "error_threshold": 0.125,
           "geolocation": "PIXEL_CENTER, or TOP_LEFT_CORNER for cell edges", "cases": []}
    grids = {}
    for c in CASES:
        ny, nx = c["shape"]
        edge = c.get("cell") == "edge"
        key = (c["grid"].__name__, ny, nx, edge)
        if key not in grids:
            lon, lat = c["grid"](ny, nx)
            v = field(lon, lat)
            base = "%s-%dx%d" % key[:3] + ("-edges" if edge else "")
            if edge:
                lon, lat = c["grid"](ny + 1, nx + 1, -0.5)
            src_tif = os.path.join(OUT, base + ".tif")
            lonlat_tif = os.path.join(OUT, base + "-lonlat.tif")
            write_tif(src_tif, [v], gdal.GDT_Float32)
            write_tif(lonlat_tif, [lon, lat], gdal.GDT_Float64)
            vrt = os.path.join(OUT, base + ".vrt")
            write_geoloc_vrt(vrt, src_tif, os.path.abspath(lonlat_tif), ny, nx,
                             "TOP_LEFT_CORNER" if edge else "PIXEL_CENTER")
            lon.astype("<f8").tofile(os.path.join(OUT, base + "-lon.bin"))
            lat.astype("<f8").tofile(os.path.join(OUT, base + "-lat.bin"))
            v.astype("<f4").tofile(os.path.join(OUT, base + "-values.bin"))
            grids[key] = base
        base = grids[key]
        w, h = c["size"]
        gt = geotransform(c["bounds"], w, h)
        rec = {
            "name": c["name"],
            "source": base,
            "shape": [ny, nx],
            "cell": "edge" if edge else "centre",
            "source_crs": "EPSG:4326",
            "crs": c["crs"],
            "geotransform": gt,
            "width": w,
            "height": h,
            "dtype": "float32",
            "results": [],
        }
        src = gdal.Open(os.path.join(OUT, base + ".vrt"))
        for m in METHODS:
            # a zero-filled MEM target of the grid; no geotransform on the
            # source, so gdal.Warp uses its GEOLOCATION arrays
            dst = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Float32)
            dst.SetGeoTransform(gt)
            dst.SetProjection(srs_wkt(c["crs"]))
            gdal.Warp(dst, src, options=gdal.WarpOptions(
                resampleAlg=m, errorThreshold=0.125, multithread=False,
                transformerOptions=["METHOD=GEOLOC_ARRAY"],
                warpOptions=["INIT_DEST=0"]))
            name = "%s-%s.bin" % (c["name"], m)
            dst.ReadAsArray().astype("<f4").tofile(os.path.join(OUT, name))
            rec["results"].append({"resampling": m, "file": name})
        src = None
        out["cases"].append(rec)
    # the VRTs point at the lon/lat files by absolute path; they are for
    # regenerating here, not for use elsewhere
    for f in os.listdir(OUT):
        if f.endswith(".vrt"):
            os.remove(os.path.join(OUT, f))
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
