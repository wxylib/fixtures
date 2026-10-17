"""Curvilinear warp oracles (wxylib stage S5b): a source grid placed by 2D
longitude and latitude arrays (cell centres, as CF coordinates are), and
what GDAL's warper writes for it through its geolocation-array
transformer.

The source is a GeoTIFF with no geotransform, wrapped in a VRT whose
GEOLOCATION metadata points at a two-band longitude/latitude GeoTIFF with
GEOREFERENCING_CONVENTION=PIXEL_CENTER (what GDAL's netCDF driver sets
for CF coordinates). rasterio.warp.reproject then warps it as gdalwarp
does, with the approximate transformer at 0.125 pixels.

For wxylib, each case's arrays are also written raw (little-endian
float64 longitude and latitude, the source values in their own type).

Run from the repository root:

    python scripts/make_curvilinear_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "curvilinear")

METHODS = {"near": Resampling.nearest, "bilinear": Resampling.bilinear, "cubic": Resampling.cubic}


def swath(ny, nx):
    """A curved, rotated swath over south-eastern Australia and the
    Southern Ocean: rows bend, columns shear, spacing varies."""
    j, i = np.mgrid[0:ny, 0:nx].astype("float64")
    lon = 140.0 + 0.11 * i + 0.03 * j + 0.0004 * (j - ny / 2) ** 2
    lat = -40.0 - 0.09 * j + 0.02 * i + 1.5 * np.sin(i / nx * np.pi)
    return lon, lat


def polar(ny, nx):
    """A swath running over the south pole region, to polar stereographic:
    rows are circles of latitude, columns meridians, slightly skewed."""
    j, i = np.mgrid[0:ny, 0:nx].astype("float64")
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
]


def write_geoloc_vrt(path, src_tif, lonlat_tif, ny, nx):
    wkt = CRS.from_epsg(4326).to_wkt()
    md = {
        "X_DATASET": lonlat_tif, "X_BAND": "1",
        "Y_DATASET": lonlat_tif, "Y_BAND": "2",
        "PIXEL_OFFSET": "0", "LINE_OFFSET": "0",
        "PIXEL_STEP": "1", "LINE_STEP": "1",
        "SRS": wkt,
        "GEOREFERENCING_CONVENTION": "PIXEL_CENTER",
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
    out = {"gdal": rasterio.__gdal_version__, "rasterio": rasterio.__version__,
           "error_threshold": 0.125, "geolocation": "PIXEL_CENTER", "cases": []}
    grids = {}
    for c in CASES:
        ny, nx = c["shape"]
        key = (c["grid"].__name__, ny, nx)
        if key not in grids:
            lon, lat = c["grid"](ny, nx)
            v = field(lon, lat)
            base = "%s-%dx%d" % key
            src_tif = os.path.join(OUT, base + ".tif")
            lonlat_tif = os.path.join(OUT, base + "-lonlat.tif")
            prof = dict(driver="GTiff", width=nx, height=ny, count=1, dtype="float32")
            with rasterio.open(src_tif, "w", **prof) as d:
                d.write(v, 1)
            prof.update(count=2, dtype="float64")
            with rasterio.open(lonlat_tif, "w", **prof) as d:
                d.write(lon, 1)
                d.write(lat, 2)
            vrt = os.path.join(OUT, base + ".vrt")
            write_geoloc_vrt(vrt, src_tif, os.path.abspath(lonlat_tif), ny, nx)
            lon.astype("<f8").tofile(os.path.join(OUT, base + "-lon.bin"))
            lat.astype("<f8").tofile(os.path.join(OUT, base + "-lat.bin"))
            v.astype("<f4").tofile(os.path.join(OUT, base + "-values.bin"))
            grids[key] = base
        base = grids[key]
        w, h = c["size"]
        dst_transform = from_bounds(*c["bounds"], w, h)
        rec = {
            "name": c["name"],
            "source": base,
            "shape": [ny, nx],
            "source_crs": "EPSG:4326",
            "crs": c["crs"],
            "geotransform": list(dst_transform.to_gdal()),
            "width": w,
            "height": h,
            "dtype": "float32",
            "results": [],
        }
        with rasterio.open(os.path.join(OUT, base + ".vrt")) as src:
            for m, alg in METHODS.items():
                dst = np.zeros((h, w), dtype="float32")
                reproject(
                    source=rasterio.band(src, 1),
                    destination=dst,
                    dst_transform=dst_transform,
                    dst_crs=c["crs"],
                    resampling=alg,
                    num_threads=1,
                )
                name = "%s-%s.bin" % (c["name"], m)
                dst.astype("<f4").tofile(os.path.join(OUT, name))
                rec["results"].append({"resampling": m, "file": name})
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
