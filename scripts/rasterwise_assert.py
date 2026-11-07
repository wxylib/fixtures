#!/usr/bin/env python3
"""Write and check Assertion documents for rasterwise files.

usage: python3 scripts/rasterwise_assert.py [--check] [ids...]

An Assertion supplies what a file lacks, or corrects what it gets wrong,
without copying a byte. Each document here is keyed by the file's sha256,
points at the file's published url, and carries the assertion in forms GDAL
opens directly: a `vrt://` connection string, or a VRT document written
beside it (rasterwise/assertions/<id>.vrt). Each comes from the verified
recipe in the rasterwise case (`how_to_read`), restated so that it needs no
read-time configuration where possible.

--check opens every emitted form with osgeo.gdal and tests the claim:

  coords   pixel centres of the asserted grid, transformed to lon/lat,
           against the file's own longitude and latitude arrays read
           through the same reader path (so row order cancels); reports
           max and RMS great-circle distance in metres
  equals   the asserted view against the stored array (read with the
           multidimensional API in storage order) after a stated numpy
           rearrangement; exact
  warp     the asserted view warped to a lon/lat window against a
           reference file warped to the same window; RMS and median

Results go into each document's `check` field, with the GDAL version.
Design: wxylib/design docs/m30-legacy-corpus.md section 5.
"""
import datetime
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADIR = os.path.join(ROOT, "rasterwise", "assertions")
BLOB = "https://github.com/rgdal-dev/rasterwise/releases/download/blobs-1/"
CASES = "https://github.com/rgdal-dev/rasterwise/blob/main/cases/"


def url(sha, ext):
    return BLOB + sha + "." + ext


def nc(sha, var):
    return 'NETCDF:"/vsicurl/%s":%s' % (url(sha, "nc"), var)


def raw_vrt(src, nx, ny, dtype, srs, gt, image_offset, pixel_offset, line_offset, order):
    return """<VRTDataset rasterXSize="%d" rasterYSize="%d">
  <SRS>%s</SRS>
  <GeoTransform>%s</GeoTransform>
  <VRTRasterBand dataType="%s" band="1" subClass="VRTRawRasterBand">
    <SourceFilename relativeToVRT="0">/vsicurl/%s</SourceFilename>
    <ImageOffset>%d</ImageOffset>
    <PixelOffset>%d</PixelOffset>
    <LineOffset>%d</LineOffset>
    <ByteOrder>%s</ByteOrder>
  </VRTRasterBand>
</VRTDataset>
""" % (nx, ny, srs, ", ".join(repr(v) for v in gt), dtype, src, image_offset,
       pixel_offset, line_offset, order)


def swapped_vrt(raw_name):
    """-180..180 view of a 0..360 grid: the two halves of a raw companion VRT
    swapped. GDAL 3.8.4's warper cannot invert a 0..360 EPSG:3857 extent
    ("Cannot find coordinate operations"), so this is the form to warp."""
    gt = [-20037508.342789244, MERC_PX, 0.0, 16030006.674231391, 0.0, -MERC_PX]
    src = """    <SimpleSource>
      <SourceFilename relativeToVRT="1">%s</SourceFilename>
      <SourceBand>1</SourceBand>
      <SrcRect xOff="%d" yOff="0" xSize="10800" ySize="17280"/>
      <DstRect xOff="%d" yOff="0" xSize="10800" ySize="17280"/>
    </SimpleSource>
"""
    return """<VRTDataset rasterXSize="21600" rasterYSize="17280">
  <SRS>EPSG:3857</SRS>
  <GeoTransform>%s</GeoTransform>
  <VRTRasterBand dataType="Int16" band="1">
%s%s  </VRTRasterBand>
</VRTDataset>
""" % (", ".join(repr(v) for v in gt), src % (raw_name, 10800, 0), src % (raw_name, 0, 10800))


def vrturl(dsn, srs, gt):
    return "vrt://%s?a_srs=%s&a_gt=%s" % (dsn, srs, ",".join(repr(v) for v in gt))


# ------------------------------------------------------------ the assertions

CS2 = "13cc4140fb59f5a42dd525303e0ca42ed50cb64018b9ac9808241724cf5afe70"
NSIDC = "549a1b01dda6aa14aaa58037d5458d7d257834185f1e875d7f4fbc41e7e40477"
NHSCE = "f43715f25daf39e50cc896be98ccee7f12a792ed7f4e359396321e01246d131d"
GPM = "df20ee33704c86265fdd85de3254f5061649fd95ccb1cb7d5851cf824ac534e7"
FI_IMG = "710a226e331a2a74e04c32c8c0252dd8ad784d452c2e97e7cefe84bffdc45dd9"
FI_LAT = "9ea0d2cfeeea9b0f88567f6db08bfc4cc17b7b17a2c5d36065973c809117b981"
FI_LON = "57c853e812682947fbe9e4dc521c2795f43d6cfd843d29480943843b29d79104"
TOPO19 = "0804ca23ec395440ebbc7ee98ac14ec8055c51e9442650c72542afb13aeee0bc"
TOPO19NC = "a8744b2763eaf3256a14985343098fa813e85a882e5d8e77ffe1267796a2a008"
TOPO27 = "477d7e804c65e3ac4d67eff7a1fca27be7d3c1e4553cab5ca597ebe726825f66"
TOPO27NC = "254ed798d28354174e299ae62439981d41eda444564d64d1336f2af38ebc556f"

# GDAL's netCDF driver assumes bottom-up rows when a file has no 1D y
# coordinate, so its default view of these top-down files is south-up. A
# south-up geotransform (origin at the south edge, positive y step) places
# that default view correctly, with no configuration option.
NSIDC_S_25 = [-3950000.0, 25000.0, 0.0, -3950000.0, 0.0, 25000.0]
NSIDC_S_12 = [-3950000.0, 12500.0, 0.0, -3950000.0, 0.0, 12500.0]
NHSCE_PROJ = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +ellps=WGS84 +units=m +no_defs"
NHSCE_GT = [-8388172.4, 16776344.8 / 88, 0.0, -8388172.4, 0.0, 16776344.8 / 88]
CEA = "+proj=cea +lat_ts=-65 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
FI_GT = [-472755.311, 2000.0, 0.0, -13426012.821, 0.0, -2000.0]
MERC_PX = 1855.3248465545596
TOPO_GT = [0.0, MERC_PX, 0.0, 16030006.674231391, 0.0, -MERC_PX]
RAW_CONFIG = {"GDAL_VRT_RAWRASTERBAND_ALLOWED_SOURCE": "ALL"}

GPM_MDIM = """<VRTDataset>
  <Group name="/">
    <Dimension name="lat" size="536" type="HORIZONTAL_Y" indexingVariable="lat"/>
    <Dimension name="lon" size="1440" type="HORIZONTAL_X" indexingVariable="lon"/>
    <Array name="lat">
      <DataType>Float32</DataType>
      <DimensionRef ref="lat"/>
      <Source>
        <SourceFilename relativeToVRT="0">/vsicurl/%s</SourceFilename>
        <SourceArray>/nlat</SourceArray>
        <SourceView>[::-1]</SourceView>
      </Source>
    </Array>
    <Array name="lon">
      <DataType>Float32</DataType>
      <DimensionRef ref="lon"/>
      <Source>
        <SourceFilename relativeToVRT="0">/vsicurl/%s</SourceFilename>
        <SourceArray>/nlon</SourceArray>
      </Source>
    </Array>
    <Array name="precipRateNearSurfMean">
      <DataType>Float32</DataType>
      <DimensionRef ref="lat"/>
      <DimensionRef ref="lon"/>
      <SRS dataAxisToSRSAxisMapping="2,1">EPSG:4326</SRS>
      <NoDataValue>-9999.900390625</NoDataValue>
      <Source>
        <SourceFilename relativeToVRT="0">/vsicurl/%s</SourceFilename>
        <SourceArray>/precipRateNearSurfMean</SourceArray>
        <SourceTranspose>0,1,3,2</SourceTranspose>
        <SourceView>[0,0,::-1,:]</SourceView>
      </Source>
    </Array>
  </Group>
</VRTDataset>
""" % (url(GPM, "nc"), url(GPM, "nc"), url(GPM, "nc"))

SPECS = [
    {"id": CS2[:12], "sha256": CS2, "ext": "nc",
     "case": "cryosat2-sea-ice-cryptic-curvilinear", "array": "sea_ice_concentration",
     "asserts": {"crs": "EPSG:3412", "geotransform": NSIDC_S_25,
                 "row_order": "south-up as GDAL's netCDF driver presents it"},
     "evidence": "fitted grid: the file's 2D lat/lon are the NSIDC 25 km south grid",
     "vrt_url": vrturl(nc(CS2, "sea_ice_concentration"), "EPSG:3412", NSIDC_S_25),
     "check": {"kind": "coords", "lon": nc(CS2, "lon"), "lat": nc(CS2, "lat"), "tol_m": 50}},
    {"id": NSIDC[:12], "sha256": NSIDC, "ext": "nc",
     "case": "nsidc-polar-grid-no-coords-flipped", "array": "concentration",
     "asserts": {"crs": "EPSG:3412", "geotransform": NSIDC_S_12, "nodata": -128,
                 "row_order": "south-up as GDAL's netCDF driver presents it"},
     "evidence": "documented: global attribute grid = NSIDC; extent of the NSIDC south grid at 12.5 km",
     "vrt_url": vrturl(nc(NSIDC, "concentration"), "EPSG:3412", NSIDC_S_12),
     "check": {"kind": "equals", "array": "/concentration", "index": [0], "op": "flipud"}},
    {"id": NHSCE[:12], "sha256": NHSCE, "ext": "nc",
     "case": "nhsce-polar-stereo-labelled-longlat", "array": "snow_cover_extent",
     "asserts": {"crs": NHSCE_PROJ, "geotransform": NHSCE_GT,
                 "row_order": "south-up as GDAL's netCDF driver presents it"},
     "evidence": "fitted grid: the file's 2D latitude/longitude are this polar stereographic grid",
     "vrt_url": vrturl(nc(NHSCE, "snow_cover_extent"), NHSCE_PROJ, NHSCE_GT),
     "check": {"kind": "coords", "lon": nc(NHSCE, "longitude"), "lat": nc(NHSCE, "latitude"),
               "tol_m": 10, "nodata": -999.0}},
    {"id": GPM[:12], "sha256": GPM, "ext": "nc",
     "case": "gpm-dpr-swapped-axes", "array": "precipRateNearSurfMean",
     "asserts": {"crs": "EPSG:4326", "dims": "transpose (nlon, nlat) to (lat, lon); reverse lat",
                 "slice": {"AD": 0, "chn": 0}},
     "evidence": "documented: 1D nlat and nlon coordinate variables",
     "mdim_vrt": GPM_MDIM,
     "check": {"kind": "equals", "array": "/precipRateNearSurfMean", "index": [0, 0],
               "op": "transpose_flipud", "mdim": True}},
    {"id": FI_IMG[:12], "sha256": FI_IMG, "ext": "img",
     "case": "aad-fast-ice-raw-cea-grid-as-lonlat-arrays", "array": None,
     "asserts": {"crs": CEA, "geotransform": FI_GT, "dtype": "Byte", "shape": [425, 4300],
                 "layout": "raw, rows stored south to north; flipped by a negative line offset"},
     "evidence": "fitted grid: the companion lats.img/lons.img are this grid",
     "vrt": raw_vrt(url(FI_IMG, "img"), 4300, 425, "Byte", CEA, FI_GT, 424 * 4300, 1, -4300, "LSB"),
     "config": RAW_CONFIG,
     "companions": {
         FI_LAT[:12]: raw_vrt(url(FI_LAT, "img"), 4300, 425, "Float32", CEA, FI_GT, 424 * 17200, 4, -17200, "LSB"),
         FI_LON[:12]: raw_vrt(url(FI_LON, "img"), 4300, 425, "Float32", CEA, FI_GT, 424 * 17200, 4, -17200, "LSB")},
     "check": {"kind": "coords", "lon": FI_LON[:12] + ".vrt", "lat": FI_LAT[:12] + ".vrt", "tol_m": 1}},
    {"id": TOPO19[:12], "sha256": TOPO19, "ext": "img",
     "case": "ucsd-topo-raw-mercator-headers-disagree", "array": None,
     "asserts": {"crs": "EPSG:3857", "geotransform": TOPO_GT, "dtype": "Int16",
                 "shape": [17280, 21600], "layout": "raw big-endian, rows north to south",
                 "view": "the companion .raw.vrt is the stored 0..360 grid; the main VRT swaps its halves to -180..180"},
     "evidence": "documented layout; grid checked against the resampled topo_19.1.nc",
     "vrt": swapped_vrt(TOPO19[:12] + ".raw.vrt"),
     "companions": {TOPO19[:12] + ".raw": raw_vrt(url(TOPO19, "img"), 21600, 17280, "Int16",
                                                "EPSG:3857", TOPO_GT, 0, 2, 43200, "MSB")},
     "config": RAW_CONFIG,
     "check": {"kind": "warp", "ref": "/vsicurl/" + url(TOPO19NC, "nc"),
               "te": [140, 35, 145, 40], "tr": 1.0 / 60}},
    {"id": TOPO27[:12], "sha256": TOPO27, "ext": "img",
     "case": "ucsd-topo-raw-mercator-headers-disagree", "array": None,
     "asserts": {"crs": "EPSG:3857", "geotransform": TOPO_GT, "dtype": "Int16",
                 "shape": [17280, 21600], "layout": "raw big-endian, rows north to south",
                 "view": "the companion .raw.vrt is the stored 0..360 grid; the main VRT swaps its halves to -180..180"},
     "evidence": "documented layout; grid checked against the resampled topo_27.1.nc",
     "vrt": swapped_vrt(TOPO27[:12] + ".raw.vrt"),
     "companions": {TOPO27[:12] + ".raw": raw_vrt(url(TOPO27, "img"), 21600, 17280, "Int16",
                                                "EPSG:3857", TOPO_GT, 0, 2, 43200, "MSB")},
     "config": RAW_CONFIG,
     "check": {"kind": "warp", "ref": "/vsicurl/" + url(TOPO27NC, "nc"),
               "te": [140, 35, 145, 40], "tr": 1.0 / 60}},
]


# ---------------------------------------------------------------- writing

def write(spec):
    doc = {"format": "rasterwise-assertion/0",
           "target": {"sha256": spec["sha256"], "url": url(spec["sha256"], spec["ext"]),
                      "array": spec["array"]},
           "case": spec["case"], "recipe": CASES + spec["case"] + "/case.yaml",
           "asserts": spec["asserts"], "evidence": spec["evidence"], "emit": {}}
    if "vrt_url" in spec:
        doc["emit"]["vrt_url"] = spec["vrt_url"]
    for key, name in (("vrt", ".vrt"), ("mdim_vrt", ".mdim.vrt")):
        if key in spec:
            fn = spec["id"] + name
            with open(os.path.join(ADIR, fn), "w") as f:
                f.write(spec[key])
            doc["emit"][key] = fn
    for cid, text in spec.get("companions", {}).items():
        with open(os.path.join(ADIR, cid + ".vrt"), "w") as f:
            f.write(text)
        doc["emit"].setdefault("companions", []).append(cid + ".vrt")
    if spec.get("config"):
        doc["emit"]["config"] = spec["config"]
        doc["emit"]["config_note"] = ("GDAL 3.10 and later refuse a raw source outside the "
                                      "VRT's directory unless this option is set")
    path = os.path.join(ADIR, spec["id"] + ".json")
    if os.path.exists(path):
        old = json.load(open(path))
        if "check" in old:
            doc["check"] = old["check"]
    return doc, path


def dump(doc, path):
    with open(path, "w") as f:
        json.dump(doc, f, indent=1)
        f.write("\n")


# ---------------------------------------------------------------- checking

def haversine(lon1, lat1, lon2, lat2):
    import numpy as np
    r = 6371008.8
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def emitted(doc):
    e = doc["emit"]
    if "vrt_url" in e:
        return e["vrt_url"]
    return os.path.join(ADIR, e.get("vrt") or e.get("mdim_vrt"))


def check(doc, spec, gdal):
    import numpy as np
    from osgeo import osr
    c = spec["check"]
    out = {"kind": c["kind"], "gdal": gdal.__version__,
           "date": datetime.date.today().isoformat()}
    for k, v in spec.get("config", {}).items():
        gdal.SetConfigOption(k, v)
    dsn = emitted(doc)
    if c["kind"] == "coords":
        ds = gdal.Open(dsn)
        nx, ny = ds.RasterXSize, ds.RasterYSize
        gt = ds.GetGeoTransform()
        jj, ii = np.meshgrid(np.arange(nx) + 0.5, np.arange(ny) + 0.5)
        x = gt[0] + jj * gt[1] + ii * gt[2]
        y = gt[3] + jj * gt[4] + ii * gt[5]
        src = ds.GetSpatialRef()
        dst = osr.SpatialReference()
        dst.SetWellKnownGeogCS("WGS84")
        dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        tr = osr.CoordinateTransformation(src, dst)
        pts = np.array(tr.TransformPoints(np.column_stack([x.ravel(), y.ravel()])))
        lon_a, lat_a = pts[:, 0].reshape(ny, nx), pts[:, 1].reshape(ny, nx)

        def opened(d):
            return gdal.Open(d if "/" in d else os.path.join(ADIR, d))
        lon_f = opened(c["lon"]).ReadAsArray().astype("float64")
        lat_f = opened(c["lat"]).ReadAsArray().astype("float64")
        ok = np.isfinite(lon_f) & np.isfinite(lat_f) & (np.abs(lat_f) <= 90)
        if "nodata" in c:
            ok &= (lon_f != c["nodata"]) & (lat_f != c["nodata"])
        d = haversine(lon_a[ok], lat_a[ok], lon_f[ok], lat_f[ok])
        out.update({"cells": int(ok.sum()), "max_m": round(float(d.max()), 3),
                    "rms_m": round(float(np.sqrt(np.mean(d ** 2))), 3),
                    "tol_m": c["tol_m"], "pass": bool(d.max() <= c["tol_m"])})
    elif c["kind"] == "equals":
        src = "/vsicurl/" + doc["target"]["url"]
        stored = gdal.OpenEx(src, gdal.OF_MULTIDIM_RASTER).GetRootGroup() \
            .OpenMDArrayFromFullname(c["array"]).ReadAsArray()
        for i in c["index"]:
            stored = stored[i]
        if c["op"] == "flipud":
            want = stored[::-1, :]
        elif c["op"] == "transpose_flipud":
            want = stored.T[::-1, :]
        if c.get("mdim"):
            md = gdal.OpenEx(dsn, gdal.OF_MULTIDIM_RASTER)
            arr = md.GetRootGroup().OpenMDArray(doc["target"]["array"])
            got = arr.ReadAsArray()
            # the classic 2D view of the asserted array, with its geotransform
            cls = arr.AsClassicDataset(1, 0)
            out["geotransform"] = list(cls.GetGeoTransform())
            out["crs"] = arr.GetSpatialRef().GetAuthorityCode(None)
            got_c = cls.ReadAsArray()
            out["classic_equals"] = bool(np.array_equal(got_c, want, equal_nan=True))
        else:
            got = gdal.Open(dsn).ReadAsArray()
        out["shape"] = list(got.shape)
        out["pass"] = bool(np.array_equal(got, want, equal_nan=True))
        if "classic_equals" in out:
            out["pass"] = out["pass"] and out["classic_equals"]
    elif c["kind"] == "warp":
        te, r = c["te"], c["tr"]
        kw = dict(format="MEM", dstSRS="EPSG:4326", outputBounds=te, xRes=r, yRes=r)
        a = gdal.Warp("", dsn, resampleAlg="cubic", **kw).ReadAsArray().astype("float64")
        b = gdal.Warp("", c["ref"], resampleAlg="near", **kw).ReadAsArray().astype("float64")
        ok = np.isfinite(a) & np.isfinite(b)
        d = a[ok] - b[ok]
        out.update({"window": te, "cells": int(ok.sum()),
                    "rms_m": round(float(np.sqrt(np.mean(d ** 2))), 3),
                    "median_abs_m": round(float(np.median(np.abs(d))), 3),
                    "pass": bool(np.sqrt(np.mean(d ** 2)) < 10)})
    for k in spec.get("config", {}):
        gdal.SetConfigOption(k, None)
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_check = "--check" in sys.argv
    os.makedirs(ADIR, exist_ok=True)
    gdal = None
    if do_check:
        from osgeo import gdal
        gdal.UseExceptions()
    for spec in SPECS:
        if args and spec["id"] not in args:
            continue
        doc, path = write(spec)
        if do_check:
            try:
                doc["check"] = check(doc, spec, gdal)
            except Exception as e:
                doc["check"] = {"kind": spec["check"]["kind"], "gdal": gdal.__version__,
                                "error": str(e).strip(), "pass": False}
            print("%s %-44s %s" % (spec["id"], spec["case"], json.dumps(
                {k: v for k, v in doc["check"].items() if k not in ("gdal", "date")})),
                file=sys.stderr)
        dump(doc, path)


if __name__ == "__main__":
    main()
