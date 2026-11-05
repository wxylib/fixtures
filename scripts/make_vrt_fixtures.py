# Writes the VRT conformance corpus (design m25, first steps 1): small
# sources, VRT documents and vrt:// strings, and what GDAL reads for each.
# Run from the repository root:
#
#   python scripts/make_vrt_fixtures.py
#
# GDAL is reached through osgeo.gdal when it is installed, otherwise through
# rasterio's bundled GDAL (same RasterIO path, nearest resampling); the
# manifest records which, and the GDAL version.
#
# Output (vrt/):
#   src/*.tif    uncompressed GeoTIFF sources, no overviews
#   src/*.bin    the same pixels, little-endian, C order, band-major
#   *.vrt        classic raster VRT documents (relative source paths)
#   oracle/*.bin GDAL's read of each case and request, little-endian, C order,
#                band-major, in the VRT band data type
#   manifest.json
# Content is ASCII only.

import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "vrt")
SRC = os.path.join(OUT, "src")
ORACLE = os.path.join(OUT, "oracle")
for d in (SRC, ORACLE):
    os.makedirs(d, exist_ok=True)

try:
    from osgeo import gdal, gdal_array
    gdal.UseExceptions()
    BACKEND = "osgeo.gdal"
    GDAL_VERSION = gdal.__version__
except ImportError:
    import rasterio
    from rasterio.windows import Window
    from rasterio.enums import Resampling
    from rasterio.transform import Affine
    BACKEND = "rasterio " + rasterio.__version__
    GDAL_VERSION = rasterio.__gdal_version__


def write_tif(path, data, gt, epsg=4326):
    bands, h, w = data.shape
    if BACKEND == "osgeo.gdal":
        etype = gdal_array.NumericTypeCodeToGDALTypeCode(data.dtype)
        ds = gdal.GetDriverByName("GTiff").Create(path, w, h, bands, etype)
        ds.SetGeoTransform(gt)
        ds.SetProjection("EPSG:%d" % epsg)
        for b in range(bands):
            ds.GetRasterBand(b + 1).WriteArray(data[b])
        ds = None
    else:
        tr = Affine.from_gdal(*gt)
        with rasterio.open(path, "w", driver="GTiff", width=w, height=h,
                           count=bands, dtype=data.dtype, crs="EPSG:%d" % epsg,
                           transform=tr) as ds:
            ds.write(data)


def read(dsn, window, buf):
    """GDAL's read of window (col, row, width, height) into buf (w, h),
    nearest resampling, all bands, in the dataset's data type."""
    c, r, w, h = window
    bw, bh = buf
    if BACKEND == "osgeo.gdal":
        ds = gdal.Open(dsn)
        a = ds.ReadAsArray(c, r, w, h, buf_xsize=bw, buf_ysize=bh,
                           resample_alg=gdal.GRIORA_NearestNeighbour)
        info = dict(size=[ds.RasterXSize, ds.RasterYSize], bands=ds.RasterCount,
                    dtype=gdal.GetDataTypeName(ds.GetRasterBand(1).DataType),
                    geotransform=list(ds.GetGeoTransform()),
                    nodata=ds.GetRasterBand(1).GetNoDataValue())
        if a.ndim == 2:
            a = a[None]
        return a, info
    with rasterio.open(dsn) as ds:
        a = ds.read(window=Window(c, r, w, h), out_shape=(ds.count, bh, bw),
                    resampling=Resampling.nearest)
        dt = {"uint8": "Byte", "int16": "Int16", "float32": "Float32",
              "float64": "Float64", "uint16": "UInt16", "int32": "Int32"}[ds.dtypes[0]]
        info = dict(size=[ds.width, ds.height], bands=ds.count, dtype=dt,
                    geotransform=list(ds.transform.to_gdal()),
                    nodata=ds.nodata)
        return a, info


def field(h, w, bands, dtype, seed):
    y, x = np.mgrid[0:h, 0:w].astype(np.float64)
    out = []
    for b in range(bands):
        v = 60 * np.sin(x / (5 + 3 * b) + seed) * np.cos(y / (4 + 2 * b)) + 3 * x + 2 * y + 40 * b + 10 * seed
        out.append(v)
    a = np.stack(out)
    if np.issubdtype(np.dtype(dtype), np.integer):
        info = np.iinfo(dtype)
        a = np.clip(np.round(a), info.min, info.max)
    return a.astype(dtype)


# ---------------------------------------------------------------------------
# Sources. Pixel size 1 degree; tiles of one 40 x 30 world.

SOURCES = {}

def source(name, data, gt):
    path = os.path.join(SRC, name + ".tif")
    write_tif(path, data, gt)
    data.astype(data.dtype.newbyteorder("<")).tofile(os.path.join(SRC, name + ".bin"))
    SOURCES[name] = dict(file="src/%s.tif" % name, raw="src/%s.bin" % name,
                         size=[data.shape[2], data.shape[1]], bands=data.shape[0],
                         dtype=str(data.dtype), geotransform=list(gt))

world = field(30, 40, 1, np.uint8, 1)
for i, (c, r) in enumerate([(0, 0), (20, 0), (0, 15), (20, 15)]):
    source("tile%d" % i, world[:, r:r + 15, c:c + 20].copy(),
           (100.0 + c, 1.0, 0.0, -10.0 - r, 0.0, -1.0))
source("small", field(10, 10, 1, np.int16, 2), (0.0, 1.0, 0.0, 10.0, 0.0, -1.0))
source("odd", field(30, 30, 1, np.int16, 3), (0.0, 1.0, 0.0, 30.0, 0.0, -1.0))
source("rgb", field(12, 16, 3, np.uint8, 4), (0.0, 1.0, 0.0, 12.0, 0.0, -1.0))
fl = field(12, 16, 1, np.float32, 5)
fl[0, 2:5, 3:9] = -9999.0
source("withnd", fl, (0.0, 1.0, 0.0, 12.0, 0.0, -1.0))


# ---------------------------------------------------------------------------
# VRT documents.

def src_xml(kind, name, band, srcrect, dstrect, extra=""):
    sw = SOURCES[name]
    return (
        "    <%s>\n"
        "      <SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n"
        "      <SourceBand>%s</SourceBand>\n"
        "      <SourceProperties RasterXSize=\"%d\" RasterYSize=\"%d\" DataType=\"%s\" BlockXSize=\"%d\" BlockYSize=\"1\"/>\n"
        "      <SrcRect xOff=\"%s\" yOff=\"%s\" xSize=\"%s\" ySize=\"%s\"/>\n"
        "      <DstRect xOff=\"%s\" yOff=\"%s\" xSize=\"%s\" ySize=\"%s\"/>\n"
        "%s"
        "    </%s>\n"
    ) % ((kind, sw["file"], band, sw["size"][0], sw["size"][1],
          {"uint8": "Byte", "int16": "Int16", "float32": "Float32"}[sw["dtype"]],
          sw["size"][0]) + tuple(srcrect) + tuple(dstrect) + (extra, kind))


def vrt_xml(size, gt, bands):
    """bands: list of (dtype, nodata or None, [source xml])"""
    s = "<VRTDataset rasterXSize=\"%d\" rasterYSize=\"%d\">\n" % tuple(size)
    s += "  <SRS dataAxisToSRSAxisMapping=\"2,1\">EPSG:4326</SRS>\n"
    s += "  <GeoTransform>%s</GeoTransform>\n" % ", ".join(repr(float(v)) for v in gt)
    for i, (dt, nd, srcs) in enumerate(bands):
        s += "  <VRTRasterBand dataType=\"%s\" band=\"%d\">\n" % (dt, i + 1)
        if nd is not None:
            s += "    <NoDataValue>%s</NoDataValue>\n" % nd
        s += "".join(srcs)
        s += "  </VRTRasterBand>\n"
    return s + "</VRTDataset>\n"


WORLD_GT = (100.0, 1.0, 0.0, -10.0, 0.0, -1.0)
CASES = []

def case(name, desc, text, size, requests=None):
    path = os.path.join(OUT, name)
    with open(path, "w") as f:
        f.write(text)
    CASES.append(dict(name=name, desc=desc, dsn=name, size=size, requests=requests))

tiles = [src_xml("SimpleSource", "tile%d" % i, 1, (0, 0, 20, 15), (c, r, 20, 15))
         for i, (c, r) in enumerate([(0, 0), (20, 0), (0, 15), (20, 15)])]
case("c01-mosaic.vrt", "2 x 2 mosaic of SimpleSources, as gdalbuildvrt writes",
     vrt_xml((40, 30), WORLD_GT, [("Byte", None, tiles)]), (40, 30))
case("c02-overlap.vrt", "two sources overlap; the later one wins",
     vrt_xml((40, 30), WORLD_GT, [("Byte", None, [
         src_xml("SimpleSource", "tile0", 1, (0, 0, 20, 15), (0, 0, 20, 15)),
         src_xml("SimpleSource", "tile3", 1, (0, 0, 20, 15), (10, 8, 20, 15))])]), (40, 30))
case("c03-upsample.vrt", "SrcRect 10 x 10 into DstRect 20 x 20 (implicit 2x nearest)",
     vrt_xml((20, 20), (0, 0.5, 0, 10, 0, -0.5), [("Int16", None, [
         src_xml("SimpleSource", "small", 1, (0, 0, 10, 10), (0, 0, 20, 20))])]), (20, 20))
case("c04-downsample.vrt", "SrcRect 30 x 30 into DstRect 20 x 20 (factor 1.5, as in COP30)",
     vrt_xml((20, 20), (0, 1.5, 0, 30, 0, -1.5), [("Int16", None, [
         src_xml("SimpleSource", "odd", 1, (0, 0, 30, 30), (0, 0, 20, 20))])]), (20, 20))
case("c05-fractional.vrt", "fractional DstRect offset and fractional SrcRect window",
     vrt_xml((24, 24), (0, 1, 0, 24, 0, -1), [("Int16", None, [
         src_xml("SimpleSource", "odd", 1, (2.5, 3.25, 20, 18), (1.5, 2.75, 20, 18))])]), (24, 24))
case("c06-complex.vrt", "ComplexSource: NODATA, ScaleOffset, ScaleRatio; band NoDataValue",
     vrt_xml((16, 12), (0, 1, 0, 12, 0, -1), [("Float32", -1, [
         src_xml("ComplexSource", "withnd", 1, (0, 0, 16, 12), (0, 0, 16, 12),
                 "      <NODATA>-9999</NODATA>\n      <ScaleOffset>273.15</ScaleOffset>\n      <ScaleRatio>0.5</ScaleRatio>\n")])]), (16, 12))
case("c07-lut.vrt", "ComplexSource with a LUT (linear interpolation between entries)",
     vrt_xml((16, 12), (0, 1, 0, 12, 0, -1), [("Float32", None, [
         src_xml("ComplexSource", "rgb", 2, (0, 0, 16, 12), (0, 0, 16, 12),
                 "      <LUT>0:0,50:100,100:150,255:155</LUT>\n")])]), (16, 12))
case("c08-bands.vrt", "band selection and reordering from a 3-band source",
     vrt_xml((16, 12), (0, 1, 0, 12, 0, -1), [
         ("Byte", None, [src_xml("SimpleSource", "rgb", 3, (0, 0, 16, 12), (0, 0, 16, 12))]),
         ("Byte", None, [src_xml("SimpleSource", "rgb", 1, (0, 0, 16, 12), (0, 0, 16, 12))])]), (16, 12))
nested = (
    "    <SimpleSource>\n"
    "      <SourceFilename relativeToVRT=\"1\">c01-mosaic.vrt</SourceFilename>\n"
    "      <SourceBand>1</SourceBand>\n"
    "      <SrcRect xOff=\"10\" yOff=\"5\" xSize=\"20\" ySize=\"20\"/>\n"
    "      <DstRect xOff=\"0\" yOff=\"0\" xSize=\"20\" ySize=\"20\"/>\n"
    "    </SimpleSource>\n")
case("c09-nested.vrt", "a VRT whose source is c01-mosaic.vrt, cropped across all four tiles",
     vrt_xml((20, 20), (110, 1, 0, -15, 0, -1), [("Byte", None, [nested])]), (20, 20))

# Requests: full read, a window at full resolution, and a downsampled read.
for c in CASES:
    w, h = c["size"]
    c["requests"] = [
        dict(window=[0, 0, w, h], buf=[w, h]),
        dict(window=[3, 2, w - 7, h - 5], buf=[w - 7, h - 5]),
        dict(window=[0, 0, w, h], buf=[(w + 1) // 2, (h + 2) // 3]),
    ]

# vrt:// strings over the sources (GDAL lowers them through gdal_translate).
STRINGS = [
    ("v01-bands", "vrt://src/rgb.tif?bands=3,1", (16, 12)),
    ("v02-srcwin", "vrt://src/odd.tif?srcwin=4,5,12,9", (12, 9)),
    ("v03-outsize", "vrt://src/odd.tif?outsize=50%,50%", (15, 15)),
    ("v04-srcwin-outsize", "vrt://src/odd.tif?srcwin=3,3,21,21&outsize=14,14", (14, 14)),
    ("v05-nodata", "vrt://src/withnd.tif?a_nodata=-9999", (16, 12)),
]
for name, s, size in STRINGS:
    w, h = size
    CASES.append(dict(name=name, desc="vrt:// string", dsn=s, size=list(size),
                      requests=[dict(window=[0, 0, w, h], buf=[w, h])]))

# ---------------------------------------------------------------------------
# Oracle.

cwd = os.getcwd()
os.chdir(OUT)  # vrt:// strings use paths relative to vrt/
try:
    for c in CASES:
        for k, rq in enumerate(c["requests"]):
            a, info = read(c["dsn"], rq["window"], rq["buf"])
            f = "oracle/%s-r%d.bin" % (c["name"].replace(".vrt", ""), k)
            a.astype(a.dtype.newbyteorder("<")).tofile(f)
            rq["file"] = f
            rq["shape"] = list(a.shape)
        c.update(info)
finally:
    os.chdir(cwd)

manifest = dict(gdal=GDAL_VERSION, backend=BACKEND,
                note="Requests use nearest resampling; window = [col, row, width, height], buf = [width, height].",
                sources=SOURCES, cases=CASES)
with open(os.path.join(OUT, "manifest.json"), "w") as f:
    json.dump(manifest, f, indent=1)
print("%d cases, %d requests, GDAL %s via %s" % (
    len(CASES), sum(len(c["requests"]) for c in CASES), GDAL_VERSION, BACKEND))
