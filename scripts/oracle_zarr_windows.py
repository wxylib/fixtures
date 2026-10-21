"""GDAL's reading of windows of one array in a Zarr hierarchy, written as a
manifest for wxylib/core (`cargo run -p wxy-zarr --example zarr_windows --
<url> --manifest <file>`). Not used in CI.

Opens `ZARR:"/vsicurl/<url>"` (or a local path) with osgeo.gdal in
multidimensional mode and records the SHA-256 of the raw values
(little-endian, the array's own type) in each window:

    python scripts/oracle_zarr_windows.py <url> <array> <out.json> \\
        [NAME:START:COUNT ...]

START and COUNT are comma-separated indices, one per dimension. With no
window given, two are taken: the first chunk at index 0 of every
dimension, and one row of the last-but-one dimension's middle at the last
index of the first dimension.

Written for the Copernicus Marine ARCO stores (`timeChunked.zarr`, Zarr
v2 with consolidated metadata, buckets that refuse listing), for example

    https://s3.waw3-1.cloudferro.com/mdl-arco-time-025/arco/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083deg_P1D-m_202311/timeChunked.zarr
"""

import hashlib
import json
import sys

import numpy as np
from osgeo import gdal

NUMPY = {
    gdal.GDT_Byte: "u1", gdal.GDT_Int8: "i1", gdal.GDT_UInt16: "u2",
    gdal.GDT_Int16: "i2", gdal.GDT_UInt32: "u4", gdal.GDT_Int32: "i4",
    gdal.GDT_UInt64: "u8", gdal.GDT_Int64: "i8", gdal.GDT_Float32: "f4",
    gdal.GDT_Float64: "f8",
}


def dsn(location):
    if location.startswith(("http://", "https://")):
        location = "/vsicurl/" + location
    return f'ZARR:"{location}"'


def default_windows(shape, chunks):
    first = ([0] * len(shape), [min(c, s) for c, s in zip(chunks, shape)])
    start = [shape[0] - 1] + [0] * (len(shape) - 1)
    count = [1] * len(shape)
    if len(shape) >= 2:
        start[-2] = shape[-2] // 2
        count[-1] = shape[-1]
    return {"chunk-first": first, "row-tlast-ymid": (start, count)}


def parse_window(text):
    name, start, count = text.split(":")
    return name, ([int(v) for v in start.split(",")],
                  [int(v) for v in count.split(",")])


def main(location, array, out, specs):
    gdal.UseExceptions()
    ds = gdal.OpenEx(dsn(location), gdal.OF_MULTIDIM_RASTER)
    a = ds.GetRootGroup().OpenMDArray(array.strip("/"))
    shape = [d.GetSize() for d in a.GetDimensions()]
    chunks = a.GetBlockSize()
    dt = a.GetDataType().GetNumericDataType()
    windows = dict(map(parse_window, specs)) if specs else default_windows(shape, chunks)
    print(f"GDAL {gdal.__version__}: {array} shape {shape} chunks {chunks}")
    oracle = {}
    for name, (start, count) in windows.items():
        buf = a.Read(array_start_idx=start, count=count,
                     buffer_datatype=gdal.ExtendedDataType.Create(dt))
        v = np.frombuffer(buf, "=" + NUMPY[dt]).astype("<" + NUMPY[dt])
        h = hashlib.sha256(v.tobytes()).hexdigest()
        oracle[name] = {"start": start, "count": count, "sha256": h}
        print(f"{name}: sha256 {h}")
    with open(out, "w") as f:
        json.dump({"array": array.strip("/"), "location": location,
                   "gdal": gdal.__version__, "oracle": oracle}, f, indent=1)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
