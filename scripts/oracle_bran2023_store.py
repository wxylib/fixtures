"""GDAL's reading of a BRAN2023 `temp` reference store, for comparison with
wxylib/core (`cargo run -p wxy-virtual --example bran_refs -- <store>/flat
--direct`). Not used in CI.

Opens the store's multidimensional VRT (`mdim.vrt`, the plan) with
osgeo.gdal and prints the SHA-256 of the raw int16 values (little-endian)
in the two windows the example reads when it has no manifest: the first
day's chunk at row 1200, and row 750 of depth 10 on the last day.

    python scripts/oracle_bran2023_store.py \\
        /vsis3/aadc-transfer/michae_sum/e2e/bran2023-daily-temp/mdim.vrt

S3 settings come from the environment as GDAL reads them
(AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_S3_ENDPOINT,
AWS_VIRTUAL_HOSTING=NO); the chunks are read from the netCDF files the
VRT names, so this also needs access to thredds.nci.org.au.
"""

import hashlib
import sys

import numpy as np
from osgeo import gdal

VAR = "temp"


def main(vrt):
    gdal.UseExceptions()
    ds = gdal.OpenEx(vrt, gdal.OF_MULTIDIM_RASTER)
    root = ds.GetRootGroup()
    a = root.OpenMDArray(VAR) or root.ResolveMDArray(VAR, "/")
    shape = [d.GetSize() for d in a.GetDimensions()]
    print(f"GDAL {gdal.__version__}: {VAR} shape {shape}")
    windows = {
        "chunk-t0-z0-y4-x0": ([0, 0, 1200, 0], [1, 1, 300, 300]),
        "row-tlast-z10-y750": ([shape[0] - 1, 10, 750, 0], [1, 1, 1, 3600]),
    }
    for name, (start, count) in windows.items():
        buf = a.Read(array_start_idx=start, count=count,
                     buffer_datatype=gdal.ExtendedDataType.Create(gdal.GDT_Int16))
        v = np.frombuffer(buf, "=i2").astype("<i2")
        print(f"{name}: sha256 {hashlib.sha256(v.tobytes()).hexdigest()}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         "/vsis3/aadc-transfer/michae_sum/e2e/bran2023-daily-temp/mdim.vrt")
