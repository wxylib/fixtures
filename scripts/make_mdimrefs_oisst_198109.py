"""The mdim-refs/0.2 table for 30 days of real NOAA OISST v2.1 (September
1981), scanned over HTTP from the public NOAA bucket, for the P3 demo
(wxylib/core, `cargo run -p wxy-virtual --example oisst_refs`). Not used
in CI.

Same steps as make_mdimrefs_fixtures.py, on the real files: `gdal mdim
mosaic` over `/vsicurl/` URLs, block info from stock GDAL, the SQLite table,
its Parquet plus arrays.json export. Sources are recorded verbatim as the
`/vsicurl/https://...` URIs GDAL opened (spec, source rule); the `public`
remap target rewrites them to plain `https://`.

The oracle is not stored as bytes (2 MB per array): manifest.json records
the SHA-256 of what GDAL reads from the source netCDF for day 30 of each
array (little-endian, C order), which the demo compares with what it reads.

Output: mdimrefs/oisst-198109/ (oisst.sqlite, refs_*.parquet, arrays.json,
manifest.json). Needs network access to the bucket, and the same Python
environment as make_mdimrefs_fixtures.py. Content is ASCII only.
"""

import hashlib
import json
import os
import shutil
import sys
import xml.etree.ElementTree as ET

from osgeo import gdal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_mdimrefs_fixtures as mr  # noqa: E402

BUCKET = "https://noaa-cdr-sea-surface-temp-optimum-interpolation-pds.s3.amazonaws.com"
PREFIX = "/vsicurl/" + BUCKET + "/data/v2.1/avhrr/"
OUT = os.path.join(mr.ROOT, "mdimrefs", "oisst-198109")
DATA = ["sst", "anom", "err", "ice"]
DAYS = [f"198109{d:02d}" for d in range(1, 31)]


def main():
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    urls = [f"{PREFIX}{d[:6]}/oisst-avhrr-v02r01.{d}.nc" for d in DAYS]

    vrt = os.path.join(OUT, "oisst.vrt")
    gdal.Run("mdim", "mosaic", input=urls, output=vrt, array=DATA)
    root = ET.parse(vrt).getroot()

    coords = mr.coord_rows(vrt, root)
    data = mr.scan(vrt, root, data=DATA, source_of=lambda f: f)
    sqlite_path = os.path.join(OUT, "oisst.sqlite")
    mr.write_sqlite(sqlite_path, coords, data, {u: u for u in urls},
                    [("public", "/vsicurl/https://", "https://")], "oisst.vrt")
    mr.export(sqlite_path, OUT)
    os.remove(vrt)  # the plan is in the table (placements)

    last = gdal.OpenEx(urls[-1], gdal.OF_MULTIDIM_RASTER).GetRootGroup()
    oracle = {}
    for name in DATA:
        v = last.OpenMDArray(name).ReadAsArray()
        oracle[name] = dict(
            start=[len(DAYS) - 1, 0, 0, 0],
            count=list(v.shape),
            sha256=hashlib.sha256(v.astype(v.dtype.newbyteorder("<")).tobytes()).hexdigest(),
        )
    manifest = dict(
        gdal=gdal.__version__,
        format_version=mr.FORMAT,
        remap_target="public",
        days=[DAYS[0], DAYS[-1]],
        oracle_source=urls[-1],
        oracle=oracle,
    )
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)


if __name__ == "__main__":
    main()
