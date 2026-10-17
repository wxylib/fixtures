"""An mdim-refs/0.2 table for 32 days of BRAN2023 `temp` (CSIRO Bluelink
reanalysis, 2010-01-01 to 2010-02-01, all 51 depths), for the demo in
wxylib/core (`cargo run -p wxy-virtual --example bran_refs`). Not used in CI.

Unlike the OISST tables this one is not scanned with GDAL: it is converted
from an existing kerchunk Parquet index (fsspec layout, one positional row
per chunk) in mdsumner/virtualized, built in situ on NCI with R, rhdf5 and
GDAL multidim, with paths already rewritten to NCI's public THREDDS server.
It is the many-chunks-per-file case: 300 x 300 chunks, 3060 per file per
day, 97,920 in all, in two files.

The index is read from a pinned commit. Coordinates come from the inline
`raw` bytes of the index; an axis is stored as affine when its steps agree
to 1e-6 of a step (the index holds float32 values widened to float64).
Attributes are copied as the index has them (its `scale_factor` is rounded
to 0.0078; the netCDF says 0.007782219909131527).

With `--oracle` (GDAL 3.12 or later with the netCDF driver, and network
access to thredds.nci.org.au), manifest.json also records the SHA-256 of
what GDAL reads from the source netCDF for each oracle window (raw int16,
little-endian, C order), which the demo compares with what it reads.

Output: mdimrefs/bran2023-temp-201001/ (refs_temp.parquet, arrays.json,
manifest.json). Needs numpy and pyarrow. Content is ASCII only.
"""

import hashlib
import io
import json
import os
import sys
import urllib.request

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "mdimrefs", "bran2023-temp-201001")
COMMIT = "d0b9c23317d2de07c3851eb9dfc968be1761d5b8"
INDEX = f"https://raw.githubusercontent.com/mdsumner/virtualized/{COMMIT}/zarr/ocean_temp.zarr"
VAR = "temp"
NTIME = 32
THREDDS = "https://thredds.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2023/daily/"

# (start, count) in the subset, and the source file and its time index.
ORACLE = {
    "chunk-t0-z0-y4-x0": dict(start=[0, 0, 1200, 0], count=[1, 1, 300, 300],
                              file="ocean_temp_2010_01.nc", t=0),
    "row-t31-z10-y750": dict(start=[31, 10, 750, 0], count=[1, 1, 1, 3600],
                             file="ocean_temp_2010_02.nc", t=0),
}


def get(path):
    with urllib.request.urlopen(f"{INDEX}/{path}") as r:
        return r.read()


def main():
    meta = json.loads(get(".zmetadata"))
    rec, md = meta["record_size"], meta["metadata"]
    J = lambda v: v if isinstance(v, dict) else json.loads(v)

    def rows(name, lo, hi):
        parts = []
        for f in range(lo // rec, (hi - 1) // rec + 1):
            t = pq.read_table(io.BytesIO(get(f"{name}/refs.{f}.parq")))
            a, b = max(lo - f * rec, 0), min(hi - f * rec, rec)
            parts.append(t.slice(a, b - a))
        return pa.concat_tables(parts)

    def cf_attrs(attrs, drop):
        # R's index stores scalar attributes as length-one lists
        return {k: (x[0] if isinstance(x, list) and len(x) == 1 else x)
                for k, x in attrs.items() if k not in drop}

    def coord(name, n):
        za, attrs = J(md[f"{name}/.zarray"]), J(md.get(f"{name}/.zattrs", "{}"))
        v = np.frombuffer(rows(name, 0, 1).column("raw")[0].as_py(), za["dtype"])[:n].astype("<f8")
        step = (v[-1] - v[0]) / (len(v) - 1) if len(v) > 1 else 0.0
        regular = len(v) > 2 and np.allclose(np.diff(v), step, rtol=0, atol=abs(step) * 1e-6)
        row = dict(array=name, kind="coord", dims=[name], shape=[len(v)], chunks=[len(v)],
                   dtype="<f8", fill=None, codec=None,
                   attrs=cf_attrs(attrs, ("_ARRAY_DIMENSIONS", "_FillValue")))
        if regular:
            row.update(storage="affine", affine=dict(start=float(v[0]), step=float(step)), values=None)
        else:
            row.update(storage="inline", affine=None, values=[float(x) for x in v])
        return row

    za, attrs = J(md[f"{VAR}/.zarray"]), J(md[f"{VAR}/.zattrs"])
    dims = attrs["_ARRAY_DIMENSIONS"]
    shape = [NTIME] + za["shape"][1:]
    grid = [-(-s // c) for s, c in zip(shape, za["chunks"])]
    t = rows(VAR, 0, int(np.prod(grid)))
    keys = np.array(np.unravel_index(np.arange(t.num_rows), grid)).T
    path = t.column("path").to_pylist()
    cols = {f"dim_{i}": pa.array(keys[:, i].astype("i8")) for i in range(len(dims))}
    cols.update(
        present=pa.array([1 if p else 0 for p in path], pa.int8()),
        path=pa.array(path, pa.string()),
        offset=t.column("offset").cast(pa.int64()),
        size=t.column("size").cast(pa.int64()),
        info=pa.nulls(t.num_rows, pa.string()),
        data=t.column("raw"),
    )
    os.makedirs(OUT, exist_ok=True)
    pq.write_table(pa.table(cols), os.path.join(OUT, f"refs_{VAR}.parquet"))

    fill = za.get("fill_value")
    data = dict(array=VAR, kind="data", storage="referenced", dims=dims, shape=shape,
                chunks=za["chunks"], dtype=za["dtype"], fill=None if fill is None else str(fill),
                codec=dict(compressor=za.get("compressor"), filters=za.get("filters")),
                attrs=cf_attrs(attrs, ("_ARRAY_DIMENSIONS",)), affine=None, values=None)
    arrays = [data] + [coord(d, n) for d, n in zip(dims, shape)]
    with open(os.path.join(OUT, "arrays.json"), "w") as f:
        json.dump(dict(format_version="mdim-refs/0.2", arrays=arrays, remap=[]), f, indent=1)

    manifest = dict(format_version="mdim-refs/0.2", source_index=INDEX,
                    files=sorted(set(p for p in path if p)), chunks=t.num_rows,
                    oracle={k: dict(start=o["start"], count=o["count"]) for k, o in ORACLE.items()})
    if "--oracle" in sys.argv:
        from osgeo import gdal
        gdal.UseExceptions()
        manifest["gdal"] = gdal.__version__
        for k, o in ORACLE.items():
            ds = gdal.OpenEx("/vsicurl/" + THREDDS + o["file"], gdal.OF_MULTIDIM_RASTER)
            a = ds.GetRootGroup().OpenMDArray(VAR)
            start = [o["t"]] + o["start"][1:]
            buf = a.Read(array_start_idx=start, count=o["count"],
                         buffer_datatype=gdal.ExtendedDataType.Create(gdal.GDT_Int16))
            manifest["oracle"][k]["sha256"] = hashlib.sha256(
                np.frombuffer(buf, "=i2").astype("<i2").tobytes()).hexdigest()
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"{OUT}: {t.num_rows} chunks, shape {shape}")


if __name__ == "__main__":
    main()
