"""An mdim-refs/0.2 byte-reference table (spec: wxylib/design,
specs/mdim-refs-0.2.md) over three tiny OISST-shaped
netCDF4 files, and what GDAL reads for a few windows of the mosaic, so a
table reader can be checked against GDAL (mdimrefs/oracle/*.bin).

Steps, as the full OISST archive table is built:

1. write three daily netCDF4 files shaped like OISST v2.1 (time 1, zlev 1,
   lat 36, lon 72, int16, shuffle then deflate level 1, one chunk each);
2. `gdal mdim mosaic` them to a multidimensional VRT;
3. scan each source chunk's byte range with stock GDAL
   (MDArray.GetRawBlockInfo), keyed within its source (physical layer),
   copy the VRT's SourceSlab and DestSlab into `placements` (plan layer),
   and define each `refs_<array>` view from the two (view layer), in an
   mdim-refs/0.2 SQLite table;
4. export each `refs_<array>` view to Parquet plus `arrays.json`, the
   form the Rust reader (wxy-virtual) takes;
5. read windows of the mosaic with GDAL into oracle/*.bin (little-endian,
   C order over the array's dims).

The `ice` chunk of the second day is never written, so it has no bytes in
its file: its table row has present = 0 and every reader must give the
fill value there.

Source paths are recorded verbatim as if the files sat in a local mirror
at /mirror/oisst/; the `fixture` remap target rewrites that prefix to the
`nc/` folder next to the table.

Needs GDAL >= 3.12 with its Python bindings and netCDF driver, netCDF4,
numpy and pyarrow, for example:

    micromamba create -n wxy -c conda-forge "gdal>=3.12" libgdal-netcdf netcdf4 numpy pyarrow
    python scripts/make_mdimrefs_fixtures.py

Content is ASCII only.
"""

import datetime
import json
import os
import shutil
import sqlite3
import xml.etree.ElementTree as ET

import netCDF4
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from osgeo import gdal

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "mdimrefs")
NC = os.path.join(OUT, "nc")
ORACLE = os.path.join(OUT, "oracle")
MIRROR = "/mirror/oisst/"

FORMAT = "mdim-refs/0.2"
DAYS = ["19810901", "19810902", "19810903"]
NLAT, NLON = 36, 72
FILL = -999
CODEC = {
    "compressor": {"id": "zlib", "level": 1},
    "filters": [{"id": "shuffle", "elementsize": 2}],
}
DATA = ["sst", "ice"]

# Windows as (array, start, count) over (time, zlev, lat, lon).
WINDOWS = [
    ("sst", [0, 0, 0, 0], [1, 1, NLAT, NLON]),
    ("sst", [0, 0, 5, 10], [3, 1, 20, 33]),
    ("sst", [2, 0, 30, 60], [1, 1, 6, 12]),
    ("ice", [0, 0, 0, 0], [3, 1, NLAT, NLON]),
    ("ice", [1, 0, 28, 0], [2, 1, 8, 20]),
]


def source_name(day):
    return f"{day[:6]}/oisst-avhrr-v02r01.{day}.nc"


def write_day(path, d):
    lat = -87.5 + 5.0 * np.arange(NLAT)
    lon = 2.5 + 5.0 * np.arange(NLON)
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing="ij")
    sst = 1500 * np.cos(la) + 300 * np.sin(2 * lo + 0.3 * d) - 200
    sst = np.where(np.abs(np.degrees(la)) > 80, FILL, np.round(sst))
    ice = np.where(np.abs(np.degrees(la)) > 60, np.round(50 + 40 * np.sin(lo + d)), FILL)

    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.Conventions = "CF-1.6"
    ds.title = "wxylib fixture shaped like NOAA OISST v2.1 (synthetic values)"
    for name, size in [("time", None), ("zlev", 1), ("lat", NLAT), ("lon", NLON)]:
        ds.createDimension(name, size)
    t = ds.createVariable("time", "f4", ("time",))
    t.units = "days since 1978-01-01 12:00:00"
    t.calendar = "standard"
    t[:] = [1339 + d]
    z = ds.createVariable("zlev", "f4", ("zlev",))
    z.units = "meters"
    z[:] = [0.0]
    v = ds.createVariable("lat", "f4", ("lat",))
    v.units = "degrees_north"
    v[:] = lat
    v = ds.createVariable("lon", "f4", ("lon",))
    v.units = "degrees_east"
    v[:] = lon
    for name, values, attrs in [
        ("sst", sst, dict(long_name="Daily sea surface temperature", units="Celsius",
                          scale_factor=np.float32(0.01), add_offset=np.float32(0.0))),
        ("ice", ice, dict(long_name="Sea ice concentration", units="%",
                          scale_factor=np.float32(0.01), add_offset=np.float32(0.0))),
    ]:
        var = ds.createVariable(name, "i2", ("time", "zlev", "lat", "lon"), zlib=True,
                                complevel=1, shuffle=True, chunksizes=(1, 1, NLAT, NLON),
                                fill_value=FILL)
        var.set_auto_maskandscale(False)
        for k, a in attrs.items():
            var.setncattr(k, a)
        if name == "ice" and d == 1:
            continue  # never written: the chunk has no bytes
        var[:] = values.astype("i2").reshape(1, 1, NLAT, NLON)
    ds.close()


def attrs_of(arr):
    out = {}
    for a in arr.GetAttributes():
        v = a.Read()
        if isinstance(v, (list, tuple)) and len(v) == 1:
            v = v[0]
        if isinstance(v, float):
            v = float(np.float32(v)) if a.GetDataType().GetNumericDataType() == gdal.GDT_Float32 else v
        out[a.GetName()] = v
    # GDAL lifts these CF attributes out of the attribute list.
    if arr.GetUnit():
        out["units"] = arr.GetUnit()
    if arr.GetScale() is not None:
        out["scale_factor"] = arr.GetScale()
        out["add_offset"] = arr.GetOffset() or 0.0
    return out


def coord_rows(vrt_path, vrt_root):
    """Coordinate arrays from the mosaic VRT: regular spacing as affine,
    anything else inline. Attributes come from the first source."""
    first = vrt_root.find("Group/Array/Source/SourceFilename").text
    src = gdal.OpenEx(os.path.join(os.path.dirname(vrt_path), first), gdal.OF_MULTIDIM_RASTER)
    rows = []
    group = vrt_root.find("Group")
    for arr in group.findall("Array"):
        name = arr.get("name")
        if name in DATA:
            continue
        dim = arr.find("DimensionRef").get("ref")
        size = int(group.find(f"Dimension[@name='{dim}']").get("size"))
        reg = arr.find("RegularlySpacedValues")
        inline = arr.find("InlineValues")
        row = dict(array=name, kind="coord", dims=[dim], shape=[size], chunks=[size],
                   dtype="<f8", fill=None, codec=None,
                   attrs=attrs_of(src.GetRootGroup().OpenMDArray(name)) or None,
                   affine=None, values=None)
        if reg is not None:
            row["storage"] = "affine"
            row["affine"] = {"start": float(reg.get("start")), "step": float(reg.get("increment"))}
        else:
            row["storage"] = "inline"
            row["values"] = [float(x) for x in inline.text.split()]
        rows.append(row)
    return rows


def scan(vrt_path, vrt_root):
    """The physical and plan layers of each data array: one row per chunk of
    each source, keyed on its position within the source, from the source's
    raw block info; and one placement per VRT <Source>."""
    out = {}
    group = vrt_root.find("Group")
    for arr in group.findall("Array"):
        name = arr.get("name")
        if name not in DATA:
            continue
        chunks = [int(x) for x in arr.find("BlockSize").text.split(",")]
        dims = [d.get("ref") for d in arr.findall("DimensionRef")]
        gshape = [int(group.find(f"Dimension[@name='{d}']").get("size")) for d in dims]
        rows, placements, infos = [], [], set()
        for src in arr.findall("Source"):
            fname = src.find("SourceFilename").text
            source = MIRROR + fname[len("nc/"):]
            slab = src.find("SourceSlab")
            src_offset = [int(x) for x in slab.get("offset").split(",")]
            count = [int(x) for x in slab.get("count").split(",")]
            step = [int(x) for x in slab.get("step").split(",")]
            dest = [int(x) for x in src.find("DestSlab").get("offset").split(",")]
            # Chunk-aligned plan (spec rule).
            if (any(x != 1 for x in step)
                    or any(o % c for o, c in zip(src_offset + dest, chunks + chunks))
                    or any((d + n) % c and d + n != g for d, n, c, g in zip(dest, count, chunks, gshape))):
                raise ValueError(f"{fname}: placement is not on the chunk grid {chunks}")
            placements.append(dict(source=source, source_array=src.find("SourceArray").text,
                                   src_offset=src_offset, count=count,
                                   dest=dict(offset=dest, count=count)))
            ds = gdal.OpenEx(os.path.join(os.path.dirname(vrt_path), fname), gdal.OF_MULTIDIM_RASTER)
            a = ds.GetRootGroup().OpenMDArray(name)
            if a.GetBlockSize() != chunks:
                raise ValueError(f"{fname}: source blocks {a.GetBlockSize()} differ from {chunks}")
            shape = [d.GetSize() for d in a.GetDimensions()]
            counts = [-(-s // c) for s, c in zip(shape, chunks)]
            for local in np.ndindex(*counts):
                try:
                    info = a.GetRawBlockInfo(list(local))
                except RuntimeError:
                    # GDAL 3.13 raises for an unallocated HDF5 chunk. Call it
                    # absent only if GDAL reads that block as all fill.
                    start = [i * c for i, c in zip(local, chunks)]
                    n = [min(c, s - o) for c, s, o in zip(chunks, shape, start)]
                    block = a.ReadAsArray(array_start_idx=start, count=n)
                    if not np.all(block == FILL):
                        raise
                    info = None
                if info is None or info.GetOffset() is None or info.GetSize() == 0:
                    rows.append(dict(source=source, key=list(local), present=0,
                                     offset=None, size=None, info=None))
                else:
                    text = json.dumps(list(info.GetInfo() or []))
                    infos.add(text)
                    rows.append(dict(source=source, key=list(local), present=1,
                                     offset=info.GetOffset(), size=info.GetSize(), info=text))
        if len(infos) > 1:
            raise ValueError(f"{name}: chunks disagree on codec info {infos}")
        first = gdal.OpenEx(os.path.join(os.path.dirname(vrt_path), arr.find("Source/SourceFilename").text),
                            gdal.OF_MULTIDIM_RASTER).GetRootGroup().OpenMDArray(name)
        out[name] = dict(
            meta=dict(array=name, kind="data", storage="referenced", dims=dims, shape=gshape,
                      chunks=chunks, dtype="<i2", fill=str(FILL), codec=CODEC,
                      attrs={k: v for k, v in attrs_of(first).items() if k != "_FillValue"},
                      affine=None, values=None),
            rows=rows,
            placements=placements,
        )
    return out


def view_sql(name, n, chunks):
    """refs_<name>: chunk table joined to its placement, keys made global,
    restricted to source chunks inside the slab."""
    jx = lambda col, i: f"json_extract(p.{col}, '$[{i}]')"
    jd = lambda i: f"json_extract(p.dest, '$.offset[{i}]')"
    dims = ", ".join(f"t.dim_{i} - {jx('src_offset', i)} / {c} + {jd(i)} / {c} AS dim_{i}"
                     for i, c in enumerate(chunks))
    inside = " AND ".join(f"t.dim_{i} * {c} >= {jx('src_offset', i)} AND "
                          f"t.dim_{i} * {c} < {jx('src_offset', i)} + {jx('count', i)}"
                          for i, c in enumerate(chunks))
    return (f'CREATE VIEW "refs_{name}" AS SELECT {dims}, t.present, '
            f'CASE WHEN t."offset" IS NULL THEN NULL ELSE s.source END AS path, '
            f't."offset", t.size, t.info, t.data '
            f'FROM "{name}" t JOIN placements p ON p.array = \'{name}\' AND p.source_id = t.source_id '
            f'JOIN sources s ON s.source_id = t.source_id WHERE {inside}')


def write_sqlite(path, coords, data, files):
    if os.path.exists(path):
        os.remove(path)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    con.executemany("INSERT INTO meta VALUES(?, ?)",
                    [("format_version", FORMAT), ("plan", "oisst.vrt")])
    con.execute("CREATE TABLE remap(target TEXT, prefix TEXT, replacement TEXT, "
                "PRIMARY KEY(target, prefix)) WITHOUT ROWID")
    con.execute("INSERT INTO remap VALUES('fixture', ?, 'nc/')", (MIRROR,))
    con.execute('CREATE TABLE arrays(array TEXT PRIMARY KEY, kind TEXT, storage TEXT, dims TEXT, '
                'shape TEXT, chunks TEXT, dtype TEXT, fill TEXT, codec TEXT, attrs TEXT, '
                'affine TEXT, "values" BLOB)')
    j = lambda v: None if v is None else json.dumps(v)
    for r in coords + [d["meta"] for d in data.values()]:
        blob = None if r["values"] is None else np.asarray(r["values"], "<f8").tobytes()
        con.execute("INSERT INTO arrays VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (r["array"], r["kind"], r["storage"], j(r["dims"]), j(r["shape"]),
                     j(r["chunks"]), r["dtype"], r["fill"], j(r["codec"]), j(r["attrs"]),
                     j(r["affine"]), blob))
    con.execute("CREATE TABLE sources(source_id INTEGER PRIMARY KEY, source TEXT UNIQUE, "
                "size BIGINT, mtime BIGINT, etag TEXT)")
    sid = {}
    for i, (source, local) in enumerate(sorted(files.items())):
        st = os.stat(local)
        con.execute("INSERT INTO sources VALUES(?, ?, ?, ?, NULL)",
                    (i + 1, source, st.st_size, int(st.st_mtime)))
        sid[source] = i + 1
    con.execute("CREATE TABLE placements(array TEXT, dest TEXT, source_id INTEGER, "
                "source_array TEXT, src_offset TEXT, count TEXT, "
                "PRIMARY KEY(array, dest)) WITHOUT ROWID")
    for name, d in data.items():
        for p in d["placements"]:
            con.execute("INSERT INTO placements VALUES(?,?,?,?,?,?)",
                        (name, j(p["dest"]), sid[p["source"]], p["source_array"],
                         j(p["src_offset"]), j(p["count"])))
        n = len(d["meta"]["dims"])
        dimcols = ", ".join(f"dim_{i} INTEGER" for i in range(n))
        key = ", ".join(["source_id"] + [f"dim_{i}" for i in range(n)])
        con.execute(f'CREATE TABLE "{name}"(source_id INTEGER, {dimcols}, present INTEGER, '
                    f'"offset" BIGINT, size BIGINT, info TEXT, data BLOB, '
                    f'PRIMARY KEY({key})) WITHOUT ROWID')
        for r in d["rows"]:
            con.execute(f'INSERT INTO "{name}" VALUES({",".join("?" * (n + 6))})',
                        (sid[r["source"]], *r["key"], r["present"], r["offset"], r["size"],
                         r["info"], None))
        con.execute(view_sql(name, n, d["meta"]["chunks"]))
    con.commit()
    con.execute("VACUUM")
    con.execute("PRAGMA journal_mode=DELETE")
    con.close()


def export(sqlite_path):
    """Parquet of each refs_ view plus arrays.json, read back from the
    SQLite table so the export is the table, not a second copy of the scan."""
    con = sqlite3.connect(sqlite_path)
    arrays = []
    for row in con.execute('SELECT array, kind, storage, dims, shape, chunks, dtype, fill, codec, '
                           'attrs, affine, "values" FROM arrays ORDER BY array'):
        r = dict(zip(["array", "kind", "storage", "dims", "shape", "chunks", "dtype", "fill",
                      "codec", "attrs", "affine", "values"], row))
        for k in ["dims", "shape", "chunks", "codec", "attrs", "affine"]:
            r[k] = None if r[k] is None else json.loads(r[k])
        if r["values"] is not None:
            r["values"] = np.frombuffer(r["values"], "<f8").tolist()
        arrays.append(r)
        if r["storage"] != "referenced":
            continue
        n = len(r["dims"])
        cur = con.execute(f'SELECT * FROM "refs_{r["array"]}" ORDER BY '
                          + ", ".join(f"dim_{i}" for i in range(n)))
        cols = [c[0] for c in cur.description]
        vals = list(zip(*cur.fetchall()))
        types = {**{f"dim_{i}": pa.int64() for i in range(n)}, "present": pa.int8(),
                 "path": pa.string(), "offset": pa.int64(), "size": pa.int64(), "info": pa.string(),
                 "data": pa.binary()}
        table = pa.table({c: pa.array(v, types[c]) for c, v in zip(cols, vals)})
        pq.write_table(table, os.path.join(OUT, f"refs_{r['array']}.parquet"), compression="zstd")
    remap = [dict(target=t, prefix=p, replacement=s)
             for t, p, s in con.execute("SELECT target, prefix, replacement FROM remap")]
    fmt = con.execute("SELECT value FROM meta WHERE key = 'format_version'").fetchone()[0]
    con.close()
    with open(os.path.join(OUT, "arrays.json"), "w") as f:
        json.dump(dict(format_version=fmt, arrays=arrays, remap=remap), f, indent=1)


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(ORACLE)
    files = []
    for d, day in enumerate(DAYS):
        p = os.path.join(NC, source_name(day))
        os.makedirs(os.path.dirname(p), exist_ok=True)
        write_day(p, d)
        # A fixed mtime (the day itself), so sources.mtime is reproducible.
        t = int(datetime.datetime(int(day[:4]), int(day[4:6]), int(day[6:]),
                                  tzinfo=datetime.timezone.utc).timestamp())
        os.utime(p, (t, t))
        files.append(p)

    vrt = os.path.join(OUT, "oisst.vrt")
    gdal.Run("mdim", "mosaic", input=files, output=vrt)
    root = ET.parse(vrt).getroot()

    coords = coord_rows(vrt, root)
    data = scan(vrt, root)
    sqlite_path = os.path.join(OUT, "oisst.sqlite")
    write_sqlite(sqlite_path, coords, data,
                 {MIRROR + source_name(day): f for day, f in zip(DAYS, files)})
    export(sqlite_path)

    ds = gdal.OpenEx(vrt, gdal.OF_MULTIDIM_RASTER)
    windows = []
    for i, (name, start, count) in enumerate(WINDOWS):
        a = ds.GetRootGroup().OpenMDArray(name)
        v = a.ReadAsArray(array_start_idx=start, count=count).astype("<i2")
        f = f"oracle/{name}-w{i}.bin"
        v.tofile(os.path.join(OUT, f))
        windows.append(dict(array=name, start=start, count=count, file=f))

    manifest = dict(
        gdal=gdal.__version__,
        netCDF4=netCDF4.__version__,
        format_version=FORMAT,
        sources=[source_name(d) for d in DAYS],
        remap_target="fixture",
        arrays={name: dict(shape=d["meta"]["shape"], chunks=d["meta"]["chunks"],
                           dtype="int16", fill=FILL,
                           present=sum(r["present"] for r in d["rows"]),
                           absent=sum(1 - r["present"] for r in d["rows"]))
                for name, d in data.items()},
        coords={c["array"]: c["affine"] if c["storage"] == "affine" else c["values"] for c in coords},
        windows=windows,
    )
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)


if __name__ == "__main__":
    main()
