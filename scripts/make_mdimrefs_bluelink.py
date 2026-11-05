"""A synthetic, Bluelink-shaped mdim-refs/0.2 table (m06 sec 8 point 3):
the many-chunks-per-file case of BRAN2023, small enough for CI.

Like BRAN2023 daily: one netCDF4 file per variable per month, with
time, st_ocean, yt_ocean, xt_ocean dims, int16 values packed with a
float64 scale_factor and add_offset, shuffle then deflate, and many
chunks per file. Here:

- two variables, `temp` and `salt`, each in two monthly files
  (`ocean_<var>_2010_01.nc` with 3 days, `ocean_<var>_2010_02.nc` with 2),
  so a time chunk's source file depends on where the month boundary falls;
- 4 depths (not evenly spaced, so `st_ocean` is inline), 40 x 60 cells,
  chunks of 1 x 1 x 10 x 15: 16 chunks per depth and day, 192 in the
  January file;
- chunks never written, so their rows have present = 0 and readers give
  the fill value (-32768): in `temp`, one chunk (day 2, depth 3, chunk
  row 1, column 0); in `salt`, a whole depth of a day (day 4, depth 2);
- `scale_factor` and `add_offset` as BRAN's netCDF has them (float64, not
  rounded), so a reader that rounds them is caught by the attribute check.

Unlike the OISST fixture this one needs no GDAL: chunk byte ranges come
from HDF5 itself (h5py `get_chunk_info_by_coord`), and the oracle windows
are what netCDF-C (netCDF4-python, mask and scale off) reads from the
source files, assembled over the global time axis. Both are independent
of the table reader under test.

Output: mdimrefs/bluelink-synthetic/ (nc/, refs_<var>.parquet,
arrays.json, manifest.json, oracle/*.bin). Source paths are recorded as
if mirrored at /mirror/bran/; the `fixture` remap target rewrites that
prefix to `nc/`, relative to the table.

Needs netCDF4, h5py, numpy and pyarrow (pip is enough). Content is ASCII
only.
"""

import datetime
import json
import os
import shutil

import h5py
import netCDF4
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "mdimrefs", "bluelink-synthetic")
NC = os.path.join(OUT, "nc")
MIRROR = "/mirror/bran/"
FORMAT = "mdim-refs/0.2"

MONTHS = [("2010_01", 3), ("2010_02", 2)]  # days per file (a cut-down month)
DEPTHS = [2.5, 7.5, 17.5, 40.0]
NY, NX = 40, 60
CHUNKS = [1, 1, 10, 15]
FILL = -32768
VARS = {
    # BRAN2023's packing, as its netCDF holds it
    "temp": dict(long_name="Potential temperature", units="degrees C",
                 scale_factor=0.007782219909131527, add_offset=210.0),
    "salt": dict(long_name="Practical Salinity", units="psu",
                 scale_factor=0.0012207403779029846, add_offset=35.0),
}
# never written, as (global day, depth, chunk row, chunk column); None is
# the whole depth
UNWRITTEN = {"temp": [(2, 3, 1, 0)], "salt": [(4, 2, None, None)]}

# Windows as (array, start, count) over (time, st_ocean, yt_ocean, xt_ocean).
WINDOWS = [
    ("temp", [0, 0, 0, 0], [5, 4, NY, NX]),       # everything, both files
    ("temp", [1, 1, 5, 7], [3, 3, 22, 40]),       # across the month boundary
    ("temp", [2, 3, 8, 0], [1, 1, 15, 20]),       # over the unwritten chunk
    ("salt", [3, 0, 0, 0], [2, 4, NY, NX]),       # February, unwritten depth
    ("salt", [0, 2, 33, 52], [5, 1, 7, 8]),       # a corner, all days
]


def values(var, day, depth):
    yy, xx = np.meshgrid(np.linspace(-1, 1, NY), np.linspace(-1, 1, NX), indexing="ij")
    if var == "temp":
        phys = 26 - 18 * yy ** 2 - 0.25 * DEPTHS[depth] + np.sin(3 * xx + 0.4 * day)
    else:
        phys = 35 + 0.6 * yy - 0.01 * DEPTHS[depth] + 0.2 * np.cos(2 * xx - 0.3 * day)
    a = VARS[var]
    packed = np.round((phys - a["add_offset"]) / a["scale_factor"])
    land = (xx > 0.7) & (yy > 0.5 - 0.1 * depth)  # land grows with depth
    return np.where(land, FILL, packed).astype("i2")


def file_name(var, month):
    return f"ocean_{var}_{month}.nc"


def write_file(path, var, month, day0, ndays):
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.title = "wxylib fixture shaped like BRAN2023 daily (synthetic values)"
    for name, size in [("Time", None), ("st_ocean", len(DEPTHS)),
                       ("yt_ocean", NY), ("xt_ocean", NX)]:
        ds.createDimension(name, size)
    t = ds.createVariable("Time", "f8", ("Time",))
    t.units = "days since 1979-01-01 00:00:00"
    t.calendar = "GREGORIAN"
    t[:] = 11323.5 + day0 + np.arange(ndays)
    z = ds.createVariable("st_ocean", "f8", ("st_ocean",))
    z.units = "meters"
    z.positive = "down"
    z[:] = DEPTHS
    y = ds.createVariable("yt_ocean", "f8", ("yt_ocean",))
    y.units = "degrees_N"
    y[:] = -50.0 + 0.5 * np.arange(NY) + 0.25
    x = ds.createVariable("xt_ocean", "f8", ("xt_ocean",))
    x.units = "degrees_E"
    x[:] = 140.0 + 0.5 * np.arange(NX) + 0.25
    v = ds.createVariable(var, "i2", ("Time", "st_ocean", "yt_ocean", "xt_ocean"),
                          zlib=True, complevel=1, shuffle=True, chunksizes=CHUNKS,
                          fill_value=np.int16(FILL))
    v.set_auto_maskandscale(False)
    for k, a in VARS[var].items():
        v.setncattr(k, a)
    skip = UNWRITTEN[var]
    for lt in range(ndays):
        g = day0 + lt
        for k in range(len(DEPTHS)):
            if (g, k, None, None) in skip:
                continue
            a = values(var, g, k)
            # write chunk by chunk, so a skipped chunk is never allocated
            for cy in range(NY // CHUNKS[2]):
                for cx in range(NX // CHUNKS[3]):
                    if (g, k, cy, cx) in skip:
                        continue
                    ys = slice(cy * CHUNKS[2], (cy + 1) * CHUNKS[2])
                    xs = slice(cx * CHUNKS[3], (cx + 1) * CHUNKS[3])
                    v[lt, k, ys, xs] = a[ys, xs]
    ds.close()


def scan(var, files):
    """One row per chunk of the global grid, from each file's HDF5 chunk
    index: the chunk's byte range in its file, or present = 0."""
    rows = []
    ntime = sum(n for _, n in MONTHS)
    grid = [ntime, len(DEPTHS), NY // CHUNKS[2], NX // CHUNKS[3]]
    day0 = 0
    for (month, ndays), path in zip(MONTHS, files):
        with h5py.File(path, "r") as h:
            dsid = h[var].id
            for lt in range(ndays):
                for k in range(grid[1]):
                    for cy in range(grid[2]):
                        for cx in range(grid[3]):
                            key = [day0 + lt, k, cy, cx]
                            coord = (lt, k, cy * CHUNKS[2], cx * CHUNKS[3])
                            try:
                                info = dsid.get_chunk_info_by_coord(coord)
                            except Exception:
                                info = None
                            if info is None or info.byte_offset is None or info.size == 0:
                                rows.append(dict(key=key, present=0, path=None, offset=None,
                                                 size=None, info=None))
                                continue
                            assert info.filter_mask == 0, (path, coord)
                            rows.append(dict(
                                key=key, present=1,
                                path=MIRROR + os.path.basename(path),
                                offset=int(info.byte_offset), size=int(info.size),
                                info="COMPRESSION=DEFLATE; FILTER=SHUFFLE; ENDIANNESS=LITTLE"))
        day0 += ndays
    rows.sort(key=lambda r: r["key"])
    return rows


def write_refs(var, rows):
    n = 4
    cols = {f"dim_{i}": pa.array([r["key"][i] for r in rows], pa.int64()) for i in range(n)}
    cols.update(
        present=pa.array([r["present"] for r in rows], pa.int8()),
        path=pa.array([r["path"] for r in rows], pa.string()),
        offset=pa.array([r["offset"] for r in rows], pa.int64()),
        size=pa.array([r["size"] for r in rows], pa.int64()),
        info=pa.array([r["info"] for r in rows], pa.string()),
        data=pa.array([None] * len(rows), pa.binary()),
    )
    pq.write_table(pa.table(cols), os.path.join(OUT, f"refs_{var}.parquet"),
                   compression="zstd")


def oracle(var, files, start, count):
    """What netCDF-C reads for the window, over the global time axis."""
    parts = []
    day0 = 0
    for (month, ndays), path in zip(MONTHS, files):
        lo, hi = max(start[0], day0), min(start[0] + count[0], day0 + ndays)
        if lo < hi:
            with netCDF4.Dataset(path) as ds:
                v = ds[var]
                v.set_auto_maskandscale(False)
                parts.append(np.asarray(v[
                    lo - day0:hi - day0,
                    start[1]:start[1] + count[1],
                    start[2]:start[2] + count[2],
                    start[3]:start[3] + count[3]]))
        day0 += ndays
    return np.concatenate(parts).astype("<i2")


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(NC)
    os.makedirs(os.path.join(OUT, "oracle"))
    files = {}
    for var in VARS:
        files[var] = []
        day0 = 0
        for month, ndays in MONTHS:
            p = os.path.join(NC, file_name(var, month))
            write_file(p, var, month, day0, ndays)
            t = int(datetime.datetime(2010, int(month[5:]), 1,
                                      tzinfo=datetime.timezone.utc).timestamp())
            os.utime(p, (t, t))
            files[var].append(p)
            day0 += ndays

    ntime = sum(n for _, n in MONTHS)
    shape = [ntime, len(DEPTHS), NY, NX]
    arrays = []
    counts = {}
    for var in VARS:
        rows = scan(var, files[var])
        write_refs(var, rows)
        counts[var] = dict(present=sum(r["present"] for r in rows),
                           absent=sum(1 - r["present"] for r in rows),
                           per_file=[sum(1 for r in rows if r["path"] == MIRROR + os.path.basename(f))
                                     for f in files[var]])
        arrays.append(dict(
            array=var, kind="data", storage="referenced",
            dims=["Time", "st_ocean", "yt_ocean", "xt_ocean"], shape=shape, chunks=CHUNKS,
            dtype="<i2", fill=str(FILL),
            codec=dict(compressor=dict(id="zlib", level=1),
                       filters=[dict(id="shuffle", elementsize=2)]),
            attrs=dict(VARS[var]), affine=None, values=None))
    with netCDF4.Dataset(files["temp"][0]) as ds:
        tunits = ds["Time"].units
        tcal = ds["Time"].calendar
    time = [11323.5 + d for d in range(ntime)]
    arrays += [
        dict(array="Time", kind="coord", storage="affine", dims=["Time"], shape=[ntime],
             chunks=[ntime], dtype="<f8", fill=None, codec=None,
             attrs=dict(units=tunits, calendar=tcal),
             affine=dict(start=time[0], step=1.0), values=None),
        dict(array="st_ocean", kind="coord", storage="inline", dims=["st_ocean"],
             shape=[len(DEPTHS)], chunks=[len(DEPTHS)], dtype="<f8", fill=None, codec=None,
             attrs=dict(units="meters", positive="down"), affine=None, values=DEPTHS),
        dict(array="yt_ocean", kind="coord", storage="affine", dims=["yt_ocean"], shape=[NY],
             chunks=[NY], dtype="<f8", fill=None, codec=None, attrs=dict(units="degrees_N"),
             affine=dict(start=-49.75, step=0.5), values=None),
        dict(array="xt_ocean", kind="coord", storage="affine", dims=["xt_ocean"], shape=[NX],
             chunks=[NX], dtype="<f8", fill=None, codec=None, attrs=dict(units="degrees_E"),
             affine=dict(start=140.25, step=0.5), values=None),
    ]
    arrays.sort(key=lambda a: a["array"])
    with open(os.path.join(OUT, "arrays.json"), "w") as f:
        json.dump(dict(format_version=FORMAT, arrays=arrays,
                       remap=[dict(target="fixture", prefix=MIRROR, replacement="nc/")]),
                  f, indent=1)

    windows = []
    for i, (var, start, count) in enumerate(WINDOWS):
        v = oracle(var, files[var], start, count)
        assert list(v.shape) == count, (var, v.shape, count)
        fname = f"oracle/{var}-w{i}.bin"
        v.tofile(os.path.join(OUT, fname))
        windows.append(dict(array=var, start=start, count=count, file=fname,
                            fill_cells=int((v == FILL).sum())))

    manifest = dict(
        oracle="netCDF-C via netCDF4-python (mask and scale off); byte ranges from HDF5",
        netCDF4=netCDF4.__version__, netcdf_c=netCDF4.__netcdf4libversion__,
        h5py=h5py.__version__, hdf5=h5py.version.hdf5_version,
        format_version=FORMAT, remap_target="fixture",
        sources={var: [os.path.basename(f) for f in fs] for var, fs in files.items()},
        arrays={var: dict(shape=shape, chunks=CHUNKS, dtype="int16", fill=FILL,
                          scale_factor=VARS[var]["scale_factor"],
                          add_offset=VARS[var]["add_offset"], **counts[var])
                for var in VARS},
        unwritten={var: [list(u) for u in us] for var, us in UNWRITTEN.items()},
        windows=windows,
    )
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)


if __name__ == "__main__":
    main()
