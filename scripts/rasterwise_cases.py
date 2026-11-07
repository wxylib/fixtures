#!/usr/bin/env python3
"""Write rasterwise/cases.json: oracle-case/0 questions over the published
files of the rgdal-dev/rasterwise problem-file archive.

usage: python3 scripts/rasterwise_cases.py <rasterwise checkout>

Questions are keyed by the file's sha256 (first 12 hex digits) and read the
published bytes in place, through /vsicurl/ on the manifest's url. Nothing
is copied into this repo except small answers. Discovery (which subdatasets
and arrays a file has) uses osgeo.gdal, so the question list reflects what
the discovering GDAL could see; files it could not open still get an open
question, whose answer is the refusal.

Per data file:
  <id>/open           raster.info on the file itself
  <id>/sds/<name>     raster.info per classic subdataset
  <id>/sds/<name>/r   raster.read of band 1, whole if at most 16384 cells,
                      else the centred 128 x 128 window
  <id>/md/<array>     mdim.info per multidimensional array
  <id>/md/<array>/r   mdim.read of every 1D and 2D array of at most 262144
                      elements (coordinates are the evidence for assertions)

Design: wxylib/design docs/m30-legacy-corpus.md (rungs R1 and R2).
"""
import csv
import json
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SET = os.path.join(ROOT, "rasterwise")
DATA_FORMATS = {"netcdf3", "netcdf4/hdf5", "hdf4", "geotiff", "tiff", "binary", "vrt", "grib"}
READ_CELLS = 16384
WIN = 128
MD_READ = 262144


def clean(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "root"


def main():
    rw = sys.argv[1]
    from osgeo import gdal
    gdal.UseExceptions()
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    commit = subprocess.run(["git", "-C", rw, "rev-parse", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    rows = [r for r in csv.DictReader(open(os.path.join(rw, "manifest.csv")))
            if r["url"] and r["format"] in DATA_FORMATS]
    cases, files = [], []
    for r in sorted(rows, key=lambda r: r["canonical_path"]):
        fid = r["sha256"][:12]
        dsn = "/vsicurl/" + r["url"]
        note = "%s (%s)" % (r["canonical_path"], r["case"] or "no case")
        files.append({"id": fid, "sha256": r["sha256"], "size": int(r["size"]),
                      "format": r["format"], "canonical_path": r["canonical_path"],
                      "case": r["case"], "url": r["url"]})
        cases.append({"id": fid + "/open", "op": "raster.info", "dsn": dsn, "note": note})
        try:
            ds = gdal.Open(dsn)
        except Exception:
            ds = None
        if ds is not None:
            sds = ds.GetMetadata("SUBDATASETS") or {}
            names = [sds[k] for k in sorted(sds) if k.endswith("_NAME")]
            if not names and ds.RasterCount:
                names = [None]
            for nm in names:
                sid = fid + "/sds/" + (clean(nm.split(":")[-1]) if nm else "main")
                sdsn = nm or dsn
                if nm:
                    cases.append({"id": sid, "op": "raster.info", "dsn": sdsn})
                try:
                    s = gdal.Open(sdsn)
                    nx, ny = s.RasterXSize, s.RasterYSize
                except Exception:
                    continue
                if nx * ny <= READ_CELLS:
                    win = [0, 0, nx, ny]
                else:
                    w, h = min(WIN, nx), min(WIN, ny)
                    win = [(nx - w) // 2, (ny - h) // 2, w, h]
                cases.append({"id": sid + "/r", "op": "raster.read", "dsn": sdsn,
                              "window": win, "bands": [1]})
        try:
            md = gdal.OpenEx(dsn, gdal.OF_MULTIDIM_RASTER)
        except Exception:
            md = None
        if md is not None:
            stack = [md.GetRootGroup()]
            while stack:
                g = stack.pop()
                for an in g.GetMDArrayNames() or []:
                    a = g.OpenMDArray(an)
                    full = a.GetFullName()
                    mid = fid + "/md/" + clean(full)
                    cases.append({"id": mid, "op": "mdim.info", "dsn": dsn, "array": full})
                    dims = a.GetDimensions()
                    n = 1
                    for d in dims:
                        n *= d.GetSize()
                    numeric = a.GetDataType().GetClass() == gdal.GEDTC_NUMERIC
                    if numeric and 1 <= len(dims) <= 2 and 1 <= n <= MD_READ:
                        cases.append({"id": mid + "/r", "op": "mdim.read", "dsn": dsn,
                                      "array": full})
                for gn in g.GetGroupNames() or []:
                    stack.append(g.OpenGroup(gn))
        print("%-14s %4d %s" % (fid, sum(c["id"].startswith(fid) for c in cases),
                                r["canonical_path"]), file=sys.stderr)
    os.makedirs(SET, exist_ok=True)
    doc = {"format": "oracle-case/0", "set": "rasterwise",
           "description": "Published files of the rgdal-dev/rasterwise problem-file archive "
                          "(manifest at %s), read in place over /vsicurl/. Questions found by "
                          "GDAL %s." % (commit[:7], gdal.__version__),
           "source": {"repo": "rgdal-dev/rasterwise", "commit": commit},
           "cases": cases}
    with open(os.path.join(SET, "cases.json"), "w") as f:
        json.dump(doc, f, indent=1)
        f.write("\n")
    with open(os.path.join(SET, "files.json"), "w") as f:
        json.dump({"source": doc["source"], "files": files}, f, indent=1)
        f.write("\n")
    print("%d files, %d cases" % (len(files), len(cases)), file=sys.stderr)


if __name__ == "__main__":
    main()
