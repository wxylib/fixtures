"""oraclecase: ask a reference library the questions in a fixture set, keep
its answers, and compare answer sets (format: ORACLE.md, design m27).

A set is a directory with a cases.json (the questions). Each oracle run
writes answers/<oracle-id>/answers.json plus one file per array or feature
answer, so GDAL 3.8, 3.12 and 3.14dev (or another library) can answer the
same questions side by side.

    python scripts/oraclecase.py record SET [--id ORACLE_ID] [--case ID ...]
    python scripts/oraclecase.py diff SET ORACLE_A ORACLE_B
    python scripts/oraclecase.py adopt-cog          # cog/manifest.json -> cog/cases.json
    python scripts/oraclecase.py adopt-vrt          # vrt/manifest.json -> vrt/cases.json
    python scripts/oraclecase.py list SET

record needs osgeo.gdal and numpy (GDAL 3.8 or later); diff, adopt and
list need only the standard library and numpy. Content is ASCII only.
"""

import argparse
import contextlib
import datetime
import hashlib
import json
import os
import platform
import sys

FORMAT = "oracle-case/0"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# numpy names are the format's type names; GDAL names map onto them.
GDAL_TYPES = {
    "Byte": "uint8", "Int8": "int8", "UInt16": "uint16", "Int16": "int16",
    "UInt32": "uint32", "Int32": "int32", "UInt64": "uint64", "Int64": "int64",
    "Float32": "float32", "Float64": "float64",
}


def set_dir(name):
    d = name if os.path.isabs(name) else os.path.join(ROOT, name)
    if not os.path.exists(os.path.join(d, "cases.json")):
        sys.exit("no cases.json in %s" % d)
    return d


def load(path):
    with open(path) as f:
        return json.load(f)


def dump(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=False)
        f.write("\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def resolve_dsn(dsn, sdir):
    """{set} is the set directory, absolute. Everything else is passed as
    written: the recorder runs with the set directory as working directory,
    so relative paths inside vrt://, NETCDF:"..." and friends resolve too."""
    return dsn.replace("{set}", sdir)


# ---------------------------------------------------------------- record

class Messages:
    """Collects CPLError output for one case, as the answer's messages."""

    LEVELS = {0: "none", 1: "debug", 2: "warning", 3: "failure", 4: "fatal"}

    def __init__(self):
        self.items = []

    def __call__(self, cls, code, msg):
        self.items.append({"level": self.LEVELS.get(cls, str(cls)), "code": code,
                           "msg": msg.strip()})


@contextlib.contextmanager
def config(gdal, opts):
    old = {k: gdal.GetConfigOption(k) for k in opts}
    for k, v in opts.items():
        gdal.SetConfigOption(k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            gdal.SetConfigOption(k, v)


def write_array(arr, sdir, rel):
    import numpy as np
    path = os.path.join(sdir, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    a = np.ascontiguousarray(arr)
    a.astype(a.dtype.newbyteorder("<"), copy=False).tofile(path)
    return {"file": rel, "dtype": a.dtype.name, "shape": list(a.shape), "order": "C",
            "endian": "little", "sha256": sha256(path)}


def crs_of(srs):
    if srs is None:
        return None
    out = {"wkt": srs.ExportToWkt()}
    name, code = srs.GetAuthorityName(None), srs.GetAuthorityCode(None)
    if name and code:
        out["code"] = "%s:%s" % (name, code)
    return out


def raster_info(gdal, ds):
    bands = []
    for i in range(1, ds.RasterCount + 1):
        b = ds.GetRasterBand(i)
        bx, by = b.GetBlockSize()
        bands.append({
            "dtype": GDAL_TYPES.get(gdal.GetDataTypeName(b.DataType),
                                    gdal.GetDataTypeName(b.DataType)),
            "nodata": b.GetNoDataValue(), "block": [bx, by],
            "overviews": [[b.GetOverview(k).XSize, b.GetOverview(k).YSize]
                          for k in range(b.GetOverviewCount())],
            "scale": b.GetScale(), "offset": b.GetOffset(),
        })
    gt = ds.GetGeoTransform(can_return_null=True)
    return {"size": [ds.RasterXSize, ds.RasterYSize], "bands": bands,
            "geotransform": list(gt) if gt else None, "crs": crs_of(ds.GetSpatialRef())}


RESAMPLING = {"nearest": 0, "bilinear": 1, "cubic": 2, "cubicspline": 3, "lanczos": 4,
              "average": 5, "mode": 6, "gauss": 7, "rms": 14}  # GDALRIOResampleAlg


def op_raster(gdal, case, sdir, rel, ans):
    ds = gdal.OpenEx(resolve_dsn(case["dsn"], sdir), gdal.OF_RASTER,
                     open_options=case.get("open_options", []))
    if case["op"] == "raster.info":
        ans["info"] = raster_info(gdal, ds)
        return
    w = case.get("window") or [0, 0, ds.RasterXSize, ds.RasterYSize]
    buf = case.get("buf") or w[2:]
    bands = case.get("bands") or list(range(1, ds.RasterCount + 1))
    alg = RESAMPLING[case.get("resampling", "nearest")]
    import numpy as np
    arrs = [ds.GetRasterBand(b).ReadAsArray(w[0], w[1], w[2], w[3], buf[0], buf[1],
                                            resample_alg=alg) for b in bands]
    ans["values"] = write_array(np.stack(arrs), sdir, rel + ".bin")


def op_mdim(gdal, case, sdir, rel, ans):
    ds = gdal.OpenEx(resolve_dsn(case["dsn"], sdir), gdal.OF_MULTIDIM_RASTER,
                     open_options=case.get("open_options", []))
    arr = ds.GetRootGroup().OpenMDArrayFromFullname(case["array"])
    if arr is None:
        raise RuntimeError("no array %s" % case["array"])
    dims = arr.GetDimensions()
    ans["info"] = {
        "dims": [{"name": d.GetName(), "size": d.GetSize()} for d in dims],
        "dtype": GDAL_TYPES.get(gdal.GetDataTypeName(arr.GetDataType().GetNumericDataType()), None),
        "nodata": arr.GetNoDataValueAsDouble() if arr.GetNoDataValueAsRaw() is not None else None,
        "scale": arr.GetScale(), "offset": arr.GetOffset(), "crs": crs_of(arr.GetSpatialRef()),
    }
    if case["op"] == "mdim.info":
        return
    n = len(dims)
    start = case.get("start") or [0] * n
    step = case.get("step") or [1] * n
    count = case.get("count") or [(d.GetSize() - s + st - 1) // st
                                  for d, s, st in zip(dims, start, step)]
    a = arr.ReadAsArray(array_start_idx=start, count=count, array_step=step)
    ans["values"] = write_array(a, sdir, rel + ".bin")


def op_vector(gdal, case, sdir, rel, ans):
    ds = gdal.OpenEx(resolve_dsn(case["dsn"], sdir), gdal.OF_VECTOR,
                     open_options=case.get("open_options", []))
    names = [case["layer"]] if case.get("layer") else \
        [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
    layers = {}
    for name in names:
        lyr = ds.GetLayerByName(name)
        if lyr is None:
            raise RuntimeError("no layer %s" % name)
        if case.get("bbox"):
            lyr.SetSpatialFilterRect(*case["bbox"])
        if case.get("where"):
            lyr.SetAttributeFilter(case["where"])
        defn = lyr.GetLayerDefn()
        fields = [{"name": defn.GetFieldDefn(i).GetName(),
                   "type": defn.GetFieldDefn(i).GetTypeName()} for i in range(defn.GetFieldCount())]
        feats = []
        lyr.ResetReading()
        for f in lyr:
            g = f.GetGeometryRef()
            feats.append({"fid": f.GetFID(),
                          "properties": {fd["name"]: f.GetField(fd["name"]) for fd in fields},
                          "wkt": g.ExportToIsoWkt() if g is not None else None})
        layers[name] = {"fields": fields, "crs": crs_of(lyr.GetSpatialRef()),
                        "geometry_type": gdal.ogr.GeometryTypeToName(lyr.GetGeomType()),
                        "features": feats}
    path = os.path.join(sdir, rel + ".json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    dump({"layers": layers}, path)
    ans["features"] = {"file": rel + ".json", "sha256": sha256(path)}


def op_warp_grid(gdal, case, sdir, rel, ans):
    """The target grid gdalwarp settles on: a warp to VRT, so no pixels move."""
    src = gdal.OpenEx(resolve_dsn(case["dsn"], sdir), gdal.OF_RASTER,
                      open_options=case.get("open_options", []))
    out = gdal.Warp("", src, options=["-of", "VRT"] + list(case.get("args", [])))
    gt = out.GetGeoTransform()
    nx, ny = out.RasterXSize, out.RasterYSize
    ans["grid"] = {"geotransform": list(gt), "size": [nx, ny],
                   "extent": [gt[0], gt[3] + ny * gt[5], gt[0] + nx * gt[1], gt[3]],
                   "crs": crs_of(out.GetSpatialRef())}


OPS = {"raster.read": op_raster, "raster.info": op_raster, "warp.grid": op_warp_grid,
       "mdim.read": op_mdim, "mdim.info": op_mdim, "vector.read": op_vector}


def record(args):
    from osgeo import gdal, ogr
    import numpy as np
    gdal.ogr = ogr
    gdal.UseExceptions()
    sdir = set_dir(args.set)
    cases = load(os.path.join(sdir, "cases.json"))
    oid = args.id or "gdal-" + gdal.__version__
    adir = os.path.join("answers", oid)
    out = os.path.join(sdir, adir, "answers.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    prev = load(out)["answers"] if os.path.exists(out) and args.case else {}
    answers = dict(prev)
    here = os.getcwd()
    os.chdir(sdir)
    try:
        for case in cases["cases"]:
            if args.case and case["id"] not in args.case:
                continue
            fn = OPS.get(case["op"])
            ans = {"status": "ok"}
            msgs = Messages()
            gdal.PushErrorHandler(msgs)
            try:
                if fn is None:
                    raise RuntimeError("op %s not known to this recorder" % case["op"])
                with config(gdal, case.get("config", {})):
                    fn(gdal, case, sdir, os.path.join(adir, case["id"].replace("/", "__")), ans)
            except Exception as e:  # GDAL's refusal is an answer too
                ans = {"status": "error", "error": str(e).strip()}
            finally:
                gdal.PopErrorHandler()
            if msgs.items:
                ans["messages"] = msgs.items
            answers[case["id"]] = ans
            print("%-40s %s" % (case["id"], ans["status"]), file=sys.stderr)
    finally:
        os.chdir(here)
    dump({"format": FORMAT, "set": cases.get("set"),
          "oracle": {"name": "gdal", "version": gdal.__version__,
                     "release": gdal.VersionInfo("RELEASE_NAME"),
                     "via": "osgeo.gdal", "python": platform.python_version(),
                     "numpy": np.__version__, "proj": "%d.%d.%d" % _proj()},
          "recorded": datetime.date.today().isoformat(),
          "answers": answers}, out)
    print("wrote %s" % os.path.relpath(out, ROOT), file=sys.stderr)


def _proj():
    from osgeo import osr
    return (osr.GetPROJVersionMajor(), osr.GetPROJVersionMinor(), osr.GetPROJVersionMicro())


# ---------------------------------------------------------------- diff

def read_values(sdir, v):
    import numpy as np
    return np.fromfile(os.path.join(sdir, v["file"]),
                       dtype=np.dtype(v["dtype"]).newbyteorder("<")).reshape(v["shape"])


def diff(args):
    """Compares two answer sets case by case: the same bytes, the same
    values within the case's tolerance, or how far apart."""
    import numpy as np
    sdir = set_dir(args.set)
    cases = {c["id"]: c for c in load(os.path.join(sdir, "cases.json"))["cases"]}
    a = load(os.path.join(sdir, "answers", args.a, "answers.json"))["answers"]
    b = load(os.path.join(sdir, "answers", args.b, "answers.json"))["answers"]
    same = changed = unanswered = 0
    for cid, case in cases.items():
        x, y = a.get(cid), b.get(cid)
        if x is None or y is None:
            print("%-40s not answered by %s" % (cid, args.a if x is None else args.b))
            unanswered += 1
            continue
        if x["status"] != y["status"]:
            print("%-40s status %s -> %s" % (cid, x["status"], y["status"]))
            changed += 1
            continue
        line = "same"
        for k in ("values", "features"):
            if k in x and x[k]["sha256"] != y[k]["sha256"]:
                if k == "values":
                    p, q = read_values(sdir, x[k]), read_values(sdir, y[k])
                    if p.shape != q.shape:
                        line = "shape %s -> %s" % (p.shape, q.shape)
                    else:
                        d = np.abs(p.astype("f8") - q.astype("f8"))
                        both_nan = np.isnan(p.astype("f8")) & np.isnan(q.astype("f8"))
                        d[both_nan] = 0
                        line = "values: %d of %d differ, max %g" % ((d > 0).sum(), d.size, d.max())
                else:
                    line = "features differ"
        if "grid" in x and line == "same":
            line = grid_diff(x["grid"], y["grid"])
        xi, yi = comparable(x.get("info")), comparable(y.get("info"))
        if xi != yi and line == "same":
            line = "info differs: " + ", ".join(sorted(
                k for k in set(xi) | set(yi) if xi.get(k) != yi.get(k)))
        if line == "same":
            same += 1
        else:
            changed += 1
            print("%-40s %s" % (cid, line))
    print("%d same, %d differ, %d not answered by both (%s vs %s)"
          % (same, changed, unanswered, args.a, args.b))


def grid_diff(p, q):
    """Two warp.grid answers: same size and geotransform to 1e-9 relative
    (of the larger cell size, for the origin), and the same CRS code."""
    if p["size"] != q["size"]:
        return "grid size %s -> %s" % (p["size"], q["size"])
    cell = max(abs(p["geotransform"][1]), abs(p["geotransform"][5]), 1e-300)
    for i, (u, v) in enumerate(zip(p["geotransform"], q["geotransform"])):
        scale = cell if i in (0, 3) else max(abs(u), abs(v), 1e-300)
        if abs(u - v) > 1e-9 * scale:
            return "grid geotransform[%d] %r -> %r" % (i, u, v)
    if comparable({"crs": p["crs"]}) != comparable({"crs": q["crs"]}):
        return "grid crs differs"
    return "same"


def comparable(info):
    """Info as compared across oracle versions: a CRS by its code when it
    has one, since WKT spelling changes between GDAL and PROJ releases."""
    info = dict(info or {})
    crs = info.get("crs")
    if isinstance(crs, dict) and crs.get("code"):
        info["crs"] = crs["code"]
    return info


# ---------------------------------------------------------------- adopt

def adopt_cog(args):
    """Describes the existing cog/ set in the format without touching its
    files: its windows become raster.read cases, its oracle/*.bin files the
    answers of the run that made them (rasterio's GDAL, per manifest)."""
    sdir = os.path.join(ROOT, "cog")
    m = load(os.path.join(sdir, "manifest.json"))
    cases, answers = [], {}
    for f in m["cases"]:
        stem = f["file"][:-4]
        cases.append({"id": stem + "/info", "op": "raster.info", "dsn": f["file"]})
        for k, w in enumerate(f["windows"]):
            cid = "%s/w%d" % (stem, k)
            cases.append({"id": cid, "op": "raster.read", "dsn": f["file"],
                          "window": [w["col"], w["row"], w["width"], w["height"]],
                          "compare": {"kind": "exact"}})
            answers[cid] = {"status": "ok", "values": {
                "file": w["file"], "dtype": f["dtype"],
                "shape": [f["count"], w["height"], w["width"]], "order": "C",
                "endian": "little", "sha256": sha256(os.path.join(sdir, w["file"]))}}
    dump({"format": FORMAT, "set": "cog",
          "description": "Window reads of the five COG fixtures (see README, cog/).",
          "cases": cases}, os.path.join(sdir, "cases.json"))
    oid = "gdal-%s" % m["gdal"]
    os.makedirs(os.path.join(sdir, "answers", oid), exist_ok=True)
    dump({"format": FORMAT, "set": "cog",
          "oracle": {"name": "gdal", "version": m["gdal"], "via": "rasterio %s" % m["rasterio"]},
          "recorded": "adopted from manifest.json", "answers": answers},
         os.path.join(sdir, "answers", oid, "answers.json"))
    print("cog: %d cases, %d answers adopted as %s" % (len(cases), len(answers), oid))


def adopt_vrt(args):
    """Describes a vrt/ corpus written as manifest.json (dsn, requests with
    window, buf and file) in the format, in place."""
    sdir = os.path.join(ROOT, args.dir)
    m = load(os.path.join(sdir, "manifest.json"))
    cases, answers = [], {}
    for c in m["cases"]:
        stem = c["name"][:-4] if c["name"].endswith(".vrt") else c["name"]
        cases.append({"id": stem + "/info", "op": "raster.info", "dsn": c["dsn"],
                      "note": c.get("desc")})
        for k, r in enumerate(c["requests"]):
            cid = "%s/r%d" % (stem, k)
            case = {"id": cid, "op": "raster.read", "dsn": c["dsn"],
                    "window": r["window"], "buf": r["buf"], "resampling": "nearest",
                    "compare": {"kind": "exact"}}
            cases.append(case)
            answers[cid] = {"status": "ok", "values": {
                "file": r["file"], "dtype": GDAL_TYPES.get(c["dtype"], c["dtype"]),
                "shape": r["shape"], "order": "C", "endian": "little",
                "sha256": sha256(os.path.join(sdir, r["file"]))}}
    dump({"format": FORMAT, "set": os.path.basename(sdir),
          "description": "VRT conformance corpus (design m25).", "cases": cases},
         os.path.join(sdir, "cases.json"))
    oid = "gdal-%s" % m["gdal"]
    os.makedirs(os.path.join(sdir, "answers", oid), exist_ok=True)
    dump({"format": FORMAT, "set": os.path.basename(sdir),
          "oracle": {"name": "gdal", "version": m["gdal"], "via": m.get("backend")},
          "recorded": "adopted from manifest.json", "answers": answers},
         os.path.join(sdir, "answers", oid, "answers.json"))
    print("%s: %d cases, %d answers adopted as %s" % (args.dir, len(cases), len(answers), oid))


def list_set(args):
    sdir = set_dir(args.set)
    cases = load(os.path.join(sdir, "cases.json"))["cases"]
    adir = os.path.join(sdir, "answers")
    runs = sorted(os.listdir(adir)) if os.path.isdir(adir) else []
    print("%d cases; answer sets: %s" % (len(cases), ", ".join(runs) or "none"))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("record")
    r.add_argument("set")
    r.add_argument("--id", help="oracle id, default gdal-<version>")
    r.add_argument("--case", action="append", help="only these case ids")
    r.set_defaults(fn=record)
    d = sub.add_parser("diff")
    d.add_argument("set")
    d.add_argument("a")
    d.add_argument("b")
    d.set_defaults(fn=diff)
    a = sub.add_parser("adopt-cog")
    a.set_defaults(fn=adopt_cog)
    v = sub.add_parser("adopt-vrt")
    v.add_argument("--dir", default="vrt")
    v.set_defaults(fn=adopt_vrt)
    ls = sub.add_parser("list")
    ls.add_argument("set")
    ls.set_defaults(fn=list_set)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
