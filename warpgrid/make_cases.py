"""Writes warpgrid/cases.json: the target grids gdalwarp infers (design m28).

    python warpgrid/make_cases.py

Sources are inline VRT documents with no pixels (a grid and a CRS only).
"""

import json
import os

SOURCES = {
    "world-0.5": ("EPSG:4326", [-180, 180, -90, 90], [720, 360]),
    "world-0.25-0to360": ("+proj=longlat +datum=WGS84 +lon_wrap=180 +no_defs", [0, 360, -90, 90], [1440, 720]),
    "antarctic-25km": ("EPSG:3031", [-6e6, 6e6, -6e6, 6e6], [480, 480]),
    "utm55s-10m": ("EPSG:32755", [499980, 609780, 5190220, 5300020], [10980, 10980]),
    "pacific-150to210": ("EPSG:4326", [150, 210, -60, -30], [240, 120]),
}
TARGETS = {
    "4326": "EPSG:4326",
    "3857": "EPSG:3857",
    "3031": "EPSG:3031",
    "laea-hobart": "+proj=laea +lat_0=-42.9 +lon_0=147.3 +datum=WGS84",
    "moll": "+proj=moll +datum=WGS84",
    "32755": "EPSG:32755",
}


def vrt(crs, ext, size):
    xmin, xmax, ymin, ymax = ext
    gt = [xmin, (xmax - xmin) / size[0], 0, ymax, 0, -(ymax - ymin) / size[1]]
    return ('<VRTDataset rasterXSize="%d" rasterYSize="%d"><SRS>%s</SRS>'
            '<GeoTransform>%s</GeoTransform>'
            '<VRTRasterBand dataType="Byte" band="1"/></VRTDataset>'
            % (size[0], size[1], crs, ", ".join(repr(float(v)) for v in gt)))


cases = []
for s, (crs, ext, size) in SOURCES.items():
    for t, tcrs in TARGETS.items():
        if s.startswith("world") and t == "32755":
            continue
        cases.append({"id": "suggest/%s/%s" % (s, t), "op": "warp.grid",
                      "dsn": vrt(crs, ext, size), "args": ["-t_srs", tcrs]})

MATRIX = {
    "utm55s-10m/4326": ("utm55s-10m", "EPSG:4326", [146.99, -43.42, 148.36, -42.4], 0.0001),
    "world-0.5/3031": ("world-0.5", "EPSG:3031", [-4e6, -4e6, 4e6, 4e6], 25000),
}
for key, (s, tcrs, te, tr) in MATRIX.items():
    src = vrt(*SOURCES[s])
    t = lambda v: repr(v)  # noqa: E731
    te_ = ["-te"] + [t(v) for v in te]
    tr1 = ["-tr", t(tr), t(tr)]
    tr2 = ["-tr", t(tr * 1.37), t(tr * 1.37)]
    rows = {
        "te": te_,
        "te-tr": te_ + tr1,
        "te-tr-not-multiple": te_ + tr2,
        "te-tr-tap": te_ + tr2 + ["-tap"],
        "te-ts": te_ + ["-ts", "1000", "1000"],
        "te-ts-w0": te_ + ["-ts", "1000", "0"],
        "ts": ["-ts", "1000", "1000"],
        "ts-w0": ["-ts", "1000", "0"],
        "tr-ts": tr1 + ["-ts", "1000", "1000"],
        "tap-alone": ["-tap"],
        "te_srs": ["-te_srs", "EPSG:4326", "-te", "140", "-70", "160", "-60"],
    }
    # -tr without -te over the whole-world source into 3031 asks for about
    # 1e39 cells, and GDAL 3.8.4 aborts the process (std::length_error):
    # not a case a recorder can survive, so only the UTM pair has it.
    if s == "utm55s-10m":
        rows["tr"] = tr1
        rows["tr-tap"] = tr1 + ["-tap"]
    for name, args in rows.items():
        cases.append({"id": "matrix/%s/%s" % (key, name), "op": "warp.grid",
                      "dsn": src, "args": ["-t_srs", tcrs] + args})

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "cases.json"), "w") as f:
    json.dump({"format": "oracle-case/0", "set": "warpgrid",
               "description": "The target grid gdalwarp infers from what it is given "
                              "(-t_srs, -te, -te_srs, -tr, -ts, -tap), over five source "
                              "grids and six target CRSs. wxylib/design m28.",
               "cases": cases}, f, indent=1)
    f.write("\n")
print(len(cases), "cases")
