"""Zarr v3 fixtures made from the COG fixtures, so a Zarr reader can be
checked against the same GDAL oracle windows (cog/oracle/*.bin).

Georeferencing follows the GeoZarr conventions (zarr-conventions/spatial and
zarr-conventions/proj, v0.1) or plain CF-style 1D coordinate arrays, one
case each. Run from the repository root after make_cog_fixtures.py:

    python scripts/make_zarr_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import rasterio
import zarr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COG = os.path.join(ROOT, "cog")
OUT = os.path.join(ROOT, "zarr")

SPATIAL = {
    "schema_url": "https://raw.githubusercontent.com/zarr-conventions/spatial/refs/tags/v0.1/schema.json",
    "spec_url": "https://github.com/zarr-conventions/spatial/blob/v0.1/README.md",
    "uuid": "689b58e2-cf7b-45e0-9fff-9cfc0883d6b4",
    "name": "spatial",
    "description": "Spatial coordinate information",
}
PROJ = {
    "schema_url": "https://raw.githubusercontent.com/zarr-conventions/proj/refs/tags/v0.1/schema.json",
    "spec_url": "https://github.com/zarr-conventions/proj/blob/v0.1/README.md",
    "uuid": "f17cb550-5864-4468-aeb7-f3180cfb622f",
    "name": "proj",
    "description": "Coordinate reference system information for geospatial data",
}


def spatial_transform(gt):
    """GDAL geotransform to spatial:transform (Affine order a, b, c, d, e, f)."""
    return [gt[1], gt[2], gt[0], gt[4], gt[5], gt[3]]


# Each case: source COG, how to georeference, layout.
CASES = [
    {
        # GeoZarr spatial + proj on the group (inherited); dimensions on the
        # array. Chunks do not line up with the COG's 128 tiles.
        "store": "world-int16-spatial.zarr",
        "source": "world-int16-deflate.tif",
        "georef": "spatial",
        "registration": "pixel",
        "where": "group",
        "dims": ["y", "x"],
        "chunks": [100, 100],
    },
    {
        # No conventions: CF-style 1D centre coordinates stored as float32,
        # which a reader must normalise to recover the affine.
        "store": "polar-f32-cf.zarr",
        "source": "polar-f32-zstd-bigtiff.tif",
        "georef": "cf",
        "dims": ["y", "x"],
        "chunks": [64, 96],
    },
    {
        # A non-spatial leading dim, node registration (transform gives cell
        # centres), conventions on the array itself.
        "store": "utm-u16-3band-node.zarr",
        "source": "utm-u16-lzw-pred2-3band.tif",
        "georef": "spatial",
        "registration": "node",
        "where": "array",
        "dims": ["band", "y", "x"],
        "chunks": [1, 64, 64],
    },
]


def main():
    cog = json.load(open(os.path.join(COG, "manifest.json")))
    by_file = {c["file"]: c for c in cog["cases"]}
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    out = {"gdal": cog["gdal"], "rasterio": cog["rasterio"], "zarr": zarr.__version__, "cases": []}

    for case in CASES:
        src = by_file[case["source"]]
        with rasterio.open(os.path.join(COG, case["source"])) as ds:
            data = ds.read()
            gt = list(ds.transform.to_gdal())
            epsg = ds.crs.to_epsg()
            nodata = ds.nodata
        if case["dims"][0] != "band":
            data = data[0]

        path = os.path.join(OUT, case["store"])
        group = zarr.open_group(path, mode="w", zarr_format=3)
        fill = nodata if nodata is not None else 0
        arr = group.create_array(
            "data",
            shape=data.shape,
            chunks=case["chunks"],
            dtype=data.dtype,
            fill_value=fill,
            dimension_names=case["dims"],
        )
        arr[...] = data

        h, w = data.shape[-2:]
        if case["georef"] == "spatial":
            t = spatial_transform(gt)
            if case["registration"] == "node":
                # move the origin from the outer corner to the first centre
                t[2] += 0.5 * (t[0] + t[1])
                t[5] += 0.5 * (t[3] + t[4])
            attrs = {
                "zarr_conventions": [SPATIAL, PROJ],
                "proj:code": "EPSG:%d" % epsg,
                "spatial:transform": t,
                "spatial:shape": [h, w],
                "spatial:registration": case["registration"],
            }
            dims = {"spatial:dimensions": ["y", "x"]}
            if case["where"] == "group":
                group.attrs.update(attrs)
                arr.attrs.update(dims)
            else:
                arr.attrs.update({**attrs, **dims})
        else:
            xc = gt[0] + gt[1] * (np.arange(w) + 0.5)
            yc = gt[3] + gt[5] * (np.arange(h) + 0.5)
            for name, v in (("x", xc), ("y", yc)):
                c = group.create_array(name, shape=v.shape, dtype="float32", dimension_names=[name])
                c[...] = v.astype("float32")
            group.attrs.update({"zarr_conventions": [PROJ], "proj:code": "EPSG:%d" % epsg})
        if nodata is not None:
            arr.attrs["_FillValue"] = nodata

        out["cases"].append(
            {
                "store": case["store"],
                "array": "data",
                "source": case["source"],
                "georef": case["georef"],
                "dims": case["dims"],
                "shape": list(data.shape),
                "chunks": case["chunks"],
                "dtype": str(data.dtype),
                "geotransform": gt,
                "epsg": epsg,
                "nodata": nodata,
                # same pixels as the COG, so the COG's GDAL windows apply
                "windows": [dict(wd, file="../cog/" + wd["file"]) for wd in src["windows"]],
            }
        )

    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
