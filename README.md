# wxylib/fixtures

Test data and GDAL oracle outputs for wxylib, kept out of
[wxylib/core](https://github.com/wxylib/core) so it stays small.

## cog/

Small Cloud Optimized GeoTIFFs with synthetic pixel values, made by
`scripts/make_cog_fixtures.py`, and what GDAL (through rasterio) reads for a
few windows of each. `cog/manifest.json` records the GDAL version, each
file's size, blocks, overviews, geotransform (GDAL order), EPSG code and
nodata, and the windows. Each window's oracle values are in
`cog/oracle/*.bin`: little-endian, C order, band-major (band, row, col).

| File | Size | Type | Compression | Notes |
|------|------|------|-------------|-------|
| world-int16-deflate.tif | 720 x 360 | int16 | DEFLATE | EPSG:4326, 3 overviews |
| polar-f32-zstd-bigtiff.tif | 500 x 400 | float32 | ZSTD | EPSG:3031, BigTIFF, nodata -9999 |
| utm-u16-lzw-pred2-3band.tif | 300 x 200 x 3 | uint16 | LZW, predictor 2 | EPSG:32755, pixel interleaved |
| planar-u8-deflate-2band.tif | 200 x 150 x 2 | uint8 | DEFLATE | EPSG:3857, band interleaved (planar), not a COG |
| rgb-u8-deflate-pixel.tif | 230 x 170 x 3 | uint8 | DEFLATE | EPSG:3577, pixel interleaved, no predictor, not a COG |

Regenerate with:

```sh
python3 -m pip install rasterio numpy
python3 scripts/make_cog_fixtures.py
```

Design and roadmap: [wxylib/design](https://github.com/wxylib/design).

License: MIT OR Apache-2.0.

## zarr/

Zarr v3 stores made from the COG fixtures by `scripts/make_zarr_fixtures.py`
(zarr-python), holding the same pixels, so the GDAL windows in `cog/oracle`
are the oracle for Zarr readers too. `zarr/manifest.json` lists each case.

| Store | Georeferencing | Layout |
|-------|----------------|--------|
| world-int16-spatial.zarr | GeoZarr `spatial:` + `proj:` on the group, `spatial:dimensions` on the array | (y, x), chunks 100 x 100 |
| polar-f32-cf.zarr | CF-style float32 `x`/`y` centre arrays, `proj:code` only | (y, x), chunks 64 x 96, nodata -9999 |
| utm-u16-3band-node.zarr | `spatial:` + `proj:` on the array, node registration | (band, y, x), chunks 1 x 64 x 64 |

## mdimrefs/

An [mdim-refs/0.2](https://github.com/wxylib/design/blob/main/specs/mdim-refs-0.2.md)
byte-reference table over three tiny netCDF4 files shaped like NOAA OISST v2.1 daily
(time 1, zlev 1, lat 36, lon 72, int16, shuffle then deflate, one chunk per
file), made by `scripts/make_mdimrefs_fixtures.py`:

1. writes `nc/198109/oisst-avhrr-v02r01.1981090{1,2,3}.nc` (netCDF4);
2. runs `gdal mdim mosaic` over them to `oisst.vrt`;
3. scans each chunk's byte range with stock GDAL (`GetRawBlockInfo`), keyed
   within its source, copies the VRT's slabs into `placements`, and defines
   the global `refs_<array>` views from the two, in `oisst.sqlite` (the spec
   form);
4. exports each `refs_<array>` view to `refs_<array>.parquet` and the
   `arrays`, `remap` and `meta` tables to `arrays.json` (the form wxy-virtual
   reads);
5. records what GDAL reads from the mosaic for five windows in
   `oracle/*.bin` (little-endian int16, C order over time, zlev, lat, lon).

Two arrays, `sst` and `ice`. The `ice` chunk of the second day is never
written, so its row has `present = 0` and readers must give the fill value
(-999) there. Coordinates: `lat` and `lon` are affine (`{start, step}`,
latitude ascending as in OISST), `time` and `zlev` are inline.

Source paths are recorded as if mirrored at `/mirror/oisst/`; the `fixture`
remap target rewrites that prefix to `nc/`, relative to the table.

Regenerate with GDAL >= 3.12 (Python bindings and netCDF driver), netCDF4,
numpy and pyarrow, for example:

```sh
micromamba create -n wxy -c conda-forge "gdal>=3.12" libgdal-netcdf netcdf4 numpy pyarrow
python scripts/make_mdimrefs_fixtures.py
```

The output is byte-for-byte reproducible with the versions in
`mdimrefs/manifest.json`.
