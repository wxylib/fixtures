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
