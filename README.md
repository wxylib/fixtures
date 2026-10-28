# wxylib/fixtures

Test data and GDAL oracle outputs for wxylib, kept out of
[wxylib/core](https://github.com/wxylib/core) so it stays small.

Every case is a small input plus what GDAL gives for it; the wxylib/core
tests must reproduce GDAL's numbers exactly. Each directory has a script
in `scripts/` that regenerates it and a `manifest.json` listing its cases.

The warp-like oracles (warp/, curvilinear/, cf/, index/) track GDAL HEAD
(3.14dev), whose per-chunk kernel scale is the one rwarp ports. The
`oracles` workflow regenerates them in hypertidy/gdal-r-ci's
`gdal-r-python:dev` image and pushes the result to a new
`oracles/gdal-<version>-<run>` branch: run it from the Actions tab and
name the sets, for example `curvilinear cf index`.

| Directory | Inputs | Oracle | Checks |
|-----------|--------|--------|--------|
| [cog/](#cog) | five small GeoTIFFs (COG and plain) | GDAL window reads | `wxy-tiff`, `wxy-virtual` |
| [zarr/](#zarr) | three of the COGs as Zarr v3 | the same windows | `wxy-zarr` |
| [warp/](#warp) | target grids over the COGs | GDAL warps (near, bilinear, cubic) | `wxy-warp`, `wxy-cli` |
| [mdimrefs/](#mdimrefs) | byte-reference tables over netCDF4 (OISST) | GDAL reads of the mosaic | `wxy-virtual` |
| [curvilinear/](#curvilinear) | swaths placed by 2D lon/lat | GDAL geolocation-array warps | `wxy-warp` |
| [cf/](#cf) | an OISST-style store as xarray writes it | GDAL warps per day | `wxy-warp`, `wxy-cli` |
| [points/](#points) | random points over each COG | GDAL point reads, by cell and interpolated | `wxy-warp` |
| [polygons/](#polygons) | random polygons with holes | the cells GDAL's rasterizer burns | `wxy-warp` |
| [index/](#index) | index space against coordinate space: a ROMS-style store, a rotated-pole store | ROMS depths, a GDAL geolocation warp, PROJ ob_tran | `wxy-zarr`, `wxy-cli` |

## cog/

Small Cloud Optimized GeoTIFFs with synthetic pixel values, made by
`scripts/make_cog_fixtures.py`, and what GDAL (through osgeo.gdal) reads for a
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
micromamba create -n gdal -c conda-forge "gdal>=3.12" numpy zarr xarray pandas
micromamba run -n gdal python scripts/make_cog_fixtures.py
```

## warp/

Warp oracles, made by `scripts/make_warp_fixtures.py` after the COGs: what
GDAL's warper (through `gdal.Warp` in osgeo.gdal, approximate transformer
at 0.125 pixels) writes for COG sources onto target grids, with nearest,
bilinear and cubic. Some targets are near 1:1 with the source; others are
coarser, where GDAL scales its kernels (stage S5), read an overview, or
need the source window split into chunks: `polar-to-4326-global` splits
into 16 where most of the target misses the source, and
`world-to-3031-wm1` into 4 under a 1 MB warp memory limit (gdalwarp
`-wm 1`). `manifest.json` gives each case's source, target CRS,
geotransform, size, overview and `warp_memory_mb`; each result is
a raw little-endian `.bin` file of the target's cell type, C order (band, row, col).

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

`oisst.sqlite` is not read by the wxylib tests; it is kept because it is the
spec form of the table (what zaro and GDAL's OGR driver read), and the
Parquet export is made from it.

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

### mdimrefs/oisst-198109/

The same table for 30 days of real NOAA OISST v2.1 (September 1981, four
arrays), scanned over HTTP from the public NOAA bucket by
`scripts/make_mdimrefs_oisst_198109.py`. Sources are the `/vsicurl/https://`
URIs GDAL opened; the `public` remap target makes them plain `https://`.
`manifest.json` holds the SHA-256 of what GDAL reads for day 30 of each
array, for the wxylib/core demo `cargo run -p wxy-virtual --example
oisst_refs`. Not used in CI; regenerating needs network access.

### mdimrefs/bran2023-temp-201001/

32 days of BRAN2023 `temp` (CSIRO Bluelink reanalysis, 2010-01-01 to
2010-02-01, 51 depths, 1500 x 3600), the many-chunks-per-file case: 300 x
300 int16 chunks, 3060 per file per day, 97,920 rows over two netCDF files
on NCI's THREDDS server. Converted by
`scripts/make_mdimrefs_bran2023_temp.py` from the kerchunk Parquet index in
[mdsumner/virtualized](https://github.com/mdsumner/virtualized) (pinned
commit; built on NCI with R, rhdf5 and GDAL multidim), not scanned here.
Attributes are as the index has them, including its rounded
`scale_factor` (0.0078; the netCDF says 0.007782219909131527).
`manifest.json` names two oracle windows (one chunk, and one row across 12
chunks); `--oracle` adds the SHA-256 of what GDAL reads for each from the
source netCDF, which needs network access to thredds.nci.org.au. For the
wxylib/core demo `cargo run -p wxy-virtual --example bran_refs`. Not used
in CI.

## curvilinear/

Curvilinear warp oracles (wxylib stage S5b), made by
`scripts/make_curvilinear_fixtures.py`: float32 grids placed by 2D
longitude and latitude arrays at cell centres (a curved swath off
south-eastern Australia, and a ring around the south pole crossing the
antimeridian), and what GDAL's warper writes for them through its
geolocation-array transformer (`GEOREFERENCING_CONVENTION=PIXEL_CENTER`,
approximate transformer at 0.125 pixels), with nearest, bilinear and cubic,
to EPSG:4326 (also coarse) and EPSG:3031. Each grid's arrays are raw
little-endian files (`*-lon.bin`, `*-lat.bin` float64, `*-values.bin`
float32) next to the GeoTIFFs GDAL read; `manifest.json` lists the cases.
Cases with `"cell": "edge"` place the swath by its cell corners instead:
`*-edges-lon.bin` and `*-edges-lat.bin` are one row and column larger
than the grid, read by GDAL with `GEOREFERENCING_CONVENTION=TOP_LEFT_CORNER`.

## cf/

CF Zarr oracles (wxylib stage S5b), made by `scripts/make_cf_fixtures.py`
(needs xarray and zarr 3): `oisst-like.zarr` is laid out as xarray writes
an OISST-style month (time x zlev x lat x lon, packed int16 with
`scale_factor` and `_FillValue`, float32 1D coordinates, latitude
ascending, longitude 0 to 360, chunked two days deep). The `.bin` files
hold what GDAL's warper writes for every day, from a north-up GeoTIFF of
that day, to EPSG:3031 with nearest and bilinear (int16, days stacked).

## points/

Point-read oracles (wxylib stage S4c), made by
`scripts/make_point_fixtures.py`: 400 random points over each COG in
`cog/`, given in EPSG:4326 and in the file's own CRS (some fall off the
file), and the values GDAL reads there as gdallocationinfo does
(reproject, inverse geotransform, floor; every band). Per case:
`*-x.bin`, `*-y.bin` (float64), `*-values.bin` (float64, bands x n) and
`*-found.bin` (uint8, 1 on the file); `manifest.json` lists the cases.
Interpolated reads are GDAL's own InterpolateAtGeolocation, called
through osgeo.gdal: `*-bilinear.bin`,
`*-cubic.bin`, `*-cubicspline.bin` (float64, bands x n, NaN where GDAL has
no value), with GDAL's (pixel, line) for each point in `*-pixel.bin` and
`*-line.bin`.

## polygons/

Polygon-to-cell oracles (wxylib stage S4c), made by
`scripts/make_polygon_fixtures.py`: random star-shaped polygons (some
with holes, both windings) over each COG in `cog/`, in EPSG:4326 and in
the file's own CRS (plus, in the native cases, a rectangle on cell edges
and one through cell centres), and the cells GDAL's rasterizer burns for
each polygon alone, by cell centre and with ALL_TOUCHED. Per case:
`*-coords.bin` (float64 x, y pairs), `*-rings.json` (ring point counts
per polygon) and `*-cells.json` (flat cell indices as [start, length]
runs, per mode and polygon).

Documentation: **https://wxylib.org/docs/**, one place for what works, how-tos,
explainers, the design documents and how to help.

License: MIT OR Apache-2.0.

## index/

Index space against coordinate space (design m07, section 5), made by
`scripts/make_index_fixtures.py`: small Zarr v3 stores as a real
producer writes them, plus the truth that producer threw away.

- **I1 `roms-fan.zarr`**: a ROMS-style coastal grid, regular in index
  space and curvilinear on the globe (playproj's `fan()`), with ROMS's own
  names and CF attributes: `temp(s_rho, eta_rho, xi_rho)` (8 x 46 x 24)
  placed by 2D `lon_rho` and `lat_rho`, and a terrain-following vertical
  (`s_rho` as CF `ocean_s_coordinate_g2` with `formula_terms`, `Cs_r`,
  `h`, `zeta`, `hc`, Vtransform 2, Vstretching 4). Oracles: every rho and
  w depth (`roms-fan-z-rho.bin`, `roms-fan-z-w.bin`, float64) from the
  NumPy port of ROMS set_scoord.F and set_depth.F in playproj's
  `test/make-roms-fixture.py`, and GDAL's geolocation warp of the surface
  layer to EPSG:4326, near and bilinear (`roms-fan-top.vrt` shows how GDAL
  was given it).
- **I3 `rotated.zarr`**: a 1 degree EURO-CORDEX window on a CF
  `rotated_latitude_longitude` grid: `tas(rlat, rlon)`, 1D `rlon`, `rlat`,
  the `rotated_pole` grid mapping and float32 2D `lon`, `lat`. The
  manifest gives the generating affine and the PROJ `ob_tran` string;
  osgeo.osr agrees with the stored 2D arrays to the float32 rounding
  recorded there.

Vertical positions are recorded but not yet read by wxylib: a reader of
`formula_terms` waits on a numbered design document.
