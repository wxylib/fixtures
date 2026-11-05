# Oracle cases: format `oracle-case/0`

A fixture set asks a reference library (GDAL first) questions and keeps
its answers, so that any implementation, in any language, can be checked
against them. Questions and answers are kept apart: one set of questions,
one answer set per oracle run. GDAL 3.8 (the version in the browser
builds), 3.12 and 3.14dev can then answer the same questions side by
side, and the difference between two answer sets is a finding in itself.
The design is in wxylib/design `docs/m26-oracle-harness.md`.

```
<set>/
  cases.json                        the questions
  answers/<oracle-id>/answers.json  one oracle run's answers
  answers/<oracle-id>/*.bin, *.json array and feature answers
```

`scripts/oraclecase.py` records (`record`), compares two answer sets
(`diff`) and describes older sets in this format without touching their
files (`adopt-cog`, `adopt-vrt`). The Rust reader is the `oracle-case`
crate in wxylib/core.

## cases.json

```json
{
 "format": "oracle-case/0",
 "set": "cog",
 "description": "...",
 "cases": [
  {"id": "world-int16-deflate/w1", "op": "raster.read",
   "dsn": "world-int16-deflate.tif",
   "window": [97, 61, 130, 77], "buf": [130, 77], "bands": [1],
   "resampling": "nearest", "compare": {"kind": "exact"}}
 ]
}
```

Every case has an `id` (unique in the set; `/` is allowed), an `op` and a
`dsn`. Optional on every case: `open_options` (a list of `KEY=VALUE`),
`config` (configuration options for the duration of the case), `note`,
and `compare`.

**dsn.** Written as the implementation under test should receive it: a
path, a VSI chain, a driver-prefixed subdataset, a `vrt://` string or an
inline document. Relative names are resolved with the set directory as
the working directory, so relative paths inside `vrt://...` or
`NETCDF:"..."` resolve the same way as plain ones. `{set}` is replaced by
the set directory's absolute path, for chains such as
`/vsizip/{set}/src/tiles.zip/mosaic.vrt`.

| op | Fields | Answer |
|---|---|---|
| `raster.info` | | `info`: size, per band dtype, nodata, block, overview sizes, scale, offset; geotransform; crs (`wkt`, `code`) |
| `raster.read` | `window` [xoff, yoff, xsize, ysize] (default: all), `buf` [xsize, ysize] (default: window size), `bands` (1-based, default: all), `resampling` (`nearest`, `bilinear`, `cubic`, `cubicspline`, `lanczos`, `average`, `mode`, `gauss`, `rms`) | `values`, shape (band, row, col) |
| `mdim.info` | `array` (full name, such as `/sst`) | `info`: dims (name, size), dtype, nodata, scale, offset, crs |
| `mdim.read` | `array`, `start`, `count`, `step` (default: all) | `values`, shape in the array's dimension order |
| `vector.read` | `layer` (default: all), `bbox`, `where` | `features`: a JSON file of layers, each with fields, crs, geometry type and features (fid, properties, ISO WKT) |

**compare.** What an implementation must reach: `{"kind": "exact"}` (the
default; NaN equals NaN), `{"kind": "abs", "tol": 1e-6}` or
`{"kind": "rel", "tol": 1e-9}`, with `max_mismatch` (default 0) for a
count of cells allowed to differ. A known difference (for example GDAL's
float mode, which an implementation departs from by design) is written as
a `compare` with a `note`, not left out.

## answers.json

```json
{
 "format": "oracle-case/0",
 "set": "cog",
 "oracle": {"name": "gdal", "version": "3.8.4", "via": "osgeo.gdal",
            "release": "3.8.4", "python": "3.12.3", "numpy": "1.26.4",
            "proj": "9.4.0"},
 "recorded": "2026-11-05",
 "answers": {
  "world-int16-deflate/w1": {"status": "ok", "values": {
    "file": "answers/gdal-3.8.4/world-int16-deflate__w1.bin",
    "dtype": "int16", "shape": [1, 77, 130], "order": "C",
    "endian": "little", "sha256": "..."}},
  "v05-nodata/r0": {"status": "error", "error": "Unknown option: a_nodata"}
 }
}
```

`oracle-id` is `<name>-<version>` by default (`gdal-3.8.4`,
`gdal-3.14.0dev-739fcf03bf`). Files are named relative to the set
directory, so an adopted set can point at its existing `oracle/*.bin`.
Type names are numpy's (`uint8`, `int16`, `float32`, ...). Arrays are raw,
little-endian, C order. `status` is `ok` or `error`: a refusal is an
answer, and an implementation that accepts what GDAL refuses is reported.
`messages` lists the warnings and debug lines GDAL emitted for the case
(level, code, text).

Content is ASCII only.
