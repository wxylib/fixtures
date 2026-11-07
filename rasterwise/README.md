# rasterwise/

Questions, GDAL answers and Assertion documents over the published files of
the rgdal-dev/rasterwise problem-file archive. The files are read in place,
over `/vsicurl/`, from the archive's `blobs-1` release assets (named by
sha256); nothing is copied here except small answers. Design: wxylib/design
`docs/m30-legacy-corpus.md`.

| Path | What |
|---|---|
| `files.json` | the 58 published data files (netCDF, HDF4, TIFF, raw binary, VRT) at the rasterwise commit named inside, by sha256, with url and case |
| `cases.json` | `oracle-case/0` questions: open, info and a small read per classic subdataset, info per multidimensional array, and a read of every 1D or 2D array up to 262144 elements (coordinates are the evidence for assertions) |
| `answers/gdal-3.8.4/` | GDAL 3.8.4's answers |
| `assertions/` | Assertion documents (`<sha12>.json`) and the VRTs they emit |

Regenerate:

```sh
python3 scripts/rasterwise_cases.py path/to/rasterwise   # cases.json, files.json
python3 scripts/oraclecase.py record rasterwise          # answers/gdal-<version>/
python3 scripts/rasterwise_assert.py --check             # assertions/, with checks
```

Set `CURL_CA_BUNDLE` behind a proxy, and `GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR`
to stop GDAL listing the release directory.

## What GDAL 3.8.4 says (rasterwise 091f2fc)

861 questions over 58 files; 832 answered, 29 refused.

| Files | Open | Notes |
|---|---|---|
| netCDF3 (21), netCDF4 (16) | 37 of 37 | the WOD single-cast profile opens only as multidim (not a raster); 9 char-variable subdatasets (Argo flags, OSCAR `date`) are listed by the classic view but cannot be read |
| GeoTIFF and TIFF (3) | 3 of 3 | |
| VRT (10) | 5 of 10, 0 read | the hand-written UCSD topography VRTs name their sources by relative path (`../topo_18.1.img`), which content-addressed storage does not keep |
| raw binary (7) | 0 of 7 | headerless by design: an assertion is the only way in |
| HDF4 (1) | 0 of 1 | the HDF4 library needs a local file |

Opening is not reading correctly; the cases say which views are wrong.

## Assertion documents

An Assertion supplies what a file lacks, or corrects what it gets wrong,
without copying a byte. Each document (`rasterwise-assertion/0`) is keyed
by sha256, points at the published url, restates the case's verified
recipe as a `vrt://` string or a VRT beside it, and records a check run by
`rasterwise_assert.py --check`:

| File | Case | Emits | Check (GDAL 3.8.4) |
|---|---|---|---|
| 13cc4140fb59 | cryosat2-sea-ice-cryptic-curvilinear | `vrt://`, EPSG:3412, south-up geotransform | asserted cell centres vs the file's 2D lon/lat: max 43 m, RMS 14 m |
| 549a1b01dda6 | nsidc-polar-grid-no-coords-flipped | `vrt://`, EPSG:3412 at 12.5 km | view equals stored rows flipped |
| f43715f25daf | nhsce-polar-stereo-labelled-longlat | `vrt://`, polar stereographic | vs 2D lon/lat: max 7.5 m |
| df20ee33704c | gpm-dpr-swapped-axes | mdim VRT: transpose, reverse lat, CRS, coordinates from the file | equals stored slice transposed; classic view origin (-180, 67), 0.25 degree |
| 710a226e331a | aad-fast-ice-raw-cea-grid-as-lonlat-arrays | raw VRT (flip by negative line offset), CEA | vs lats.img/lons.img through the same raw layout: max 0.56 m |
| 0804ca23ec39 | ucsd-topo-raw-mercator-headers-disagree | raw VRT (0..360) + halves-swapped VRT (-180..180), EPSG:3857 | warped vs topo_19.1.nc, 140-145 E 35-40 N: RMS 4.2 m |
| 477d7e804c65 | ucsd-topo-raw-mercator-headers-disagree | same, topo_27.1 | vs topo_27.1.nc: RMS 3.8 m |

Two things the assertions show:

- **No configuration option is needed for the netCDF flips.** The recipes
  use `--config GDAL_NETCDF_BOTTOMUP NO`, which a VRT cannot remember. GDAL's
  default view of these files is south-up, and a south-up geotransform
  (origin at the south edge, positive y step) places that view correctly,
  so the whole assertion fits in one `vrt://` string.
- **GDAL 3.8.4 cannot warp the 0..360 Mercator grid.** The warper fails to
  invert a 0 to 40075016 m EPSG:3857 extent ("Cannot find coordinate
  operations") and writes zeros. The halves-swapped -180..180 view warps.
  Raw sources also need `GDAL_VRT_RAWRASTERBAND_ALLOWED_SOURCE=ALL` on GDAL
  3.10 and later (3.8.4 does not check); each document carries it under
  `emit.config`.
