# Producers checked against these tables

A producer counts for commit criterion 2 (wxylib/design
`docs/commit-criteria.md`) once its rows equal a fixture table on
presence, offset and size. `scripts/check_producer_rows.py` makes the
comparison; each run is recorded here.

| Date | Fixture | Producer | Table | Result |
|---|---|---|---|---|
| 2026-11-21 | `bran2023-temp-201001/refs_temp.parquet` (97,920 rows, two source files, built from the VirtualiZarr index in mdsumner/virtualized `d0b9c23`) | raad.virtualize `19f31e3` (store built 2026-10-28) | `s3://aadc-transfer/michae_sum/e2e/bluelink-daily-temp/flat/refs_temp.parquet` (37,530,900 rows, 1993 on) | 97,920 of 97,920 rows present in both with identical offset and size; none present on one side only |
| 2026-10-30 | `bran2023-temp-201001`, `oisst-198109` | wxylib/refharvest (async-hdf5) | OISST in its tests (`tests/fixtures.rs`); the BRAN files harvested live over HTTPS (its README) | identical rows; wxylib owns this producer, so it is a cross-check, not the criterion |

Not yet exercised by an outside producer: a chunk that was never written
(`present = 0`). The BRAN month has none. `bluelink-synthetic` and the
small `oisst.sqlite` carry that case.

To repeat the first row (read-only key for the endpoint in the
environment):

```
AWS_S3_ENDPOINT=transfer.data.aad.gov.au \
  python scripts/check_producer_rows.py \
  mdimrefs/bran2023-temp-201001/refs_temp.parquet \
  s3://aadc-transfer/michae_sum/e2e/bluelink-daily-temp/flat/refs_temp.parquet
```
