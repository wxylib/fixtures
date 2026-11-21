"""Row-level check of an outside producer's refs table against a fixture
table (commit criterion 2, wxylib/design docs/commit-criteria.md).

Compares, for every chunk of the fixture, present, offset and size with
the producer's row for the same source file and the same chunk index
within that file. The first-dimension index is taken relative to each
table's first row for the fixture's first source, so a fixture cut from
the middle of a long series lines up with the series.

    python scripts/check_producer_rows.py FIXTURE.parquet PRODUCER.parquet

PRODUCER may be a local path or an s3:// URL; for S3 the usual AWS_*
variables apply, plus AWS_S3_ENDPOINT for a custom endpoint (path style).
Sources are matched by file name, since each table spells the location
its own way (local path, /vsicurl/ or https).

Needs pyarrow only. Exit status 0 when every row agrees.
"""

import os
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.fs as fs
import pyarrow.parquet as pq


def open_table(url):
    if url.startswith("s3://"):
        kw = {}
        ep = os.environ.get("AWS_S3_ENDPOINT")
        if ep:
            kw = {"endpoint_override": ep, "force_virtual_addressing": False}
        return url[len("s3://"):], fs.S3FileSystem(**kw)
    return url, fs.LocalFileSystem()


def base(col):
    return pc.replace_substring_regex(col, r"^.*/", "")


def main(fixture, producer):
    want = pq.read_table(fixture)
    dims = [n for n in want.column_names if n.startswith("dim_")]
    want = want.append_column("file", base(want["path"]))
    files = pc.unique(pc.drop_null(want["file"])).to_pylist()

    path, filesystem = open_table(producer)
    dset = ds.dataset(path, filesystem=filesystem, format="parquet")
    cols = dims + ["present", "path", "offset", "size"]
    got = dset.to_table(
        columns=cols,
        filter=pc.match_substring_regex(ds.field("path"), "(" + "|".join(files) + ")$"),
    )
    got = got.append_column("file", base(got["path"]))
    print("fixture rows %d in %d sources; producer rows for those sources %d"
          % (want.num_rows, len(files), got.num_rows))

    # first-dimension origin: the smallest dim_0 over the shared sources
    w0 = pc.min(want["dim_0"]).as_py()
    g0 = pc.min(got["dim_0"]).as_py()
    w1 = pc.max(want["dim_0"]).as_py()
    want = want.set_column(want.schema.get_field_index("dim_0"), "dim_0",
                           pc.subtract(want["dim_0"], w0))
    got = got.set_column(got.schema.get_field_index("dim_0"), "dim_0",
                         pc.subtract(got["dim_0"], g0))
    # the fixture may stop part way through its last source
    got = got.filter(pc.less_equal(got["dim_0"], w1 - w0))

    keys = ["file"] + dims
    # absent chunks have no path: key them by index alone
    wp = want.filter(pc.equal(want["present"], 1)).select(keys + ["offset", "size"])
    gp = got.filter(pc.equal(got["present"], 1)).select(keys + ["offset", "size"])
    j = wp.join(gp, keys=keys, join_type="full outer", right_suffix="_p")
    both = pc.and_(pc.is_valid(j["offset"]), pc.is_valid(j["offset_p"]))
    same = pc.and_(pc.equal(j["offset"], j["offset_p"]), pc.equal(j["size"], j["size_p"]))
    n_both = pc.sum(both).as_py() or 0
    n_same = pc.sum(pc.and_(both, pc.fill_null(same, False))).as_py() or 0
    only_fixture = pc.sum(pc.and_(pc.is_valid(j["offset"]), pc.is_null(j["offset_p"]))).as_py() or 0
    only_producer = pc.sum(pc.and_(pc.is_null(j["offset"]), pc.is_valid(j["offset_p"]))).as_py() or 0
    n_absent_w = want.num_rows - wp.num_rows
    print("present in both: %d, identical offset and size: %d" % (n_both, n_same))
    print("present only in fixture: %d; present only in producer: %d" % (only_fixture, only_producer))
    print("absent in fixture: %d" % n_absent_w)
    ok = n_both == wp.num_rows and n_same == n_both and only_fixture == 0 and only_producer == 0
    print("AGREE" if ok else "DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
