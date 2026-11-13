"""Read the allele-frequency cache: the single per-site source after QC.

The cache (``alleleflux-cache-allele-freq``) holds, per MAG and timepoint, every
QC-passing sample's rows after the site-depth filter, in Parquet.  Everything
that needs per-site counts after QC (baseline presence, dN/dS, pairwise ANI,
the plotting tools) reads it through this module, so "QC-passing samples only"
and "site-filtered" hold by construction and the filter is applied exactly once.
Only the cache writer reads raw profiles.

Files are streamed in batches and filtered per batch.  The writer appends one
row group per sample, so each group spans the whole genome and Parquet
statistics cannot skip any of them by position; a whole-file read of a large
cache (measured: 1.1 GB, 180 M rows for one genome at one timepoint) does not
fit in memory, a batch does.  Label columns arrive dictionary-encoded (a
per-row-group codebook, "0 = s1_pre"); pyarrow filters them as they are and
pandas receives them as categories, which are cast to plain ``str`` on the way
out.

Example: candidate sites of MAG_A at c1:1 across pre and end ->
``load_cache_sites(cache_dir, "MAG_A", ["pre", "end"], sites)`` returns one row
per (sample, site) with the counts, plus ``timepoint`` naming the file each
row came from.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

# Columns the writer stores dictionary-encoded (pandas category); handed out
# as plain str so no consumer has to think about categories.
STR_COLUMNS = ("sample_id", "group", "subjectID", "replicate", "ref_base")
KEYS = ["contig", "position", "sample_id"]
BATCH_ROWS = 5_000_000  # ~200 MB of decoded int32/str columns per batch


def cache_path(cache_dir: str, mag_id: str, timepoint: str) -> str:
    """``{cache_dir}/{timepoint}/{mag}_{timepoint}_allele_frequency.parquet`` (the writer's layout)."""
    return os.path.join(
        cache_dir, timepoint, f"{mag_id}_{timepoint}_allele_frequency.parquet"
    )


def _open(cache_dir: str, mag_id: str, timepoint: str) -> pq.ParquetFile:
    """The cache file as a batch reader; a missing file raises, naming the path."""
    path = cache_path(cache_dir, mag_id, timepoint)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"no cache for MAG {mag_id} at timepoint {timepoint}: {path}"
        )
    return pq.ParquetFile(path)


def _site_mask(batch: pa.RecordBatch, wanted: pd.MultiIndex) -> np.ndarray:
    """True where the batch row's (contig, position) is one of ``wanted``'s."""
    rows = pd.MultiIndex.from_arrays(
        [
            batch.column("contig").to_numpy(zero_copy_only=False).astype(str),
            batch.column("position").to_numpy(),
        ]
    )
    return rows.isin(wanted)


def iter_cache_batches(
    cache_dir, mag_id, timepoint, sites=None, sample_ids=None, columns=None
):
    """Stream one cache file in filtered batches (frames with labels as ``str``).

    The building block of every cache read: ``load_cache_sites`` concatenates
    these; pairwise ANI scatters them into dense arrays as they come, so a
    1 GB file never has to fit in memory at once.  ``sites`` / ``sample_ids``
    / ``columns`` as in ``load_cache_sites``.  Yields nothing for a batch with
    no matching row.  A missing cache file raises FileNotFoundError naming it.
    """
    # The predicates, prepared once: wanted sites as a (contig, position) index,
    # wanted samples as an arrow array (ids coerced to str: numeric-looking ids
    # like 20 would otherwise never match the cache's strings).
    wanted = None
    if sites is not None:
        wanted = pd.MultiIndex.from_arrays(
            [sites["contig"].astype(str), sites["position"].astype("int64")]
        )
    keep_samples = (
        None if sample_ids is None else pa.array([str(s) for s in set(sample_ids)])
    )
    # Always read the three key columns; dedupe in case the caller named one.
    read_cols = (
        None if columns is None else KEYS + [c for c in columns if c not in KEYS]
    )
    reader = _open(cache_dir, mag_id, timepoint)
    for batch in reader.iter_batches(batch_size=BATCH_ROWS, columns=read_cols):
        mask = np.ones(batch.num_rows, dtype=bool)
        if wanted is not None:
            mask &= _site_mask(batch, wanted)
        if keep_samples is not None:
            # Works directly on the dictionary-encoded column (verified on a
            # real 230-row-group file, every batch), so no decoding step.
            mask &= pc.is_in(
                batch.column("sample_id"), value_set=keep_samples
            ).to_numpy(zero_copy_only=False)
        if not mask.any():
            continue
        # Filter in arrow BEFORE converting: only the survivors become pandas.
        df = batch.filter(pa.array(mask)).to_pandas()
        for col in STR_COLUMNS:
            if col in df.columns:
                df[col] = df[col].astype(str)  # category -> plain str
        df["timepoint"] = timepoint
        yield df


def load_cache_sites(
    cache_dir, mag_id, timepoints, sites=None, sample_ids=None, columns=None
):
    """Rows of ``mag_id``'s cache at ``timepoints``, at ``sites`` and for ``sample_ids``.

    Parameters
    ----------
    sites : DataFrame with ``contig`` and ``position`` (other columns ignored),
        or None = every position.
    sample_ids : iterable of sample ids, or None = every QC-passing sample.
    columns : extra columns to read besides ``contig, position, sample_id``, or
        None = every column.

    Returns a frame with labels as ``str`` and a ``timepoint`` column naming the
    file each row came from; empty (with the key columns) when nothing matches.
    A missing cache file raises FileNotFoundError naming it.

    Example: ``load_cache_sites(d, "MAG_A", ["pre"], sample_ids=["s2_pre"],
    columns=["A", "total_coverage"])`` -> s2_pre's rows with exactly the columns
    contig, position, sample_id, A, total_coverage, timepoint.
    """
    frames = [
        df
        for tp in timepoints
        for df in iter_cache_batches(
            cache_dir, mag_id, tp, sites=sites, sample_ids=sample_ids, columns=columns
        )
    ]
    if not frames:
        return pd.DataFrame(columns=KEYS + ["timepoint"])
    return pd.concat(frames, ignore_index=True)


def cache_samples(cache_dir, mag_id, timepoint) -> list[str]:
    """The distinct sample ids of one cache file: its QC-passing roster for that MAG.

    Reads the ``sample_id`` column only.  Example: a cache holding rows of
    s1_pre and s2_pre -> ``["s1_pre", "s2_pre"]``.
    """
    roster: set[str] = set()
    reader = _open(cache_dir, mag_id, timepoint)
    for batch in reader.iter_batches(batch_size=BATCH_ROWS, columns=["sample_id"]):
        # to_pylist() on a dictionary column yields the words, not the codes.
        roster.update(batch.column("sample_id").unique().to_pylist())
    return sorted(roster)
