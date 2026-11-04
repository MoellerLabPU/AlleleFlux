"""Inputs for the ANI feature: the QC roster, contig lengths, dense count arrays from the cache.

Renamed from profile_io.py on 2026-10-17 when the last profile reader left; the
name now says what it reads.

The inputs, concretely
----------------------
Per-site counts come from the allele-frequency cache (``alleleflux-cache-allele-freq``):
one Parquet file per MAG and timepoint holding every QC-passing sample's rows
after the site-depth filter.  One cache row, and where its numbers land::

    contig        position  sample_id  A  C  G  T
    k141_102792   0         S1         2  0  1  1

    dense["k141_102792"][0] == [2, 0, 1, 1]      (uint16, order A,C,G,T)

Two facts drive everything here.  The cache holds **only covered positions** --
an absent row means zero coverage, not missing data -- so the dense arrays start
at zero and rows are scattered into them.  And positions are **0-based**, so a
position is a direct array index.  Base columns stay in the cache's own
``A, C, G, T`` order: the set-based classifier never depends on column order.

Raw profile files are not read here (nor anywhere after QC except the cache
writer): every ANI number rests on the same rows the statistics saw.

Library module: imported by the rest of the ani package and by tests; never run
directly, hence no ``if __name__ == "__main__"`` block.  The one runnable entry
point of this feature is the ``alleleflux-pairwise-ani`` command (pairwise_ani.py).
"""

import logging

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.ani.classify import BASES
from alleleflux.scripts.utilities.cache_io import cache_samples, iter_cache_batches
from alleleflux.scripts.utilities.utilities import build_contig_length_index

logger = logging.getLogger(__name__)


# Metadata carried through to the output tables.  Everything is read as str:
# group labels can be numbers (20/40) and subject ids plain integers, and letting
# pandas sniff them as int64 makes every later string comparison match nothing.
_QC_META_COLUMNS = ["sample_id", "subjectID", "group", "time", "replicate"]
# genome_size rides along (numeric, NOT str) so the CLI can cross-check its own
# reference-derived genome length against the value QC recorded -- catching "QC
# ran against a different FASTA/mapping" loudly instead of silently.  breadth is
# copied as recorded into the per-sample table (the cache is site-filtered, so
# QC's unfiltered breadth cannot be re-derived from it).  average_coverage is
# NOT needed here: the cache this CLI reads was already filtered against it.
_QC_KEEP_COLUMNS = _QC_META_COLUMNS + ["genome_size", "breadth"]
_QC_DTYPES = {column: str for column in _QC_META_COLUMNS}

# Dense arrays use uint16 (2 bytes/cell): memory scales as genome_length x 4 x
# n_samples, ~2.2 GB for a 3 Mb MAG across 90 samples.  Counts above 65,535 are
# clipped -- presence thresholds are two digits at most, so a clip can never flip
# a verdict; it only nudges a frequency by <0.1% at absurd depth.
_DENSE_DTYPE = np.uint16
_DENSE_MAX = np.iinfo(_DENSE_DTYPE).max


def qc_passing_samples(qc_files: list[str], mag_id: str) -> pd.DataFrame:
    """Samples that cleared AlleleFlux QC for this MAG, unioned across QC files.

    One QC file exists per MAG **per timepoint combination** and lists EVERY
    sample of that combination (passing or not).  Combinations overlap on shared
    timepoints, so the same sample can be judged in several files; the union keeps
    a sample wherever it passed at least once, collapsed to one row per
    ``sample_id``.

    Worked example (combinations pre_end and pre_post share "pre")::

        QC/QC_pre_end/MAG_X_QC.tsv    60 rows, 30 pass   (pre and end samples)
        QC/QC_pre_post/MAG_X_QC.tsv   60 rows, 29 pass   (pre and post samples)

        qc_passing_samples([qc_pre_end, qc_pre_post], "MAG_X")
        -> 43 rows: 59 passing rows minus 16 "pre" samples that passed in both
           files (their duplicate rows are identical, so first-wins loses nothing)

    Note the caller chooses the scope: pass ONE file and the roster is exactly
    that combination's samples; pass all files (what the Snakemake rule does) and
    the roster spans every timepoint, which is what within-subject cross-timepoint
    pairs need.

    Parameters
    ----------
    qc_files : list of str
        Per-MAG QC TSVs (``{mag}_QC.tsv``), one per timepoint combination --
        quality_control.py's real output files.
    mag_id : str
        Only used for log/error context; the files are already per-MAG.

    Returns
    -------
    pandas.DataFrame
        One row per QC-passing sample with the metadata columns present in the
        files (``sample_id, subjectID, group, time, replicate`` as str, plus
        ``genome_size`` numeric for the CLI's reference cross-check).
        MAY BE EMPTY: a low-abundance MAG where no sample passes QC is a
        legitimate outcome -- the caller writes header-only outputs, not an error.

    Raises
    ------
    ValueError
        If a QC file lacks ``coverage_threshold_passed`` (malformed, or
        quality_control.py did not finish), or if no QC files are supplied at all.
    """
    frames = []
    for qc_file in qc_files:
        # Force the label columns to str on the way in (_QC_DTYPES note above);
        # columns not named there, like the bool verdict, keep pandas' inference.
        table = pd.read_csv(qc_file, sep="\t", dtype=_QC_DTYPES)
        if "coverage_threshold_passed" not in table.columns:
            raise ValueError(
                f"QC file {qc_file} has no 'coverage_threshold_passed' column; "
                "was quality_control.py run to completion?"
            )
        # The verdict column arrives as real bools (pandas parses TSV True/False);
        # keep only rows that passed.
        passing = table[table["coverage_threshold_passed"] == True]  # noqa: E712
        # Keep whatever metadata columns this file has -- 'time' is absent for
        # single-timepoint data, and the outputs simply omit what was never there.
        available = [c for c in _QC_KEEP_COLUMNS if c in passing.columns]
        frames.append(passing[available].copy())

    if not frames:
        raise ValueError(f"No QC files supplied for MAG {mag_id}")

    # A sample judged usable at several timepoint combinations appears in several
    # files; keep its first row. Safe because group/subject/time come from the one
    # input metadata sheet, so duplicate rows carry identical values.
    combined = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates(subset=["sample_id"])
        .reset_index(drop=True)
    )
    if combined.empty:
        logger.warning(
            f"MAG {mag_id}: no sample passed QC in any of {len(qc_files)} QC file(s)"
        )
    else:
        logger.info(f"MAG {mag_id}: {len(combined)} samples passed QC")
    return combined


def contig_lengths_for_mag(
    fasta_path: str, mag_mapping_path: str, mag_id: str
) -> dict[str, int]:
    """Contig -> length for the contigs belonging to one MAG.

    A thin per-MAG wrapper around the shared
    ``utilities.build_contig_length_index``, which prefers the samtools ``.fai``
    index beside the FASTA (milliseconds) over parsing the 0.5 GB reference
    (seconds, repeated by every per-MAG job).  What this wrapper adds is the
    ANI-specific contract: subset to ONE MAG's contigs, and fail loudly when the
    mapping and the reference disagree -- dense arrays sized from wrong lengths
    would silently corrupt every comparison downstream.

    Parameters
    ----------
    fasta_path : str
        The reference FASTA (config ``input.fasta_path``); ``{fasta_path}.fai``
        is used by the shared helper when present.
    mag_mapping_path : str
        TSV with columns ``mag_id, contig_id`` (config ``input.mag_mapping_path``).
    mag_id : str
        Which MAG's contigs to return.

    Raises
    ------
    ValueError
        If the mapping lacks its two required columns, names no contigs for this
        MAG, or names contigs the reference does not contain (mapping and FASTA
        out of sync).
    """
    # dtype=str: MAG and contig ids are labels; a numeric-looking id must not
    # silently become int64 and then fail every string comparison.
    mapping = pd.read_csv(mag_mapping_path, sep="\t", dtype=str)
    missing_columns = {"mag_id", "contig_id"} - set(mapping.columns)
    if missing_columns:
        raise ValueError(
            f"MAG mapping {mag_mapping_path} is missing columns "
            f"{sorted(missing_columns)}; found {list(mapping.columns)}"
        )

    # The mapping covers the whole assembly (all MAGs); keep this MAG's contigs.
    wanted = set(mapping.loc[mapping["mag_id"] == mag_id, "contig_id"])
    if not wanted:
        raise ValueError(f"MAG {mag_id} has no contigs listed in {mag_mapping_path}")

    # The shared helper does the actual length reading (.fai fast path when
    # present).  It returns lengths for every mapped contig genome-wide; that is
    # cheap (the index is ~1.5 MB) and keeps ONE implementation serving QC,
    # positions-QC and ANI alike.
    all_lengths = build_contig_length_index(fasta_path, mag_mapping_path)

    # Iterate the REFERENCE-derived dict (built in .fai order), not `wanted`:
    # `wanted` is a set of strings, and Python randomises string hashing per
    # process, so iterating it would give a different contig order on every run
    # -- and the engine's contig-outer loop would emit SNP-location rows in a
    # different order each time. This keeps the output reproducible.
    lengths = {
        contig: int(length)
        for contig, length in all_lengths.items()
        if contig in wanted
    }

    # Every mapped contig must exist in the reference; anything left over means
    # the mapping and the FASTA are from different assemblies.
    absent = wanted - set(lengths)
    if absent:
        raise ValueError(
            f"MAG {mag_id}: {len(absent)} contig(s) from the mapping are not in "
            f"{fasta_path}, e.g. {sorted(absent)[:3]}. Mapping and reference disagree."
        )
    return lengths


def scatter_rows(
    dense: dict, rows: pd.DataFrame, contig_lengths: dict[str, int]
) -> None:
    """Write ``rows`` (contig, position, A, C, G, T) into the dense arrays, in place.

    ``dense`` is ``{contig -> (length, 4) uint16 zeros}`` for EVERY contig of the
    MAG (``dense_counts_from_cache`` allocates it); positions absent from ``rows``
    stay zero, which is exactly what absence means.  Called once per cache batch,
    so a sample's rows may arrive in pieces.  A contig outside ``contig_lengths``
    (cache and MAG mapping out of sync) or a position beyond its contig's length
    (cache and reference out of sync) raises: contract violations, not data to drop.

    Example: rows (c1, 0, A=4) and (c1, 3, C=6) into a 5-long c1 array ->
    ``dense["c1"][0] == [4, 0, 0, 0]``, ``dense["c1"][3] == [0, 6, 0, 0]``, the
    other three positions all zero.
    """

    # sort=False: iterate contigs in file order, no pointless re-sort.
    # observed=True: if a caller ever hands us categorical contigs, don't explode
    # into every unobserved category (the regional_contrast bug of 2026-07).
    for contig, sub in rows.groupby("contig", sort=False, observed=True):
        if contig not in dense:
            raise ValueError(
                f"Cache contains contig {contig!r} that is not in this MAG's "
                f"mapping ({len(contig_lengths)} contigs). Cache and mapping disagree."
            )
        length = contig_lengths[contig]
        # 0-based positions ARE array indices; bounds-check before scattering
        # (initial= keeps max/min safe even on an empty slice).
        positions = sub["position"].to_numpy()
        if positions.max(initial=-1) >= length or positions.min(initial=0) < 0:
            raise ValueError(
                f"Cache positions for contig {contig} fall outside its length "
                f"{length}: min={positions.min()}, max={positions.max()}. "
                "Cache and reference FASTA disagree."
            )
        counts = sub[list(BASES)].to_numpy()
        # Saturate instead of overflowing on the way down to uint16 (see the
        # _DENSE_DTYPE note at the top of the module).
        np.clip(counts, 0, _DENSE_MAX, out=counts)
        # Fancy-index scatter: one vectorised write per contig instead of a Python
        # loop over ~3M rows.
        dense[contig][positions] = counts.astype(_DENSE_DTYPE)


def dense_counts_from_cache(
    cache_dir: str,
    mag_id: str,
    timepoint: str,
    sample_ids: list[str],
    contig_lengths: dict[str, int],
    min_cov: int,
) -> dict[str, tuple[dict[str, np.ndarray], dict]]:
    """Dense ``(contig_length, 4)`` arrays per sample, from one timepoint's cache file.

    The cache is the single per-site source after QC (QC-passing samples only,
    site-filtered).  The file is streamed in batches and each batch's rows are
    scattered into the owning sample's arrays as they come, so memory is the
    dense arrays themselves (genome_length x 4 x 2 bytes per sample) plus one
    batch, never the whole file.  Every requested sample must have rows in the
    cache; a requested sample with none means the cache was built from a
    different QC table than the roster came from, and that raises.

    Returns ``{sample_id: (dense_by_contig, stats)}`` with ``stats`` =
    ``{"n_positions": rows in the cache (covered after the site filter),
    "breadth_minCov": fraction of the genome with >= min_cov reads}``.

    Example: a 6 bp contig, S2 covered 20x at five positions and 2x at one,
    min_cov 5 -> n_positions 6, breadth_minCov 5/6.
    """
    wanted = [str(s) for s in sample_ids]
    genome_length = int(sum(contig_lengths.values()))
    dense = {
        sid: {
            c: np.zeros((n, len(BASES)), dtype=_DENSE_DTYPE)
            for c, n in contig_lengths.items()
        }
        for sid in wanted
    }
    n_positions = {sid: 0 for sid in wanted}
    n_deep = {sid: 0 for sid in wanted}
    for batch in iter_cache_batches(
        cache_dir, mag_id, timepoint, sample_ids=wanted, columns=list(BASES)
    ):
        for sid, rows in batch.groupby("sample_id", sort=False):
            scatter_rows(dense[sid], rows, contig_lengths)
            n_positions[sid] += len(rows)
            n_deep[sid] += int((rows[list(BASES)].sum(axis=1) >= min_cov).sum())
    missing = [sid for sid in wanted if n_positions[sid] == 0]
    if missing:
        raise ValueError(
            f"MAG {mag_id} at {timepoint}: QC-passing sample(s) {missing} have no rows in "
            f"the cache; the cache was built from a different QC table than --qc_files"
        )
    return {
        sid: (
            dense[sid],
            {
                "n_positions": n_positions[sid],
                "breadth_minCov": n_deep[sid] / genome_length,
            },
        )
        for sid in wanted
    }
