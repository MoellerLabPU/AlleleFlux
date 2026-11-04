"""Pure rules of the de novo analysis, shared by baseline presence and rising origin.

No I/O here: every function takes and returns DataFrames or numbers, so both
CLIs and the tests call exactly the same code.  Spec:
``docs/superpowers/specs/2026-10-14-site-depth-filter-and-de-novo-phases.md``
section 2.

* :func:`assign_blocks` -- the cluster rule (spec 2.3): de novo alleles that sit
  within a window of each other on one contig, in one group, are one block
  event, never counted as independent de novo mutations.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# The two de novo flavours of rising origin's ``origin`` vocabulary: a covered,
# clean baseline (``de_novo``) or a baseline with a few reads under the
# detection bar (``de_novo_below_detection``).  Both are candidate new alleles,
# so both take part in clustering, exactly as both are relabelled by the strain
# gate; standing alleles and the two non-verdicts never do.
DE_NOVO_ORIGINS = ("de_novo", "de_novo_below_detection")
BLOCK_EVENT = "block_event"


def assign_blocks(
    alleles: pd.DataFrame, window_bp: int, max_neighbours: int
) -> pd.DataFrame:
    """Group de novo alleles that sit near each other into block events (spec 2.3).

    Two de novo alleles are *neighbours* when they are on the same contig of the
    same MAG, rose in the same group (``rose_in``) and are at most ``window_bp``
    apart (inclusive); which comparison or age they rose in does not matter, so
    an imported chunk seen across two comparisons is still one chunk.  Blocks
    are the connected components of that relation (single linkage: 0, 800 and
    1600 are one block with a 1000 bp window although the ends are 1600 apart).

    An allele with more than ``max_neighbours`` neighbours is relabelled
    ``block_event`` and gets ``block_id = "{mag}:{contig}:{group}:{k}"`` (k =
    1-based block number in position order within that MAG x contig x group);
    an allele with ``max_neighbours`` or fewer keeps its origin and an empty
    block id.  Rows whose origin is not de novo are untouched and never count
    as neighbours.  The same (MAG, contig, position, allele) appearing in
    several comparisons is ONE allele: it is never its own neighbour, and all
    its rows get the same labels.

    Two different letters at the same position (only possible at a site where
    three or four bases were tested and two of them each cleared the
    significance line, e.g. baseline A, later some mice G and others T) are
    two alleles 0 bp apart, so they are neighbours and form a block.  Two
    independent mutations hitting the same base in the same group of mice is
    about as unlikely as it gets; what produces this picture is reads from a
    related strain mapping onto this genome with a mismatch at that base, or an
    imported chunk carrying a different letter there.  Neither is two clean de
    novo events, so the block label is the conservative reading (the 2026-10-07
    audit traced such same-site clumps to cross-mapped relatives).  At a
    two-letter site this cannot arise: one letter's rise is the other's fall,
    so only the riser gets a rising row.

    Parameters
    ----------
    alleles
        Rows with ``mag_id, contig, position, allele, rose_in, origin``
        (rising origin's allele table after ``load_rising``).
    window_bp, max_neighbours
        The two config keys ``cluster_window_bp`` and ``max_de_novo_neighbours``.

    Returns
    -------
    A copy of ``alleles`` in the same row order with ``block_id`` (str, ``""``
    when not clustered) and ``n_de_novo_neighbours`` (int) added and ``origin``
    relabelled where clustered.

    Example (window 1000, max 0): de novo alleles at 100, 900 and 1600 on one
    contig in ``fat`` -> one block ``M:c1:fat:1``, neighbours 1, 2, 1, all three
    ``block_event``; an allele at 5000 on the same contig stays ``de_novo`` with
    0 neighbours.  With max 2 all four stay de novo (900 has exactly 2).
    """
    out = alleles.copy()
    out["block_id"] = ""
    out["n_de_novo_neighbours"] = 0
    if out.empty:
        return out
    is_dn = out["origin"].isin(DE_NOVO_ORIGINS)
    groups = out[is_dn].groupby(["mag_id", "contig", "rose_in"], sort=True).groups
    for (mag, contig, group), idx in groups.items():
        sub = out.loc[idx]
        # One node per distinct (position, allele): rows of the same allele in
        # two comparisons share the node and are not each other's neighbour.
        nodes = (
            sub[["position", "allele"]]
            .drop_duplicates()
            .sort_values(["position", "allele"])
            .reset_index(drop=True)
        )
        pos = nodes["position"].to_numpy(dtype=np.int64)
        # Neighbours of node i: every OTHER node within the window (inclusive).
        # n nodes per MAG x contig x group is small (tens to hundreds), so the
        # n x n distance matrix is cheap.
        dist = np.abs(pos[:, None] - pos[None, :])
        adjacent = (dist <= window_bp) & ~np.eye(len(pos), dtype=bool)
        n_neighbours = adjacent.sum(axis=1)
        # Connected components by single linkage along sorted positions: a new
        # component starts wherever the gap to the previous node exceeds the
        # window (any two nodes closer than that are linked, directly or via
        # the nodes between them).
        component = np.concatenate([[0], np.cumsum(np.diff(pos) > window_bp)])
        clustered = n_neighbours > max_neighbours
        # Dense 1-based numbering over the components that actually hold a
        # clustered node, in position order.
        block_number = {
            c: k + 1 for k, c in enumerate(sorted(set(component[clustered])))
        }
        key_of_node = list(zip(nodes["position"], nodes["allele"]))
        neighbours_of = dict(zip(key_of_node, n_neighbours))
        block_of = {
            k: f"{mag}:{contig}:{group}:{block_number[c]}"
            for k, c, is_c in zip(key_of_node, component, clustered)
            if is_c
        }
        row_keys = list(zip(sub["position"], sub["allele"]))
        out.loc[idx, "n_de_novo_neighbours"] = [int(neighbours_of[k]) for k in row_keys]
        out.loc[idx, "block_id"] = [block_of.get(k, "") for k in row_keys]
        hit = [k in block_of for k in row_keys]
        out.loc[idx[hit], "origin"] = BLOCK_EVENT
    return out


# ---------------------------------------------------------------------------
# Carrier-level strain gate (spec 2.4)
# ---------------------------------------------------------------------------

CARRIER_NEW_STRAIN = "new_strain_carriers"
CARRIER_SAME_STRAIN = "same_strain_carriers"
CARRIER_UNDETERMINED = "strain_undetermined"

# How a per-mouse verdict reads after a round trip through a TSV read with
# dtype=str ("True" / "False" / blank -> NaN) or straight from the turnover
# step's nullable-boolean column (True / False / pd.NA).
_VERDICT_WORDS = {"True": True, "False": False, True: True, False: False}


def _verdict(value):
    """A turnover verdict as True / False / None (undetermined); other text raises."""
    if (
        value is None
        or (isinstance(value, float) and np.isnan(value))
        or value is pd.NA
    ):
        return None
    if value in _VERDICT_WORDS:
        return _VERDICT_WORDS[value]
    raise ValueError(
        f"dominant_strain_change must be True, False or blank, got {value!r}"
    )


def carrier_strain_verdicts(
    carriers, turnover: pd.DataFrame, period: str, min_determined: int
) -> tuple[str, int, int]:
    """Did the later-timepoint carriers of an allele keep their strain? (spec 2.4)

    Reuses the per-mouse verdicts the strain-turnover step already computed
    (``{mag}_strain_turnover.tsv``: one row per mouse and transition, the pair
    oriented by time, ``dominant_strain_change`` = conANI below the threshold,
    blank when too little of the genome was compared), so the thresholds and
    the pair selection live in that step alone and are never re-derived here.

    For each carrier mouse take its row for ``period`` (``transition``, spelled
    ``{earlier}_{later}``).  A mouse with TWO rows for that transition (sampled
    twice at one timepoint) is a data error and raises, naming the mouse: the
    designs this serves sample a mouse once per timepoint, and silently picking
    one of its rows would hide the problem.  A carrier is *determined* when its verdict is
    True or False (not blank, not missing); *changed* when it is True.  A
    carrier listed twice raises (the list is written from a set).  Fewer
    determined carriers than ``min_determined`` -> ``strain_undetermined``;
    more than half of the determined ones changed -> ``new_strain_carriers``;
    otherwise ``same_strain_carriers``.

    Returns ``(status, n_determined, n_changed)``.

    Example: carriers m1, m2, m3 with verdicts True, True, False at floor 3 ->
    ``("new_strain_carriers", 3, 2)``; with m3 blank at floor 3 ->
    ``("strain_undetermined", 2, 2)``; two of four changed ->
    ``("same_strain_carriers", 4, 2)`` (half is not more than half).
    """
    # One row per mouse for this transition; two rows for one mouse means the
    # mouse was sampled twice at a timepoint, which these designs never do.
    rows = turnover[turnover["transition"] == period]
    dup = rows["subjectID"].astype(str).duplicated(keep=False)
    if dup.any():
        raise ValueError(
            f"{period}: more than one turnover row for mouse/mice "
            f"{sorted(set(rows.loc[dup, 'subjectID'].astype(str)))} (sampled twice "
            f"at one timepoint?); one row per mouse and transition is required"
        )
    verdict_of = dict(
        zip(rows["subjectID"].astype(str), rows["dominant_strain_change"])
    )
    names = [str(m) for m in carriers]
    if len(set(names)) != len(names):
        # baseline presence writes the carriers from a set; a repeat means a
        # hand-edited or concatenated rising table and would count a mouse twice
        repeated = sorted({m for m in names if names.count(m) > 1})
        raise ValueError(f"carrier(s) listed more than once: {repeated}")
    determined = changed = 0
    for mouse in names:
        verdict = _verdict(verdict_of.get(mouse))  # missing row -> undetermined
        if verdict is None:
            continue
        determined += 1
        changed += int(verdict)
    if determined < min_determined:
        return CARRIER_UNDETERMINED, determined, changed
    if changed * 2 > determined:  # strictly more than half
        return CARRIER_NEW_STRAIN, determined, changed
    return CARRIER_SAME_STRAIN, determined, changed
