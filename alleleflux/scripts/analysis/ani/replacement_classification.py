"""Roll per-mouse strain-background calls up to one row per MAG.

``strain_turnover.call_transitions`` gives one verdict per (mouse, transition).
The hypergeometric enrichment test works per MAG, so it needs ONE answer per
(MAG, diet group, transition): "did most of the mice swap strains?"  If yes the
MAG's significant sites are more likely linked-haplotype noise from a strain
replacement than parallel evolution, and the MAG is excluded from enrichment.

Every output row carries the SAME roll-up twice, over two units:

* the **mouse block**     (``n_mice_with_call`` ... ``no_mouse_changed``): each
  mouse is one vote.
* the **replicate block** (``n_replicates_with_call`` ...): each replicate is one
  vote, and a replicate counts as changed if ANY of its mice changed -- the
  rule for designs where a replicate is a cage of several mice.  In DRiDO
  the metadata has no replicate column, ``mag_metadata.py`` fills it with the
  mouse id, and the two blocks are identical; no toggle needed, the columns
  simply agree.

Both blocks are computed for BOTH metrics -- ``strain_replacement`` (popANI) and
``dominant_strain_change`` (conANI) -- and the ``metric`` column names which one
a row is about, so no consumer can silently favour one.  Groups are NEVER pooled.

Worked example (MAG_A, group 40, 5mo -> 22mo, metric strain_replacement, ``min_voters=1``):
    cage c1 = {m1 True, m2 False, m3 False}, cage c2 = {m4 True}, cage c3 = {m5 undetermined}
    mice:       with_call 4, undetermined 1, changed 2 -> majority 2 > 4/2 False, any True, all False, none False
    replicates: with_call 2 (c3 has no called mouse), changed 2 (c1 via m1, c2 via m4)
                -> majority True, all True
Undetermined mice are reported but never in either denominator.

Three states, not two
---------------------
A key is in one of THREE situations, and a True/False column can hold only two:
the strain was replaced, it was checked and was not, or too few mice could be
checked to say.  Two safeguards keep the third from passing as the second:

* **the voter floor** (``min_voters``, default 8).  With fewer called voters
  than that, every yes/no verdict of the block is left BLANK (pd.NA, an empty
  cell on disk) -- never ``False``, which would read as "checked, and it did
  not change".  The counts are always reported, so the facts stay visible.
* **the ``strain_status`` word** -- what a consumer should read.  One of
  (shown for the default ``vote_rule="majority"``)::

      replaced         >= min_voters called, MORE than half of them changed
      not_replaced     >= min_voters called, half or fewer changed
      too_few_voters   1 .. min_voters-1 called: a hint, not a verdict
      no_voters        nobody could be called
      tie              exactly half changed -- ONLY under ``tie="unresolved"``

  A consumer keeping ``strain_status == "not_replaced"`` cannot be fooled by a
  blank or by a status it has never seen; one asking "is the flag True?" can.

How votes become "replaced" is a setting, ``vote_rule``, decided HERE and not
by each consumer, so two scripts reading one table cannot disagree about it:

    majority (default)  more than half the called voters changed
    any                 at least one called voter changed   (cautious)
    all                 every called voter changed          (lenient)

The floor applies under every rule.  The yes/no columns are facts and never
depend on ``vote_rule``; only ``strain_status`` does.

An exact tie (4 of 8) is never a majority, so ``majority_*_changed`` is False
for it under every setting.  What a tie MEANS is a policy, set by ``tie``:
``"not_replaced"`` (default), ``"replaced"``, or ``"unresolved"`` (its own
status).  It matters under ``vote_rule="majority"`` only: the other two rules
cannot tie.  ``min_voters``, ``vote_rule`` and ``tie_rule`` are stamped on
every row, so a table says which rule produced it.

Real keys under the defaults (a sparsely covered longitudinal run,
dominant_strain_change): 30 called / 20
changed -> replaced; 9 / 4 -> not_replaced; 8 / 4 -> not_replaced (a tie); 5 / 3
-> too_few_voters, verdicts blank; 0 called -> no_voters, verdicts blank.
"""

import argparse
import glob
import logging
import os

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# The two per-mouse verdict columns produced by strain_turnover.call_transitions.
METRICS = ("strain_replacement", "dominant_strain_change")

KEYS = ["MAG_ID", "group", "transition"]

# Fewest called voters (mice, or replicates) a key needs before ANY verdict is
# given.  Eight matches the pipeline's own ``min_sample_num``: below it a
# "majority" can be one or two animals out of dozens sampled.
DEFAULT_MIN_VOTERS = 8
# How the called voters' verdicts become "replaced" (module docstring).
VOTE_RULES = ("majority", "any", "all")
DEFAULT_VOTE_RULE = "majority"
# What an exact tie (changed == called / 2) means for ``strain_status``.
TIE_RULES = ("not_replaced", "replaced", "unresolved")
DEFAULT_TIE = "not_replaced"

MOUSE_COLUMNS = (
    "n_mice_with_call",
    "n_mice_undetermined",
    "n_mice_changed",
    "majority_mice_changed",
    "any_mouse_changed",
    "all_mice_changed",
    "no_mouse_changed",
    "strain_status",  # the word a consumer should read (see the module docstring)
)
# any/none at replicate level would equal the mouse-level ones by construction
# (some mouse changed <=> some replicate changed), so only the counts, the
# majority and the all-changed verdicts are reported for replicates.
REPLICATE_COLUMNS = (
    "n_replicates_with_call",
    "n_replicates_changed",
    "majority_replicates_changed",
    "all_replicates_changed",
    "strain_status_replicates",  # same word, with replicates as the voters
)
# The rule that produced the row, stamped like strain_turnover stamps its thresholds.
PROVENANCE_COLUMNS = ("min_voters", "vote_rule", "tie_rule")


def _raw_answers(called: pd.Series, changed: pd.Series) -> dict[str, pd.Series]:
    """The four yes/no questions, answered ONCE, before any floor or policy.

    This is the single place the definitions of "majority", "any", "all" and
    "none" live.  Both outputs are derived from it: the yes/no columns
    (``_blank_below_floor``) and the status word (``_strain_status``), so the
    two can never drift apart.  The first three keys are deliberately the three
    ``VOTE_RULES``: the status word picks its answer by that name.

    Parameters
    ----------
    called
        Number of voters (mice, or replicates) that had a verdict, one value per
        (MAG_ID, group, transition) key.  E.g. ``9``.
    changed
        Number of those voters flagged as changed.  E.g. ``4``.

    Returns
    -------
    ``{"majority", "any", "all", "none"}`` -> plain boolean Series aligned with
    the inputs.  called=9, changed=4 -> majority False (8 is not > 9), any True,
    all False, none False.  NOT yet safe to report: called=5, changed=3 gives
    majority True here, and called=0 gives all True (0 == 0) and none True.
    The floor deals with those; nothing downstream may skip it.
    """
    return {
        # Strictly MORE than half of the called voters (4 of 8 is not a majority).
        # Written as 2 * changed > called so the comparison stays in whole numbers.
        "majority": changed * 2 > called,
        "any": changed > 0,
        "all": changed == called,
        "none": changed == 0,
    }


def _blank_below_floor(raw: dict[str, pd.Series], enough: pd.Series) -> dict[str, pd.Series]:
    """The yes/no columns as reported: the raw answers, BLANK where the floor is not met.

    Parameters
    ----------
    raw
        ``_raw_answers`` output.
    enough
        Boolean Series, True where ``called >= min_voters``.

    Returns
    -------
    The same four Series as nullable booleans.  With ``min_voters=8``: called=9,
    changed=4 -> unchanged (False, True, False, False).  called=5, changed=3 ->
    all four pd.NA: 3 of 5 is not allowed to read as a majority.  called=0 ->
    all four pd.NA: "no change" is a positive finding that needs evidence, and
    a blank -- unlike False -- cannot be mistaken for one.
    """
    # ``.where(enough)`` keeps the answer where the floor is met and blanks the
    # rest; the nullable "boolean" dtype is what lets a blank exist at all (plain
    # bool has no missing value and would silently turn it back into False).
    return {name: flag.astype("boolean").where(enough) for name, flag in raw.items()}


def _strain_status(
    called: pd.Series,
    changed: pd.Series,
    raw: dict[str, pd.Series],
    enough: pd.Series,
    vote_rule: str,
    tie: str,
) -> pd.Series:
    """The one-word status per key -- the column consumers should read.

    Takes the two counts, the ``_raw_answers`` computed from them, the floor
    mask ``enough`` (``called >= min_voters``), the vote rule (one of
    ``VOTE_RULES``) and the tie rule (one of ``TIE_RULES``).  It does not
    recompute any answer: ``replaced`` is ``raw[vote_rule]``.  Example,
    ``min_voters=8, vote_rule="majority", tie="not_replaced"``::

        called  changed  status
        30      20       replaced         (20 > 15)
        9       4        not_replaced     (4 <= 4.5)
        8       4        not_replaced     (a tie, by the tie rule)
        5       3        too_few_voters
        0       0        no_voters

    With ``tie="replaced"`` the 8 / 4 key reads ``replaced``; with
    ``tie="unresolved"`` it reads ``tie``.  No other key is affected.

    Under ``vote_rule="any"`` the 9 / 4 and 8 / 4 keys read ``replaced`` (someone
    changed); under ``"all"`` even 30 / 20 reads ``not_replaced`` (not everyone
    did).  The two bottom rows never change: the floor comes first.
    """
    # The chosen rule's answer, taken from the shared raw answers -- never
    # recomputed here (vote_rule is validated by classify_mags, so the key exists).
    replaced = raw[vote_rule]
    # Only a majority vote can tie: exactly half is neither a majority nor a
    # minority, so the tie rule decides.  "at least one" and "every one" cannot.
    is_tie = enough & (changed * 2 == called) & (vote_rule == "majority")
    tie_word = {"not_replaced": "not_replaced", "replaced": "replaced", "unresolved": "tie"}[tie]
    # Start from the safest word and overwrite upward, so a key no branch claims
    # can never end up labelled as verified.
    status = pd.Series("no_voters", index=called.index, dtype=object)
    status[(called > 0) & ~enough] = "too_few_voters"
    status[enough & ~replaced] = "not_replaced"
    status[enough & replaced] = "replaced"
    status[is_tie] = tie_word
    return status


def _judge(
    called: pd.Series, changed: pd.Series, min_voters: int, vote_rule: str, tie: str
) -> tuple[dict[str, pd.Series], pd.Series]:
    """Two counts per key -> ``(yes/no columns, status word)``, for one block of voters.

    The one sequence both blocks run (mice as voters, then replicates as
    voters): answer the four questions once, work out where the floor is met
    once, and derive BOTH outputs from those two things.

    Example, ``min_voters=8, vote_rule="majority", tie="not_replaced"``:
    called=9, changed=4 -> ({majority False, any True, all False, none False},
    ``not_replaced``); called=5, changed=3 -> (all four pd.NA, ``too_few_voters``).
    """
    raw = _raw_answers(called, changed)     # majority / any / all / none, computed once
    enough = called >= min_voters           # the floor, computed once
    return (
        _blank_below_floor(raw, enough),                               # the facts, blank under the floor
        _strain_status(called, changed, raw, enough, vote_rule, tie),  # the decision, from the same answers
    )


def _classify_one_metric(
    turnover: pd.DataFrame, metric: str, min_voters: int, vote_rule: str, tie: str
) -> pd.DataFrame:
    """Roll one verdict column up to one row per (MAG_ID, group, transition).

    Parameters
    ----------
    turnover
        The ``strain_turnover.call_transitions`` output: one row per mouse per
        transition (970 rows for MRGM_0841 on DRiDO), in the turnover-file
        schema.  Columns used here: ``MAG_ID``, ``group``, ``replicate``,
        ``transition``, and the two nullable-boolean verdict columns
        ``strain_replacement`` and ``dominant_strain_change`` (True = changed,
        False = same strain, pd.NA = undetermined).  Everything else rides along
        unused.
    metric
        WHICH verdict column to count: ``"strain_replacement"`` (the popANI
        verdict) or ``"dominant_strain_change"`` (the conANI verdict).  Nothing
        else in the row is read as evidence.
    min_voters, vote_rule, tie
        The voter floor, the vote rule and the tie rule (module docstring).  The
        SAME three apply to both blocks, counting mice in one and replicates in
        the other.

    Returns
    -------
    One row per (MAG_ID, group, transition) with the columns ``KEYS +
    ["metric"] + MOUSE_COLUMNS + REPLICATE_COLUMNS + PROVENANCE_COLUMNS``.
    Diet groups are separate keys and never pooled.

    Worked example, metric="strain_replacement", ``min_voters=1``, MAG_A / group
    fat / pre_end, five mice in three cages::

        mouse  cage  strain_replacement
        m1     c1    True
        m2     c1    False
        m3     c1    False
        m4     c2    True
        m5     c3    <NA>

    mouse block:     n_mice_with_call 4, n_mice_undetermined 1, n_mice_changed 2
                     -> majority False, any True, all False, none False
    replicate block: c1 changed (via m1), c2 changed, c3 has no called mouse
                     -> n_replicates_with_call 2, n_replicates_changed 2
                     -> majority True, all True
    statuses:        mice 2 of 4 is a tie -> ``not_replaced`` under the default
                     tie rule; replicates 2 of 2 -> ``replaced``.
    With the production ``min_voters=8`` the same five mice give blank verdicts
    and ``too_few_voters`` in both blocks.
    """
    df = turnover
    # ``metric`` is the NAME of the verdict column, e.g. "strain_replacement";
    # ``flag`` is that column's values (True / False / <NA>), one per mouse.
    flag = df[metric].astype("boolean")
    # Per-mouse 0/1 helper columns that sum cleanly under groupby.
    df = df.assign(
        _called=flag.notna().astype(int),
        _undetermined=flag.isna().astype(int),
        _changed=flag.fillna(False).astype(int),
    )

    # --- mouse block: every row is one mouse, so sums are mouse counts.
    mice = df.groupby(KEYS, sort=True, observed=True)[
        ["_called", "_undetermined", "_changed"]
    ].sum()
    mouse_verdicts, mouse_status = _judge(
        mice["_called"], mice["_changed"], min_voters, vote_rule, tie
    )
    mice = mice.rename(
        columns={
            "_called": "n_mice_with_call",
            "_undetermined": "n_mice_undetermined",
            "_changed": "n_mice_changed",
        }
    )
    mice["majority_mice_changed"] = mouse_verdicts["majority"]
    mice["any_mouse_changed"] = mouse_verdicts["any"]
    mice["all_mice_changed"] = mouse_verdicts["all"]
    mice["no_mouse_changed"] = mouse_verdicts["none"]
    mice["strain_status"] = mouse_status

    # --- replicate block: first collapse mice to replicates (a replicate has a
    # call if ANY mouse has one, changed if ANY mouse changed), then count the
    # replicates exactly like mice above.
    per_replicate = df.groupby(KEYS + ["replicate"], sort=True, observed=True)[
        ["_called", "_changed"]
    ].max()
    reps = per_replicate.groupby(KEYS, sort=True, observed=True).sum()
    rep_verdicts, rep_status = _judge(
        reps["_called"], reps["_changed"], min_voters, vote_rule, tie
    )
    reps = reps.rename(
        columns={
            "_called": "n_replicates_with_call",
            "_changed": "n_replicates_changed",
        }
    )
    reps["majority_replicates_changed"] = rep_verdicts["majority"]
    reps["all_replicates_changed"] = rep_verdicts["all"]
    reps["strain_status_replicates"] = rep_status

    out = mice.join(reps).reset_index()
    out.insert(len(KEYS), "metric", metric)
    # Which rule made this row: a consumer can check it got the table it expects.
    out["min_voters"] = min_voters
    out["vote_rule"] = vote_rule
    out["tie_rule"] = tie
    return out[
        KEYS + ["metric"] + list(MOUSE_COLUMNS) + list(REPLICATE_COLUMNS) + list(PROVENANCE_COLUMNS)
    ]


def classify_mags(
    turnover: pd.DataFrame,
    min_voters: int = DEFAULT_MIN_VOTERS,
    vote_rule: str = DEFAULT_VOTE_RULE,
    tie: str = DEFAULT_TIE,
) -> pd.DataFrame:
    """The enrichment filter's input: both metrics rolled up per MAG, stacked.

    Parameters
    ----------
    turnover
        The ``call_transitions`` output, one row per mouse per transition (see
        ``_classify_one_metric`` for the columns read).
    min_voters
        Fewest called voters a key needs for any verdict (default 8).  Must be
        >= 1: zero would let a key nobody could check count as assessed.
    vote_rule
        One of ``VOTE_RULES``: ``"majority"`` (default), ``"any"`` or ``"all"``.
    tie
        One of ``TIE_RULES``: what an exact half-and-half vote means (used under
        ``vote_rule="majority"`` only).

    Raises ``ValueError`` on a bad ``min_voters``, ``vote_rule`` or ``tie`` -- a
    typo here would otherwise quietly change which MAGs an analysis keeps.

    Returns
    -------
    ``_classify_one_metric`` run for each of ``METRICS`` and concatenated, so
    every (MAG_ID, group, transition) key appears TWICE: once with
    ``metric == "strain_replacement"`` and once with
    ``metric == "dominant_strain_change"``.  Consumers filter on ``metric``;
    there is deliberately no way to compute only one.

    Real output for MRGM_0841, group 2D, 5mo_10mo (12 called, 12 undetermined)::

        metric                  n_mice_changed  majority_mice_changed  strain_status
        strain_replacement      12              True                   replaced
        dominant_strain_change   3              False                  not_replaced
    """
    if min_voters < 1:
        raise ValueError(f"min_voters must be >= 1, got {min_voters}")
    if vote_rule not in VOTE_RULES:
        raise ValueError(f"vote_rule must be one of {VOTE_RULES}, got {vote_rule!r}")
    if tie not in TIE_RULES:
        raise ValueError(f"tie must be one of {TIE_RULES}, got {tie!r}")
    stacked = pd.concat(
        [_classify_one_metric(turnover, metric, min_voters, vote_rule, tie) for metric in METRICS],
        ignore_index=True,
    )
    logger.info(
        f"MAG-level classification: {len(stacked) // len(METRICS)} keys x {len(METRICS)} metrics "
        f"(min_voters={min_voters}, vote_rule={vote_rule}, tie={tie})"
    )
    # How much of the table each status holds, per metric: the first thing to
    # look at, since "too_few_voters" / "no_voters" are usually the majority.
    for metric in METRICS:
        counts = stacked.loc[stacked["metric"] == metric, "strain_status"].value_counts().to_dict()
        logger.info(f"  {metric}: {counts}")
    return stacked


# ---------------------------------------------------------------------------
# The ``alleleflux-replacement-classification`` command: all MAGs -> one table
# ---------------------------------------------------------------------------

# Columns of a {mag}_strain_turnover.tsv that carry ids: read as str so DRiDO's
# numeric group names ("20", "40") never sniff to int and break equality tests.
_ID_COLUMNS = (
    "MAG_ID",
    "subjectID",
    "replicate",
    "group",
    "transition",
    "sample_t1",
    "sample_t2",
)


def load_turnover_dir(turnover_dir: str) -> pd.DataFrame:
    """Stack every ``*_strain_turnover.tsv`` in a directory into one frame.

    Parameters
    ----------
    turnover_dir
        Directory holding the ``alleleflux-strain-turnover`` outputs, one file
        per MAG.  Header-only files (a MAG with no mouse matching any
        transition) are read and contribute zero rows -- not an error.

    Returns
    -------
    The concatenated per-mouse table; the two verdict columns are restored to
    the nullable "boolean" dtype (a blank cell on disk = undetermined = pd.NA).
    Raises ``FileNotFoundError`` if the directory holds no turnover files at
    all -- an empty classification is never a legitimate outcome.

    Example: files for MAG_A (3 rows), MAG_B (1 row), MAG_C (header only)
    -> 4 rows, MAG_ID in {MAG_A, MAG_B}.
    """
    paths = sorted(glob.glob(os.path.join(turnover_dir, "*_strain_turnover.tsv")))
    if not paths:
        raise FileNotFoundError(f"no *_strain_turnover.tsv files under {turnover_dir}")
    frames = [
        pd.read_csv(path, sep="\t", dtype={c: str for c in _ID_COLUMNS})
        for path in paths
    ]
    turnover = pd.concat(frames, ignore_index=True)
    for col in METRICS:
        # Written as True/False/blank; blank reads back as NaN (float) -> nullable bool.
        turnover[col] = turnover[col].astype("boolean")
    logger.info(f"{len(paths)} turnover files -> {len(turnover)} mouse-transition rows")
    return turnover


def round_up_the_verdicts(args: argparse.Namespace) -> int:
    """Orchestrator: turnover directory -> classify_mags -> one TSV.

    The turnover table is the classifier's native input: ``call_transitions``
    already verified both sides of every pair agree on group and replicate
    before collapsing them to one column each.
    """
    turnover = load_turnover_dir(args.turnover_dir)
    classified = classify_mags(
        turnover, min_voters=args.min_voters, vote_rule=args.vote_rule, tie=args.tie
    )
    os.makedirs(os.path.dirname(os.path.abspath(args.output_path)), exist_ok=True)
    classified.to_csv(args.output_path, sep="\t", index=False)
    logger.info(
        f"wrote {len(classified)} rows ({classified['MAG_ID'].nunique()} MAGs x groups x transitions "
        f"x {len(METRICS)} metrics) to {args.output_path}"
    )
    return 0


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Roll per-mouse strain calls up to one row per MAG, group, transition and metric.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--turnover_dir",
        required=True,
        help="Directory of {mag}_strain_turnover.tsv files from alleleflux-strain-turnover",
    )
    parser.add_argument(
        "--output_path",
        required=True,
        help="Output TSV (the enrichment filter's input)",
    )
    parser.add_argument(
        "--min_voters",
        type=int,
        default=DEFAULT_MIN_VOTERS,
        help=(
            "Fewest called mice (or replicates) a MAG x group x transition needs "
            "for any verdict. Below it the yes/no columns are left blank and "
            "strain_status reads too_few_voters / no_voters"
        ),
    )
    parser.add_argument(
        "--vote_rule",
        choices=list(VOTE_RULES),
        default=DEFAULT_VOTE_RULE,
        help=(
            "When a checked key counts as replaced: majority (more than half the "
            "called voters changed), any (at least one did), all (every one did)"
        ),
    )
    parser.add_argument(
        "--tie",
        choices=list(TIE_RULES),
        default=DEFAULT_TIE,
        help=(
            "What an exact half-and-half vote means for strain_status (used with "
            "--vote_rule majority only): not_replaced, replaced, or unresolved "
            "(its own status, 'tie')"
        ),
    )
    args = parser.parse_args()
    return round_up_the_verdicts(args)


if __name__ == "__main__":
    raise SystemExit(main())
