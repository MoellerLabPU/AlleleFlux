# Strain Turnover and Baseline Presence

This guide covers the four commands that answer two questions the statistical
tests cannot: **did a mouse's strain of a MAG change between timepoints?** and
**was a significant allele already present at baseline?**

```
profiles ──► alleleflux-pairwise-ani ──► alleleflux-strain-turnover ──► alleleflux-replacement-classification
                                                    │
p_value_summary ──────────────────────────► alleleflux-baseline-presence
```

All four are wired into the workflow behind three flags (`use_pairwise_ani`,
`use_strain_turnover`, `use_baseline_presence`) and can also be run by hand on
any completed run.  Longitudinal data only.

## Why these exist

A significant site says that an allele's frequency differs between groups or
changed within a group.  Two very different histories produce the same p-value:

* **Selection on standing variation.**  The allele was already in the mouse at
  baseline and rose.
* **Strain replacement or a new mutation.**  A different strain of the same MAG
  arrived, or the allele arose, after baseline.

Strain turnover detects the second case at the whole-genome level (was the
population in this mouse the *same strain* at both timepoints?).  Baseline
presence checks it at the site level (was *this allele* there at baseline?).
An enrichment analysis over significant sites should exclude MAGs whose strain
background changed, otherwise the sites reflect strain identity rather than
evolution within a strain.

## Step 1: pairwise ANI (`alleleflux-pairwise-ani`)

One job per MAG.  For every pair of samples that passed QC it computes, over
the positions covered by at least `min_cov` reads in **both** samples:

| metric | counts a position as different when | reads |
|---|---|---|
| **conANI** | the majority base differs between the two samples | the dominant strain swapped |
| **popANI** | the two samples share **no** allele at the position | no strain in common at all |

conANI is always the lower of the two.  "Present" for popANI means at least the
null-model bar for that depth **and** at least `min_freq` of reads; the bar
comes from a per-depth error model (`fdr`, `min_base_quality`).

`pairs` chooses which pairs are compared: `within_subject` (every same-mouse
pair), `transitions` (only same-mouse pairs matching the configured timepoint
combinations; what strain turnover consumes), or `all` (every pair, for
between-mouse strain-sharing questions).  On single-timepoint data every mouse
has one sample, so there are no same-mouse pairs: `all` is the default there
and the only accepted value, and strain turnover itself (which needs a
transition) is refused at parse time together with the other
longitudinal-only features (see `data_type` in the configuration reference).

**Which samples take part.** A MAG's ANI job reads the MAG's QC tables at the
timepoints of the comparisons the MAG is eligible in, and keeps the QC-passing
samples of the groups named in `groups_combinations`.  This does not move with
the de novo settings below: samples that are read only as baseline witnesses
(another timepoint, or a group outside every comparison) never enter an ANI
pair.

The filter-by-filter reference, with a worked position, lives with the code at
`alleleflux/scripts/analysis/ani/README.md`.

## Step 2: strain turnover (`alleleflux-strain-turnover`)

One job per MAG, seconds each (a local rule in the workflow).  For every mouse
and every configured transition (`EARLIER:LATER`) it takes that mouse's pair
from the ANI table and issues two verdicts:

| verdict | rule | reading |
|---|---|---|
| `strain_replacement` | popANI < `pop_threshold` (default 0.99999) | the strain was replaced, or an invisible-at-baseline strain took over |
| `dominant_strain_change` | conANI < `con_threshold` (default 0.999) | the dominant strain swapped even if the old one lingers |

A verdict is only issued when the pair compared at least `min_compared` of the
genome (default 0.1) and at least one position; otherwise both verdicts are
`<NA>` (undetermined), never a silent False.  A `background` label summarises
the two: `stable`, `strain_replacement`, `dominant_strain_change`, or
`undetermined`.

Outputs per MAG: `{mag}_strain_turnover.tsv` (one row per mouse and
transition, with both ANI values, the fraction compared, and both verdicts) and
`{mag}_turnover_rollup.tsv` (counts per group and transition).

## Step 3: replacement classification (`alleleflux-replacement-classification`)

Runs once over every MAG's turnover table and writes one file,
`replacement_classification.tsv`: one row per MAG, transition and metric, with
a **mouse block** (how many mice were called, how many changed, whether all
did) and a **replicate block** (a replicate counts as changed if any of its
mice changed, for designs where a replicate is a cage of several mice; where
the metadata has no replicate column, replicate equals subject and the two
blocks agree).  Both metrics are always reported, stacked, with a `metric`
column, so the choice of threshold is made downstream.

This table is the input for a strain-aware enrichment filter.  Filter on the
`strain_status` column and **keep** the MAG × group × transition keys that read
`not_replaced`, before counting significant sites.

Do not write the filter the other way round ("drop the keys where the
background changed").  A key is in one of three situations:

| `strain_status` | meaning |
|---|---|
| `replaced` | at least `min_voters` mice had a verdict and more than half changed |
| `not_replaced` | at least `min_voters` mice had a verdict and half or fewer changed |
| `too_few_voters`, `no_voters` | not enough mice could be checked to say |

In sparsely covered data most keys are in the third row, because a mouse only
gets a verdict when both of its samples cover enough of the genome.  Dropping
the `replaced` keys keeps all of those unchecked keys as if they had passed.
For the same reason the yes/no columns are left **blank** below the floor
instead of `False`: a blank cannot be read as "checked, and it did not change".

Three settings under `strain_turnover` in the config shape the status:
`min_voters` (default 8), `vote_rule` (`majority` by default; `any` flags a key
when a single voter changed, `all` only when every voter did) and `tie` (what
an exact half-and-half vote means under `majority`, default `not_replaced`).
The rule is chosen here, once, and stamped on every row, so every analysis
reading the table uses the same one.  All three belong to this step alone:
changing them reruns only this roll-up, not pairwise ANI or the per-mouse
calls.

A fourth setting, `replicate_rule`, matters only when a replicate holds several
subjects, for example a cage of co-housed mice.  It decides how those subjects
become the replicate's single vote in the replicate block:

| `replicate_rule` | The replicate counts as changed when |
|---|---|
| `average` (default) | the mean ANI of its subjects with a verdict is below the threshold: conANI against `con_threshold`, popANI against `pop_threshold` |
| `any` | at least one of its subjects with a verdict changed |
| `majority` | more than half of them changed (1 of 2 is not a majority) |

Subjects without a verdict never take part, under any rule.  Their rows still
carry ANI numbers, computed from too little of the genome to trust, and
`average` leaves those out of the mean.  A replicate with a single subject gives
the same answer under all three rules.

## Step 4: baseline presence (`alleleflux-baseline-presence`)

One job per `{timepoints}-{groups}` comparison, after the statistics.  It reads
the comparison's `p_value_summary` file for one test family, keeps the
significant sites (`threshold_column` ≤ `threshold`), picks the candidate
allele(s) at each site (at a site where two bases were tested, both: their
frequencies are complementary, so the two tests are one test and which base
rose is decided later from the frequency change; at a site with three or four
tested bases, only those whose own p-value is at or below the comparison's
significance line, reported with `n_candidate_alleles`), then opens **every**
sample's
profile for both groups and both timepoints and gives each (site, allele,
sample) one of four statuses:

| status | meaning |
|---|---|
| `present` | covered, the allele clears the bar and `min_freq` |
| `below_detection` | covered, some reads of the allele but under the bar |
| `absent` | covered, zero reads of the allele |
| `not_covered` | fewer than `min_cov` reads, or the site-depth filter removed the position for that sample; the cache cannot tell the two apart (a `filtered` verdict is planned). Samples that failed QC for the MAG have no rows at all |

The presence rule is the same one pairwise ANI uses, so the two agree by
construction; the workflow passes the same `min_cov`, `min_freq` and `fdr`.

### Which samples are the baseline (`baseline_scope`)

"Absent at baseline" is only meaningful once you say *which* earlier samples
count.  `analysis.baseline_presence.baseline_scope` (CLI `--baseline_scope`)
picks one of four answers.  In plain words, for an allele G that rose in the
fat group of a fat-vs-control comparison, where mice m1 and m2 (fat) carry G
at the later timepoint and m3 (control) does not:

| value | look for G in | intended design |
|---|---|---|
| `individual` | only the carriers' own earlier samples, m1's and m2's. If neither had any G, it is new. m3's earlier sample is not consulted, even if it had a whisper of G, because m3 is not a carrier. | isolated populations that never share alleles |
| `own_group` | every fat mouse's earlier sample (within-group test families only: `single_sample`, `lmm_across_time`; an error otherwise) | groups housed apart |
| `comparison` | every earlier sample of both groups, fat and control | one population split across groups |
| `cohort` (default) | every QC-passing sample of the whole run, any group, any age, as long as it was collected on an earlier **date** than the first carrier's sample. If a completely unrelated mouse had G a year earlier, G was already in the building: not new. | shared facility or colony |

**What `cohort` needs.** A `collection_date` column in the metadata, ISO
`YYYY-MM-DD` (see [inputs](../reference/inputs.md)).  The job then loads the
comparison's samples plus every other sheet sample that may testify (by
default all of them; see `cohort_witnesses` below), and for each site keeps as
baseline those dated strictly before the earliest later-timepoint sample that
carries the allele (the *first carrier*).  Samples dated on or after it, the
carriers' own later samples included, are not baseline.  A `cohort` column
(enrollment batch), when present, is carried into the outputs.

**What `cohort` builds.** A sample can only be read as a baseline witness if
it has QC and an allele-frequency cache.  So when all four of these hold:
longitudinal data, `baseline_scope: cohort`, a `collection_date` for every
sample, and baseline presence or the gene-level de novo analysis switched on,
the workflow builds a cache for every eligible MAG at **every witness
timepoint**, not only at the timepoints of the comparisons it is eligible
in.  ("Eligible" = the eligibility table of at least one comparison admits
the MAG to at least one test.)  When one of the four does not hold, nothing
extra is built.

Which timepoints and groups are witnesses is one setting,
`analysis.baseline_presence.cohort_witnesses`:

| `cohort_witnesses` | Who may testify (if dated before the first carrier) | Extra work |
|---|---|---|
| `sheet` (default) | any sample of the metadata sheet, at any timepoint and in any group | sample metadata, QC and caches also at the sheet's timepoints that sit in no comparison, and also for the sheet's groups that sit in no comparison |
| `comparisons` | only samples at a timepoint named in `timepoints_combinations` **and** in a group named in `groups_combinations` | no extra timepoint or group; each eligible MAG still gets a cache at every compared timepoint |

Example.  The sheet has timepoints `pre`, `mid`, `post`, `end` and groups `A`,
`B`, `C`.  The config compares `pre_post` and `pre_end`, group `A` against `B`.

| | witness timepoints | witness groups | QC folders | caches of a MAG eligible in `pre_post` only |
|---|---|---|---|---|
| `sheet` | `pre`, `post`, `end`, `mid` | `A`, `B`, `C` | `pre`, `post`, `end`, `mid` | `pre`, `post`, `end`, `mid` |
| `comparisons` | `pre`, `post`, `end` | `A`, `B` | `pre`, `post`, `end` | `pre`, `post`, `end` |
| the four conditions not met | (none outside the comparison) | (none) | `pre`, `post`, `end` | `pre`, `post` |

Why groups matter: an allele rises in `A`, and a `C` sample dated earlier
carried it.  With `C` left out the allele is called de novo.  Under `sheet`
the `C` sample testifies and the allele is standing variation.

The extra timepoints and groups are only ever read as witnesses.  They are
never tested: eligibility and the statistics name their own two groups and two
timepoints, and pairwise ANI keeps to the compared timepoints and groups.  The
workflow turns the setting into two lists (witness timepoints, witness groups),
builds QC and caches for exactly those, and hands the same lists to both de
novo commands (`--witness_timepoints`, `--witness_groups`), so a command never
asks for a cache that was not built.

A MAG with no QC-passing sample at a witness timepoint still gets a cache: a
file with the normal columns and **no rows**.  The de novo commands read it as
"no witness of this MAG at this timepoint".  A cache file that is missing
altogether stops them with an error naming the path.

A timepoint label that sits in no comparison becomes a folder name, so under
`sheet` it may hold only letters, digits, `.` and `-`.  A group label that
sits in no comparison may hold letters, digits, `.`, `-` and `_`.  The run
stops while the workflow is being read, before any job, on a label that breaks
these rules and on a blank `time` or `group` cell in the sheet.

**Dates are all-or-nothing.** A `collection_date` column with *some* blank
cells is refused under `cohort`: the job stops before writing anything and
names every undated sample.  Ordering the baseline by timepoint label while
some samples are dated would silently change the answer, so it is treated as a
data error.  (Under every other scope blank dates are fine; they are only
carried along.)  A sheet with **no** dates at all, a study that never recorded
them, still runs, but nothing can then place a third timepoint before or after
the comparison, and the order of rows in the sheet says nothing.  The baseline
narrows to the comparison's earlier timepoint, exactly as under `comparison`,
the other timepoints are not loaded, and the run logs one warning saying so.

**What stays the comparison's, whatever the scope.** The later-timepoint
checks and counts, the frequency change, and the rising table use only the
comparison's own samples (long-table column `is_comparison_sample`).  Samples
loaded for `cohort` are baseline material only; they never get a rising row.

**Naming under `cohort`.** The baseline columns keep the comparison's earlier
label (`n_pre_samples_covered`, `total_reads_pre`, `pre_mice_allele_present`)
even when the baseline includes samples from other timepoints.  Four columns
give the calendar evidence: `first_carrier_date` and `first_carrier_cohort`
(the earliest-dated later sample of the row's group that carries the allele,
and its enrollment batch) and `first_detection_date` and
`first_detection_cohort` (the earliest-dated sample among those judged,
baseline and later, in which the allele is present, and its batch).  The two
pairs coincide when nobody earlier had the allele; the gap between them is
how long it was in the colony before these mice carried it.

**How hard we looked.** Two reporting-only columns on the summary and rising
tables, filtered on by nothing: `n_baseline_reads`, the reads at the position
over the covered baseline samples (a de novo call on 900 reads and one on 20
are not the same claim), and `n_carriers_own_baseline_absent`, how many of the
carriers have an earlier sample of their own that was covered and clean (the
pooled "nobody had it" can rest entirely on other mice).  Sort the lone-allele
candidates by them; they never change a verdict.

**Where the baseline is cut.** Each row is judged against the samples dated
strictly before its *cutoff*: the first carrier's date, or, for a group that
rose without any later sample clearing the presence bar (no carrier, verdict
`allele_not_present_at_{later}`), the day that group was first sampled at the
later timepoint.  So the two rising rows of one site, and the summary row,
describe comparable baselines; a sheet with no dates at all falls back to the
comparison's earlier timepoint.

### Three labels (`pre_treatment_timepoints`)

Under the dated `cohort` scope an allele is "not de novo" as soon as any
earlier-dated sample carried it.  But some earlier-dated samples were taken
**after treatment began**.  An allele first seen in one of those may have
arisen under treatment, which is exactly the event the analysis looks for.

`analysis.baseline_presence.pre_treatment_timepoints` is a list of timepoint
labels (as spelled in the sheet's `time` column) that were sampled before
treatment began.  The default is an empty list: no split, two labels as
before.  With the list filled, "an earlier sample carried the allele" splits
in two:

| Witnesses that carry the allele | Label |
|---|---|
| none | de novo (unchanged) |
| at least one pre-treatment witness | `standing_variation` |
| only witnesses sampled after treatment began | `first_seen_after_treatment_began` |

Two things do **not** change:

* **Who is a witness.**  Still the covered samples dated strictly before the
  first carrier.  The list never adds or removes a witness.
* **How a witness is classed.**  By its timepoint label only: a witness at a
  listed timepoint is pre-treatment, a witness at any other timepoint is
  on-treatment.  Its group and its date play no part.

Example.  Timepoints `t0`, `t1`, `t2`; the comparison is `t0_t2`; the config
says `pre_treatment_timepoints: ["t0"]`.  Allele G rose, and its first carrier
is dated 2021-06-10.

| sample | timepoint | collected | covered | carries G | witness? | class |
|---|---|---|---|---|---|---|
| w1 | `t0` | 2021-01-10 | yes | no | yes | pre-treatment |
| w2 | `t0` | 2021-01-12 | yes | no | yes | pre-treatment |
| w3 | `t1` | 2021-04-02 | yes | yes | yes | on-treatment |
| w4 | `t0` | 2021-07-01 | yes | yes | no: dated after the first carrier | not counted |

| question | answer |
|---|---|
| Did any witness carry G? | yes, w3: not de novo |
| Did a pre-treatment witness carry G? | no, w1 and w2 are clean |
| Label | `first_seen_after_treatment_began` |
| Witness counts | pre-treatment: 2 covered, 0 carrying, 0 below detection; on-treatment: 1 covered, 1 carrying |

Had w1 carried G, the label would be `standing_variation`.  Without the list
the label is `standing_variation` and the five count columns are empty.  w4
shows a consequence of the calendar rule: a pre-treatment sample dated after
the first carrier is not a witness, so "no pre-treatment witness carried it"
can rest on few samples.  That is why the counts are reported.

**"First seen after treatment began" is a label, not a proof.**  It also holds
alleles that were in the population all along but in no sampled pre-treatment
witness.  Read it together with `n_pre_treatment_witnesses_covered`: a call
resting on 2 pre-treatment witnesses is much weaker than one resting on 80.  A
pre-treatment witness with a few reads of the allele under the detection bar
does not carry it by the presence rule, so it does not make the allele
standing; `n_pre_treatment_witnesses_below_detection` flags those rows.

**Where the label appears.**

| Output | What changes with the list filled |
|---|---|
| Baseline-presence summary and rising table | `origin_any_mouse` can read `first_seen_after_treatment_began`; the five witness-count columns are filled. The columns are always in the header and are empty when the list is empty. The long table is not affected |
| `rising_origin_alleles.tsv` | the same five columns, and `origin` can read the new label |
| `rising_origin_summary.tsv` | a `first_seen_after_treatment_began` row per block, only when the run made the split. The label is a verdict, so it is counted in `n_called`. The cluster rule and the strain gate never relabel it |
| Gene-level events table | the label appears in the existing `origin` column. It is never an event, so gene counts, cage counts and both tests are identical with and without the list |

**What stops the run.**  Whenever baseline presence or the gene-level de novo
analysis is on, the list is checked while the workflow is being read, before
any job is submitted.  The run stops when:

* the list is filled and `baseline_scope` is not `cohort`;
* the list is filled and the sheet does not have a `collection_date` for every
  sample (without dates every witness comes from the comparison's earlier
  timepoint, so no witness is "after treatment began");
* the value is not a list (`t0` instead of `["t0"]`), or a label is blank;
* a label is not a timepoint of the sheet;
* a label is a sheet timepoint but not a witness timepoint (possible only
  under `cohort_witnesses: comparisons`, for a timepoint in no comparison).

Run by hand, `alleleflux-baseline-presence` and `alleleflux-gene-de-novo`
refuse `--pre_treatment_timepoints` in the same three cases: a scope other
than `cohort`, a sheet without a date for every sample, and a label that is no
timepoint of the samples they read.

### Naming

Every column or label that refers to a timepoint uses the comparison's own
labels, exactly as the metadata spells them: `n_pre_samples_covered`,
`allele_absent_at_end`, `de_novo_candidate_below_detection_at_5mo`.  Nothing
is hard-coded.

### Outputs

**Long table**, `{comparison}_{family}_{statistic}_baseline_presence.tsv.gz`:
one row per site × allele × sample with `allele_reads`, `total_reads`,
`detection_threshold_reads`, `allele_frequency`, `allele_status`, and
`origin_in_own_mouse`, the per-mouse verdict on later-timepoint rows:

| later status | own earlier status | label (shown for `pre`/`end`) |
|---|---|---|
| present | present | `standing_variation` |
| present | below_detection | `de_novo_candidate_below_detection_at_pre` |
| present | absent | `de_novo_candidate` |
| present | not_covered | `pre_not_covered` |
| present | no earlier sample | `no_pre_sample` |
| below_detection | any | `allele_below_detection_at_end` |
| absent | any | `allele_absent_at_end` |
| not_covered | any | `end_not_covered` |

With `--turnover_dir`, a `strain_background` column joins each mouse's strain
verdict for the same transition.

Two optional metadata columns ride along onto every row of their sample and are
blank when the metadata lacks them: `collection_date`, the calendar date the
sample was collected, accepted in **ISO `YYYY-MM-DD` form only** (`2017-01-10`;
`10/01/2017` or `2017-1-10` stop the job with an error naming the sample, and an
empty cell means unknown), and `cohort`, a free-text enrollment batch label.

**Summary**, `..._baseline_presence_summary.tsv`: one row per site × allele.
Its columns come in two kinds:

* **Sample counts, filtered by the presence rule:** `origin_any_mouse`
  (`standing_variation` if any covered earlier sample has the allele present,
  else `de_novo_candidate_below_detection_at_{earlier}`, else
  `de_novo_candidate`, else `{earlier}_not_covered`; `allele_not_present_at_{later}`
  when no covered later sample shows the allele), `n_{earlier}_samples_allele_present`,
  `n_{earlier}_samples_covered`, `n_replicates_with_allele_at_{earlier}`,
  `{earlier}_mice_allele_present`, and the same-mouse counts
  `n_mice_standing_variation`, `n_mice_de_novo_candidate`,
  `n_mice_de_novo_candidate_below_detection_at_{earlier}`.  Five witness-count
  columns (`n_pre_treatment_witnesses_covered`,
  `n_pre_treatment_witnesses_allele_present`,
  `n_pre_treatment_witnesses_below_detection`,
  `n_on_treatment_witnesses_covered`,
  `n_on_treatment_witnesses_allele_present`) are always present and are filled
  only when `pre_treatment_timepoints` is set; `origin_any_mouse` can then also
  read `first_seen_after_treatment_began` (see
  [Three labels](#three-labels-pre_treatment_timepoints)).
* **Read counts, unfiltered:** `total_reads_{tp}`, `allele_reads_{tp}`,
  `allele_frequency_{tp}` for each timepoint, summed over **every** sample at
  that timepoint, thin and profile-less samples included.  That is what makes
  "never seen in N reads at baseline" a frequency bound of 1/N.  The two kinds
  can disagree on purpose: a `de_novo_candidate` site can have a few baseline
  reads of the allele sitting in samples under the detection bar.

**Rising table**, `..._baseline_presence_rising.tsv`: Rising table: the summary's columns plus `rose_in`, `frequency_change`, `n_replicates`, `n_mice`; one row per (site, allele, group the allele rose in). At a biallelic site both alleles tie on the p-value and both are listed in the summary, but only one rose; the one that fell was the common allele at baseline and always reads standing. The rising table leaves it out, and checks end-point presence in the `rose_in` group only. Count standing vs de novo from this table.

"Rose" means `frequency_change > 0`, computed the way the within-group test
measures change: each mouse sampled at both timepoints, later minus earlier
frequency; the mean per replicate; the mean over replicates.  `n_replicates` and
`n_mice` say how many replicates and mice that rests on.  A within-group site is
asked about its own group only; a between-group site about both groups, so an
allele that rose in both gets one row per group.  Alleles with no mouse sampled
at both timepoints have no change and no row; the run log counts them.

**Across comparisons:** `alleleflux-rising-origin` tallies every rising
table into one summary (see [Outputs](../reference/outputs.md)).  Before
counting it applies two rules to the de novo alleles, in this order:

1. **The cluster rule.** De novo alleles within `cluster_window_bp` (default
   1000) of each other on one contig, in one group, at any age, are one
   `block_event`, never independent mutations; an allele with more than
   `max_de_novo_neighbours` (default 0) such neighbours is relabelled.  This
   runs first because adjacency is a fact about positions: a strain gate that
   ran first could relabel some members and hide a block member's neighbours.
2. **The strain gate.** A de novo allele most likely arrived with a new strain
   rather than by mutation when the strain changed.  Two verdicts exist.  The
   *group* verdict is the across-cage vote of Step 3 for the allele's MAG,
   group and period.  The *carrier* verdict asks only the mice carrying the
   allele at the later timepoint (`{later}_mice_allele_present` in the rising
   table): each carrier's own per-mouse `dominant_strain_change` from Step 2's
   `{mag}_strain_turnover.tsv`; more than half of the determined carriers
   changed → `new_strain_carriers`; fewer determined than
   `min_determined_carriers` (default 3) → `strain_undetermined`.
   `strain_gate` picks which decides: `carrier` (default), `group`, `both`
   (one undetermined verdict makes the allele undetermined; otherwise one
   swap verdict makes it `new_strain`) or `none`.  Under every gate but `none`, an allele whose
   strain could not be determined, by its carriers or by a group vote that is
   neither replaced nor not replaced, becomes `strain_undetermined` and is not
   counted as de novo.  Both verdicts are written to every row for comparison.

Neither rule touches an allele labelled `first_seen_after_treatment_began`:
both act on the two de novo origins only.

The workflow runs it once, after the last baseline-presence job, with the
run's classification table and turnover tables; by hand it takes one
classification table per voter floor and the turnover folder.

### Worked example (synthetic)

Four mice, two per group, sampled at `pre` and `end`; allele G at one site.
Baseline reads of G / total: m1 5/30, m2 0/30, m3 1/30, m4 0/2.

| mouse | pre status | end status | origin_in_own_mouse |
|---|---|---|---|
| m1 | present | present | standing_variation |
| m2 | absent | present | de_novo_candidate |
| m3 | below_detection | present | de_novo_candidate_below_detection_at_pre |
| m4 | not_covered | present | pre_not_covered |

Summary row: `origin_any_mouse = standing_variation` (m1 had it);
`n_pre_samples_allele_present / n_pre_samples_covered = 1 / 3` (m4 is not
covered and out of the denominator); `total_reads_pre = 92`,
`allele_reads_pre = 6` (m4's two reads **are** in the total).

A rising-table example (a different site from the one above): G and A tie, G
rose in fat by +0.667 and in control by +0.017, A fell in both.  The rising
table has two rows, both G, and none for A.  The fat row is
`standing_variation` (a fat or control mouse had G at `pre`).  The control row
is `allele_not_present_at_end`: its `end` check looks at control mice only, and
neither has G above the bar there (2/30 and 0/30 reads), although the summary,
which pools both groups' `end` samples, calls G standing.

## Running by hand

```bash
# One MAG's ANI table (QC files: one per timepoint)
alleleflux-pairwise-ani --mag MAG_A --cache_dir run/longitudinal/allele_freq_cache \
  --qc_files run/longitudinal/QC/QC_pre/MAG_A_QC.tsv run/longitudinal/QC/QC_end/MAG_A_QC.tsv \
  --fasta ref.fa --mag_mapping mapping.tsv --output_dir run/longitudinal/pairwise_ani \
  --pairs transitions --transitions pre:end

# Its strain verdicts, then the all-MAG classification
alleleflux-strain-turnover --mag MAG_A --pair_table run/longitudinal/pairwise_ani/MAG_A_pairwise_ani.tsv \
  --output_dir run/longitudinal/strain_turnover --transitions pre:end
alleleflux-replacement-classification --turnover_dir run/longitudinal/strain_turnover \
  --output_path run/longitudinal/strain_turnover/replacement_classification.tsv

# Baseline presence for one comparison and one test
alleleflux-baseline-presence --run_dir run/longitudinal --comparison pre_end-fat_control \
  --summary two_sample_paired --test_type two_sample_paired_tTest \
  --cache_dir run/longitudinal/allele_freq_cache --metadata metadata.tsv \
  --output_dir run/longitudinal/baseline_presence \
  --turnover_dir run/longitudinal/strain_turnover --cpus 8
```

Under the `cohort` scope the baseline-presence command reads, besides the
comparison's own samples, every sheet sample whose timepoint is in
`--witness_timepoints` and whose group is in `--witness_groups`.  Leave both
out, as above, and every sample of the sheet is read; on a dated sheet the
command then needs a cache file for every MAG it assesses at every timepoint
of those samples and stops on a missing one.  The workflow always passes both
lists, so it asks for exactly the caches it built.  Add
`--pre_treatment_timepoints pre` for the three-label split.

## In the workflow

```yaml
analysis:
  use_pairwise_ani: true
  pairwise_ani:
    pairs: transitions          # what strain turnover consumes
    store_snp_locations: none
  use_strain_turnover: true
  strain_turnover:
    min_compared: 0.1
    pop_threshold: 0.99999
    con_threshold: 0.999
  use_baseline_presence: true
  baseline_presence:
    summary: two_sample_paired
    test_type: two_sample_paired_tTest
    threshold_column: q_value
    threshold: 0.05
    baseline_scope: cohort          # needs collection_date in the metadata
    cohort_witnesses: sheet         # sheet | comparisons (cohort scope only)
    pre_treatment_timepoints: []    # e.g. ["pre"]; empty = no three-label split
    cluster_window_bp: 1000
    max_de_novo_neighbours: 0
    strain_gate: carrier            # carrier | group | both | none
    min_determined_carriers: 3
```

The seven de novo knobs are read once (`de_novo_rule_params` in `common.smk`)
so the baseline-presence and rising-origin rules cannot disagree; every gate
but `none` requires `use_strain_turnover` (the rising-origin job then waits
for the classification table and the per-MAG turnover tables and receives
`--classification` and `--turnover_dir`), and the run stops at parse time
otherwise.

Pairwise ANI, strain turnover and classification run for the **tested MAGs
only**: the union, over every configured comparison and enabled test, of the
MAGs the eligibility (and, when enabled, preprocessing) checkpoints admitted.
MAGs no test ran on never get an ANI job.  Strain turnover and classification
are local rules (seconds each); pairwise ANI and baseline presence are cluster
jobs whose cost is dominated by loading profiles.

## Choosing a threshold

Both verdicts are always reported so the choice is made in the analysis, not
baked into the files.  `strain_replacement` at popANI 0.99999 is the strict
same-strain line and flags many pairs in deeply sampled data;
`dominant_strain_change` at conANI 0.999 is far more conservative.  Compare
the two counts in `replacement_classification.tsv` before deciding which
column an enrichment filter reads.
