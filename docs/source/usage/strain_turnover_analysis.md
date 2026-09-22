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
between-mouse strain-sharing questions).

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
significant sites (`threshold_column` ≤ `threshold`), finds the allele(s) that
carried the minimum p-value (at a biallelic site both alleles tie, and both are
reported with `n_alleles_tied_at_min_p = 2`), then opens **every** sample's
profile for both groups and both timepoints and gives each (site, allele,
sample) one of four statuses:

| status | meaning |
|---|---|
| `present` | covered, the allele clears the bar and `min_freq` |
| `below_detection` | covered, some reads of the allele but under the bar |
| `absent` | covered, zero reads of the allele |
| `not_covered` | fewer than `min_cov` reads, or the MAG has no profile for the sample |

The presence rule is the same one pairwise ANI uses, so the two agree by
construction; the workflow passes the same `min_cov`, `min_freq` and `fdr`.

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
  `n_mice_de_novo_candidate_below_detection_at_{earlier}`.
* **Read counts, unfiltered:** `total_reads_{tp}`, `allele_reads_{tp}`,
  `allele_frequency_{tp}` for each timepoint, summed over **every** sample at
  that timepoint, thin and profile-less samples included.  That is what makes
  "never seen in N reads at baseline" a frequency bound of 1/N.  The two kinds
  can disagree on purpose: a `de_novo_candidate` site can have a few baseline
  reads of the allele sitting in samples under the detection bar.

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

## Running by hand

```bash
# One MAG's ANI table (QC files: one per timepoint combination)
alleleflux-pairwise-ani --mag MAG_A --profiles_dir run/longitudinal/profiles \
  --qc_files run/longitudinal/QC/QC_pre_end/MAG_A_QC.tsv \
  --fasta ref.fa --mag_mapping mapping.tsv --output_dir run/longitudinal/pairwise_ani \
  --pairs transitions --transitions pre:end --cpus 8

# Its strain verdicts, then the all-MAG classification
alleleflux-strain-turnover --mag MAG_A --pair_table run/longitudinal/pairwise_ani/MAG_A_pairwise_ani.tsv \
  --output_dir run/longitudinal/strain_turnover --transitions pre:end
alleleflux-replacement-classification --turnover_dir run/longitudinal/strain_turnover \
  --output_path run/longitudinal/strain_turnover/replacement_classification.tsv

# Baseline presence for one comparison and one test
alleleflux-baseline-presence --run_dir run/longitudinal --comparison pre_end-fat_control \
  --summary two_sample_paired --test_type two_sample_paired_tTest \
  --profiles_dir run/longitudinal/profiles --metadata metadata.tsv \
  --fasta ref.fa --mag_mapping mapping.tsv --output_dir run/longitudinal/baseline_presence \
  --turnover_dir run/longitudinal/strain_turnover --cpus 8
```

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
```

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
