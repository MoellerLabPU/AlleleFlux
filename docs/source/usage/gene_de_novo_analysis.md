# Gene-Level De Novo Accumulation

Which genes collect new alleles faster than the rest of their genome, and faster in one group than in the other, with the cage as the replicate. Switched on with `analysis.use_gene_de_novo: true` (longitudinal data only); one job per comparison, run by `alleleflux-gene-de-novo` right after the allele-frequency step. It builds on the site-level de novo rules of the [Strain Turnover and Baseline Presence guide](strain_turnover_analysis.md), so read that first for the presence rule, the baseline scopes and the strain gate.

## The question in plain words

Every cage keeps a diary of the new letters that appeared in a genome between the two timepoints of a comparison. Each entry is filed under the gene it landed in. Two questions are then asked of every gene:

* **Test A, within a group:** did the gene get more entries than its share of the readable genome? If a cage saw 20 new alleles across its 1,000,000 readable bases, a gene that is 1 percent of those bases deserves 0.2 of them. Every cage of the group is asked; the answer is discounted when the cages disagree with one another.
* **Test B, between the groups:** did the gene collect entries per readable base faster in one group than in the other? Treatment cages caught 4 in 400 metres of river, control caught 0 in 400 metres: could four-to-nothing happen by luck if both rivers were equally stocked, and did the four come from four cages or from one?

The cage is the replicate throughout: mice in a cage share a microbiome, so two mice of one cage carrying the same new allele are one observation, not two.

## What counts as an event

1. **Candidates** are the letters whose cage-average share went up between the two timepoints, read from the pipeline's mean-change table (`allele_freq.py`, mice sampled at both timepoints, after the zero-difference filter). Nothing else can become an event.
2. Each candidate is judged **per cage** with the site-level rules, reused as they are: a later-timepoint mouse of the cage must carry the letter by the presence rule (otherwise `not_present_later`); the letter must be absent from every covered baseline sample under `baseline_scope` (`standing_variation` when present, `de_novo_below_detection` when a whisper of it was there, `baseline_not_covered` when no baseline sample was deep enough); de novo alleles within `cluster_window_bp` of one another **in the same cage** form a block (`block_event`: one imported or recombined chunk, never independent mutations; the window never spans two cages); under `strain_gate: carrier` the cage's carriers are looked up in the strain-turnover table (`new_strain` when more than half of them changed strain, `strain_undetermined` when fewer than `min_determined_carriers` have a verdict).
3. **Counted events** are `de_novo` and `de_novo_below_detection`, plus `block_event` rows when `block_rule: once`.

## The denominator: readable bases, never gene length

A gene read over 300 of its 1,200 bases was searched over 300. A base is readable in a cage when some mouse of the cage that was sampled at both timepoints has at least `min_cov` reads there at **both** timepoints (the same `min_cov` as the presence rule). The per-gene and genome-wide readable base counts are the denominators of both tests. A cage made only of mice sampled once was never searched and contributes nothing, rather than "zero mutations in N bases". `min_gene_covered_fraction` (default 0.0) lets a thin reading of a gene keep that cage out of the gene's test; the covered fraction is reported either way.

## The two tests

Both are likelihood-ratio tests of Poisson counts with an exposure offset, in closed form, with a quasi-Poisson dispersion `phi = max(1, Pearson X² / df_resid)` that discounts a result when the cages disagree. They differ in the reference curve.

**Test A** (`gene_vs_genome`), once per group. Each cage contributes the events observed in the gene (`o`) and the events expected if the gene mutated at the cage's genome-wide rate (`e` = the cage's events × the gene's share of the cage's readable genome). Rate ratio `Σo / Σe`; `LRT = 2[Σo ln ρ − (ρ − 1) Σe]`; `p = P[χ²₁ > LRT / phi]`.

**Test B** (`group_vs_group`), once per gene. One event rate per group (events over readable bases, pooled over the group's cages) against a single pooled rate; `p = P[F(1, n_cages − 2) > LRT / phi]`. The F reference charges for the dispersion being estimated from few cages: in simulations with heterogeneous cages the chi-square reference raised false alarms about 12 percent of the time at 4 vs 4 cages, F held 5 to 6 percent; for Test A the F reference was far too conservative, so it keeps chi-square.

A gene is tested in a group only when at least `quality_control.min_sample_num` of that group's cages entered it (Test B: both groups); below that the note reads `below_cage_floor`. Benjamini-Hochberg runs within each of the three families (Test A per group, Test B) over the tested genes only, per MAG when `statistics.fdr_group_by_mag_id` is true. A gene with no events in a family is `no_events`, never p = 1.

**Both tests are also reported without the discount** (`p_value_raw_*`, with their own q-values): the same statistic looked up on the chi-square curve with the dispersion fixed at 1, the textbook Poisson test that treats every cage as an independent draw. Read the pair together: when the cages agree the two are close (the planted case: 0.0185 raw, 0.057 discounted, the gap being only the F curve's charge for estimating the dispersion); when one cage carries the events they split wide (ten events in one cage and none in three: 0.0002 raw, 0.147 discounted). The raw p uses the chi-square and never F, because no dispersion is estimated for it; in a simulation with no cage-to-cage noise the raw chi-square rejected 5 to 8 percent of null genes at 0.05 while an F lookup of the same statistic rejected 1 to 4 percent, and with strong cage noise the raw p rejected 18 to 42 percent of null genes on either curve, which is what the discount is for.

With 4 cages per group a perfect four-to-nothing split in Test B reads p = 0.057: suggestive, not significant. At 8 cages per group the same pattern is p = 0.005. The test has the power it needs at the floor, not below it.

## Worked cage

Cage A (treatment) has mice m1 and m2 at both timepoints. Position 150 of gene g1: every baseline sample reads 40 A and no T; m1's later sample reads 20 A and 20 T. The mean-change table says T rose in cage A, so (cage A, 150, T) is a candidate. m1 carries T later; no covered baseline sample has it; no other de novo allele within 1,000 bp in cage A; m1 kept its strain. One `de_novo` event, filed under g1. Cage A read all 300 bases of g1 and 1,000 bases of the genome, and caught one event in all: g1's fair share is 1 × 300 / 1000 = 0.3, observed 1. With one such event in each of two treatment cages and none in two control cages, Test B reads 4 cages, `df_resid_between` 2, rate ratio infinite, p = 0.24: two cages per group cannot say much, which is the point of the floor.

## Outputs

See [Outputs](../reference/outputs.md#gene-level-de-novo-accumulation): one tests table per comparison (one row per MAG × gene) and, per MAG, the events, the gene × cage table and the cage totals.

## Running by hand

```bash
alleleflux-gene-de-novo \
  --comparison pre_post-treatment_control \
  --cache_dir output/longitudinal/allele_freq_cache \
  --allele_analysis_dir output/longitudinal/allele_analysis/allele_analysis_pre_post-treatment_control \
  --metadata metadata/sample_metadata.tsv \
  --prodigal_fasta reference/prodigal_genes.fna \
  --fasta reference/combined_mags.fasta \
  --mag_mapping reference/mag_mapping.tsv \
  --output_dir output/longitudinal/gene_de_novo/pre_post-treatment_control \
  --strain_gate carrier --turnover_dir output/longitudinal/strain_turnover \
  --min_cages 4 --cpus 8
```

## In the workflow

```yaml
analysis:
  use_gene_de_novo: true
  gene_de_novo:
    block_rule: "exclude"
    min_gene_covered_fraction: 0.0
    min_determined_carriers: 1
  baseline_presence:          # the shared de novo rules
    baseline_scope: "cohort"
    cluster_window_bp: 1000
    max_de_novo_neighbours: 0
    strain_gate: "carrier"    # needs use_strain_turnover (and use_pairwise_ani); or none
```

The rule reads the comparison's mean-change tables and cache files, the sample sheet, the gene file, the reference and the MAG mapping, plus the strain-turnover tables under the carrier gate. `strain_gate: group` or `both` stops the run at parse time: a cage has no across-cage vote. Under `allele_analysis_only` no gene de novo job is made, like every other downstream analysis. In a permuted (null) run strain turnover is switched off, so the carrier gate cannot be used there: set `strain_gate: none` or `use_gene_de_novo: false`, the run stops at parse time otherwise and says so. The synonymous / nonsynonymous split of the events is deferred: computed against the assembly's codon it would be wrong at every fixed difference between the population and an unrelated reference genome, so it will be added against the baseline codon, the way the dN/dS step reconstructs its ancestral sequence.
