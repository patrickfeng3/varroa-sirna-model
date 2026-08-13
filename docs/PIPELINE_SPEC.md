# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.3  
**Status:** Scientific specification before canonical implementation  
**Scope:** Viral small-RNA analysis through viral spatial/transitivity-consistency analysis  
**Host transitivity:** Excluded  
**vdCHIBIN ranking:** Excluded from this build

---

## 1. Purpose

This repository will contain the canonical, reproducible downstream analysis of *Varroa destructor* viral small-RNA sequencing data.

The expensive upstream work—read preprocessing, virus discovery, sample-specific consensus generation, strict mapping, and read-level feature extraction—has already been completed and audited. The existing legacy project is therefore treated as a **read-only validated data core**.

The canonical repository must regenerate the biological analyses from frozen inputs without recreating the ~60 GB upstream project.

Conceptual workflow:

```text
Validated legacy core
        ↓
00 Validate legacy core
        ↓
01 Fixed 23/24-nt populations
        ↓
02 Terminal nucleotide enrichment
        ↓
03 Official stepRNA
        ↓
04 Dicer evidence and Dicer-conditioned features
        ↓
05 Viral spatial/transitivity-consistency analysis
```

A key reproducibility principle is that Stage 05 has two named outputs:

- **historical_v1.4.1_replication** — reproduces the uploaded v1.4.1 algorithm exactly enough to match its archived results.
- **canonical_transitivity_analysis** — preserves the biological endpoints of v1.4.1 but improves cross-dataset inference by respecting sample-level clustering.

The historical result is a regression target, not an input to the canonical calculation.

---

## 2. Frozen legacy data core

The external data location is supplied locally through:

```text
config/paths.local.yaml
```

Current local example:

```text
/Users/patrickmod/Desktop/varroa_all_samples_pipeline_v1.0.0
```

This machine-specific path must not be committed to GitHub.

The legacy directory is **read-only**.

### Validated reusable layers

The completed audit supports reuse of:

- 21 corrected processed small-RNA FASTQs;
- 21 corrected preprocessing audit records;
- 273 virus-discovery mappings;
- approved virus metadata and selected sample-virus manifest;
- sample-specific final viral consensuses;
- sample-specific depth-masked background consensuses;
- 21 competitive exact SAM files;
- 21 competitive one-mismatch SAM files;
- 21 read-level feature tables;
- eligibility and mapping-summary tables.

The canonical pipeline reads these products but does not silently reuse old downstream 23/24, Dicer, or transitivity summaries.

---

## 3. General analysis principles

### 3.1 Analysis levels

A **sample** is one sequencing library/sample identifier such as an SRR accession.

A **sample-virus unit** is one virus analysed within one sample.

A **sample-virus-contig unit** is one mapped viral reference contig within one sample-virus unit. In the archived v1.4.1 transitivity result, each analysed sample-virus unit contributed one analysed contig, but the new implementation must not assume this is universally true.

Several virus observations from the same sample share the same library context and must not automatically be treated as fully independent biological replicates.

### 3.2 Weighting modes

Two modes are retained and never mixed:

**Abundance mode** — a sequence contributes according to its observed read abundance.

**Unique-sequence mode** — each distinct RNA sequence contributes total weight 1 within the explicitly defined analysis unit, regardless of how many QNAME/read rows carry that sequence.

### 3.3 Primary uncertainty framework

For cross-dataset canonical inference, the top-level resampling cluster is the **sample**.

Where a statistic is first computed per sample-virus-contig, the canonical sample-balanced summary is:

1. calculate the statistic for each eligible sample-virus-contig;
2. take the median across eligible virus-contigs within each sample;
3. take the median across sample-level medians;
4. bootstrap samples with replacement, retaining all observations belonging to a selected sample together.

Historical pair-balanced and virus-balanced summaries are retained as sensitivity/reproduction views.

### 3.4 Randomness

All bootstrap and permutation procedures must record:

- seed;
- number of requested replicates;
- number of valid replicates;
- exact aggregation rule.

### 3.5 Multiple testing

The hypothesis family must be defined before interpreting significance.

Historical v1.4.1 FDR adjustment is reproduced exactly in the historical output. The canonical output uses a broader predefined family described in Stage 05 so that analytical choices are not hidden in separate three-test families.

---

# 00 — Validate legacy core

## Purpose

Confirm that the frozen legacy project still contains all files and schemas required by the canonical downstream workflow.

This stage performs no remapping and no biological inference.

## Required inputs

At minimum:

```text
results/descriptive/eligibility.tsv
config/virus_catalog.tsv
tables/<sample>/<sample>.read_level_features.tsv.gz
alignments/<sample>.all_viruses.exact.sam
references/consensus/<sample>.<virus>.final.fa
references/consensus/<sample>.<virus>.final.background_masked.fa
```

plus the relevant sample/reference manifests and provenance records.

## Required checks

Confirm:

- all expected libraries are present;
- corrected preprocessing provenance is present;
- required read-level columns exist;
- exact SAM files are readable and structurally valid;
- required sample-specific consensuses exist;
- sample names and virus names agree across inputs;
- no downstream output path resolves inside the legacy core;
- frozen input files have recorded identity information such as size and, where practical, checksum.

## Outputs

```text
results/00_validation/
    legacy_core_validation.tsv
    legacy_core_validation.md
```

## Failure behaviour

A missing or inconsistent required dependency causes an explicit failure. There is no silent fallback to an old downstream result.

---

# 01 — Fixed 23-nt and 24-nt viral small-RNA populations

## Biological question

What are the abundance, strand distribution, and cross-sample behaviour of the 23-nt and 24-nt viral small-RNA populations?

## Inputs

Validated read-level feature tables, eligibility table, and virus metadata.

No remapping is performed.

## Primary inclusion criteria

Use reads satisfying the validated equivalents of:

```text
mapping_mode = exact
virus_assignment = assigned
primary_eligible = true
length = 23 or 24 nt
strand = sense or antisense
```

Cross-virus ambiguous reads are excluded from virus-specific primary summaries.

## Required dimensions

Retain:

- sample;
- analysis unit / virus;
- biological virus where distinct;
- viral polarity;
- strand;
- length;
- sequence;
- abundance.

## Required summaries

For 23 nt and 24 nt separately, calculate:

- total abundance;
- distinct sequence count;
- sense abundance/count;
- antisense abundance/count;
- sense fraction;
- antisense fraction;
- pair-level 23:24 relationships;
- sample-balanced and pair-balanced across-dataset summaries.

## Outputs

```text
results/01_viral_23_24/
    23_24_counts_by_pair.tsv
    23_24_fractions_by_pair.tsv
    23_24_strand_bias_by_pair.tsv
    23_24_across_samples.tsv
    23_24_across_pairs.tsv
    figures/
```

---

# 02 — Length-matched terminal nucleotide enrichment

## Biological question

Do observed 23-nt and 24-nt Varroa viral small RNAs contain terminal nucleotide patterns more or less often than expected from the viral sequence actually available for processing?

## Inputs

- exact eligible 23/24 read-level data;
- sample-specific depth-masked viral consensuses;
- eligibility table.

## Terminal positions

For each physical sequenced RNA in its own 5′→3′ orientation:

```text
5p1 = first nucleotide
5p2 = second nucleotide
3p2 = penultimate nucleotide
3p1 = final nucleotide
```

## Observed frequency

For each sample × virus × length × strand × weighting mode, calculate the observed A/C/G/U frequency at each terminal position.

## Expected frequency

Enumerate every fully supported window of the same length from the sample-specific depth-masked viral consensus.

- Sense expectation uses reference orientation.
- Antisense expectation uses reverse-complement orientation.
- Combined-strand expectation is weighted by the observed sense/antisense mixture for the corresponding unit; it is not forced to 50:50.

## Enrichment ratio

```text
enrichment_ratio = observed_fraction / expected_fraction
```

Interpretation:

```text
1   = observed as often as sequence availability predicts
>1  = enriched
<1  = depleted
```

If the expected fraction is zero, the ratio is `NA`; no arbitrary pseudocount is introduced.

## Across-dataset summaries

Retain both:

- historical/pair-level median enrichment across eligible sample-virus units;
- sample-balanced median of within-sample medians.

Confidence intervals for the canonical sample-balanced summary use sample-clustered bootstrap resampling.

The historical design-facing field `median_enrichment_ratio` remains available for regression comparison and must be clearly labelled if a new sample-balanced field is also exported.

## 23-versus-24 comparison

Use Spearman rank correlation to compare matched 23- and 24-nt enrichment landscapes, including at least overall and antisense-specific comparisons.

## Outputs

```text
results/02_terminal_enrichment/
    fixed_length_positional_nucleotides_by_pair.tsv
    fixed_length_positional_nucleotides_across_samples.tsv
    fixed_length_positional_nucleotides_across_pairs.tsv
    ALL_VIRUSES_23nt_positional_nucleotide_ratios.tsv
    ALL_VIRUSES_24nt_positional_nucleotide_ratios.tsv
    ALL_VIRUSES_23nt_strand_specific_positional_nucleotide_ratios.tsv
    ALL_VIRUSES_24nt_strand_specific_positional_nucleotide_ratios.tsv
    23_vs_24_enrichment_correlations.tsv
    figures/
```

## Interpretation limit

This is an empirical sequence-enrichment measurement. It does not by itself identify whether the preference arose from Dicer cleavage, Argonaute loading, strand selection, RNA stability, library effects, or another process.

---

# 03 — Official stepRNA Dicer-overhang analysis

## Biological question

Do the 23-nt and/or 24-nt populations show complementary guide/passenger duplex-end geometry consistent with Dicer-like processing?

## Primary software

Use the official published **stepRNA** implementation as the primary duplex-overhang method.

The older project-native reconstruction is retained only for diagnostic comparison and historical validation.

## Why the full overhang spectrum is analysed

The analysis does not define Dicer as “2-nt overhang or nothing.” stepRNA reports the full signed 5′ and 3′ overhang/underhang distributions and passenger lengths. The previously observed Varroa +2 5′ / −2 3′ joint geometry is a pre-specified feature of interest within that broader distribution.

Official stepRNA sign convention must be preserved:

```text
negative distance = reference overhang
positive distance = reference underhang
0                 = blunt end
```

## File A — reference populations

Run separately for:

```text
23-nt sense
23-nt antisense
24-nt sense
24-nt antisense
```

within each eligible sample-virus unit.

## File B — potential passengers

Canonical Varroa project setting:

```text
15–30 nt
```

from the same sample and same virus, restricted to the opposite mapped strand relative to File A.

A pre-specified sensitivity run uses:

```text
18–28 nt
```

because passenger-range choice is an analysis parameter rather than a biological constant.

## Alignment behaviour

The canonical run uses official stepRNA exact complementary alignment behaviour. Mismatch-tolerant analyses, if ever performed, are exploratory and separately labelled.

## Collapsed and non-collapsed modes

Retain both official views:

- unique/collapsed reference sequences;
- non-collapsed/expression-weighted reference sequences.

## Outputs retained

Retain official stepRNA outputs including:

- 5′ overhang/underhang counts;
- 3′ overhang/underhang counts;
- unique-overhang counts;
- passenger number;
- passenger length;
- overhang type;
- relevant classified alignment files;
- official log-ratio/log-odds-style enrichment output as produced by the installed stepRNA version;
- official Wald Z-scores as produced by stepRNA.

## Outputs

```text
results/03_steprna/
    inputs/
    raw_outputs/
    summaries/
    sensitivity/
```

## Interpretation limit

A reproducibly enriched duplex-end geometry is evidence consistent with Dicer/Dicer-like cleavage. It does not directly observe intact duplexes in vivo, identify a specific Dicer paralogue, or prove that every RNA of that length was Dicer-generated.

---

# 04 — Dicer evidence aggregation and Dicer-conditioned sequence features

This stage has two distinct purposes and keeps them separate.

## 04A — Population-level Dicer evidence

### Question

How reproducible is the stepRNA Dicer-like signal across the Varroa dataset, and how do 23- and 24-nt populations compare?

### Primary input

Official stepRNA outputs from Stage 03.

### Required summaries

For each sample-virus unit and File-A class, retain:

- complete 5′ distance distribution;
- complete 3′ distance distribution;
- official stepRNA enrichment/Z-score outputs;
- passenger recovery fraction;
- passenger-length distribution;
- support for the pre-specified Varroa +2 5′ / −2 3′ joint geometry.

### Passenger recovery versus geometry

These are reported separately.

A low fraction of references with any recoverable passenger does not equal absence of Dicer processing. Among references with a recovered passenger, the geometry distribution is reported independently.

### Project-specific secondary statistic

The historical project statistic `Δ_Dicer` may be reproduced as a secondary validation statistic:

```text
Δ_Dicer = support at a pre-specified Dicer-compatible distance
          - mean support across a pre-specified comparison-distance set
```

The comparison-distance set must be defined in configuration before examining results.

`Δ_Dicer` is not an official stepRNA statistic.

### Aggregation

Canonical confidence intervals use sample-aware aggregation/resampling. Pair-balanced and virus-balanced views are retained as sensitivity summaries.

## 04B — Dicer-conditioned sequence features

### Question

Do RNAs with strong pre-specified Dicer-like geometry have terminal sequence features that differ from the general 23/24 population?

### Primary Dicer-supported subset

Use the pre-specified Varroa joint geometry rather than selecting whichever overhang distance happens to maximize enrichment in the same dataset.

Other reproducibly enriched geometries may be explored but are labelled exploratory until independently justified.

### Initial features

Pre-specify:

```text
5p1
5p2
3p2
3p1
guide length
reference strand
passenger length
```

### Absolute Dicer-conditioned enrichment

```text
E_Dicer_absolute = frequency among Dicer-supported RNAs
                   / matched viral-sequence expected frequency
```

### Dicer-specific contrast

```text
dicer_specific_log2_contrast
    = log2(E_Dicer_absolute / E_all_observed)
```

where `E_all_observed` is the matching general enrichment from Stage 02.

This asks whether Dicer-supported RNAs contain information beyond the general sequence preference already seen among all siRNAs.

Undefined or zero cases that make the logarithm invalid are reported as `NA`; arbitrary pseudocounts are not introduced merely to force finite values.

### Redundancy test

Quantify the relationship between general terminal enrichment and Dicer-conditioned enrichment. A strongly redundant Dicer-derived signal should not later be treated as an independent candidate-scoring dimension without justification.

## Outputs

```text
results/04_dicer_features/
    dicer_population_summary_by_pair.tsv
    dicer_population_summary_by_sample.tsv
    dicer_23_vs_24.tsv
    dicer_passenger_length_summary.tsv
    dicer_conditioned_terminal_features.tsv
    dicer_specific_contrasts.tsv
    dicer_vs_general_enrichment.tsv
    dicer_sensitivity_summary.tsv
    figures/
```

---

# 05 — Viral spatial/transitivity-consistency analysis

## 05.1 Biological question

Is the 24-nt antisense population spatially related to primary-like 23-nt processing in a pattern **consistent with** secondary/transitive amplification?

This is an observational analysis of natural viral infections. Viral replication can itself create spatially structured complementary RNA, so the result cannot by itself prove RdRP-mediated transitivity.

## 05.2 Audited historical reference implementation

The exact uploaded v1.4.1 reference package is:

```text
Varroa_vsiRNA_v1.4.1_strengthened_transitivity/
    analysis_tools/analyse_strengthened_transitivity.py
    tests/test_strengthened_transitivity.py
    docs/ACADEMIC_BASIS_AND_INTERPRETATION.md
    run_strengthened_transitivity.py
```

The corresponding archived result package is:

```text
Varroa_v1.4.1_strengthened_transitivity_results/
```

The canonical repository must preserve these packages outside the executable workflow or document their checksums/provenance so that the historical algorithm is auditable.

## 05.3 Historical v1.4.1 default parameters

Exact defaults recovered from code:

```text
bin_size_nt              = 10
max_crosscorr_lag_nt     = 500
windows_nt               = [100, 250, 500]
anchor_percentile        = 90.0
anchor_min_separation_nt = 50
minimum_anchors          = 3
permutations             = 5000
bootstrap_replicates     = 5000
random_seed              = 20260810
```

These values are historical analysis parameters, not biological constants.

## 05.4 Eligible viral units

Use sample-virus units that are:

- `primary_eligible = true`; and
- classified as positive-sense RNA viruses by the validated virus catalogue.

The historical result contained 14 samples, 19 positive-sense sample-virus units, and three analysis units/viruses contributing to the final transitivity summaries. The new workflow must derive these counts from the frozen inputs rather than hard-code them.

## 05.5 Read selection

Within eligible units, select reads satisfying the validated equivalents of:

```text
mapping_mode = exact
virus_assignment = assigned
length = 23 or 24 nt
strand = sense or antisense
```

## 05.6 Coordinate recovery from exact SAM files

Coordinates are recovered from the frozen exact competitive SAM files; reads are not remapped.

Historical parser behaviour:

- unmapped records are excluded;
- supplementary records are excluded;
- all retained exact compatible loci are available for fractional multimapper handling;
- strand mismatches are counted diagnostically.

Canonical behaviour strengthens QC:

- a metadata conflict causes explicit exclusion and reporting;
- any non-zero strand mismatch causes Stage 05 QC failure unless a documented cause is resolved;
- duplicate alignment records for the same physical locus are deduplicated before weighting.

## 05.7 Positional coordinate used for each small RNA

The historical analysis uses the alignment midpoint:

```text
midpoint_nt = alignment_start_0based + (read_length - 1) / 2
```

The midpoint is assigned to a 10-nt bin:

```text
bin_index = floor(midpoint_nt / 10)
```

This midpoint/bin convention is retained in both historical replication and canonical analysis.

## 05.8 Abundance track construction

For one QNAME/read with observed abundance `a` and `k` unique exact compatible loci:

```text
weight per locus = a / k
```

Thus the QNAME contributes total abundance `a` across all of its exact compatible positions rather than `a` at every locus.

Canonical code requires the SAM strand to match the read metadata strand before a locus is accepted.

## 05.9 True unique-sequence track construction

Define sequence identity by:

```text
virus × strand × length × sequence
```

All QNAMEs carrying that same physical RNA sequence are collapsed.

If the sequence has `k` unique exact compatible loci across the analysed reference contigs:

```text
weight per locus = 1 / k
```

Therefore each distinct sequence contributes **total weight exactly 1**, regardless of read abundance or number of QNAME rows.

This is the corrected v1.4.1 unique-sequence definition and must be explicitly tested.

## 05.10 Position-wise tracks

For each sample-virus-contig and weighting mode, construct 10-nt binned tracks:

```text
23 sense
23 antisense
24 sense
24 antisense
```

## 05.11 Primary-like 23-nt anchor scores

Two predefined anchor signals are analysed separately.

### `balanced23_anchor_score`

```text
balanced23 = sqrt(23S × 23AS)
```

A high value requires signal from both strands.

### `combined23_anchor_score`

```text
combined23 = 23S + 23AS
```

This measures total local 23-nt activity without requiring strand balance.

Neither score proves that a locus is biologically primary; they define **primary-like spatial anchors** for this analysis.

## 05.12 Exact historical anchor selection

For each binned anchor-score track:

1. keep bins with score > 0;
2. if fewer than `minimum_anchors` non-zero bins exist, the unit is ineligible for that anchor analysis;
3. calculate the 90th percentile **only across non-zero bins**;
4. candidate anchors are non-zero bins with score at or above that threshold;
5. rank candidates from strongest score to weakest;
6. greedily retain a candidate only if it is at least 50 nt from all already chosen anchors;
7. after selection, require at least three anchors.

The historical code has no separate hotspot-merging step.

Canonical code must make tie-breaking deterministic: equal-score candidates are ordered by genomic coordinate after score ranking.

## 05.13 Exact upstream/downstream window calculation

For each selected anchor and window `W ∈ {100, 250, 500}` nt:

- convert `W` to 10-nt bins;
- exclude the anchor bin itself;
- collect available downstream bins and upstream bins separately;
- truncate naturally at reference boundaries;
- calculate mean signal per valid bin.

Historical pooled mean:

```text
mean_downstream = sum(signal over all anchor-specific downstream bins)
                  / total number of valid downstream bin observations

mean_upstream   = sum(signal over all anchor-specific upstream bins)
                  / total number of valid upstream bin observations
```

If neighbourhoods around different anchors overlap, a genomic bin can contribute once for each anchor neighbourhood that contains it. Therefore the historical statistic is an **anchor-window pooled mean density**, not the mean over unique genomic territory.

This property is retained for historical replication and explicitly documented in canonical outputs.

## 05.14 Normalized 24-nt directionality

For a strand-specific 24-nt track:

```text
D = (mean_downstream - mean_upstream)
    / (mean_downstream + mean_upstream)
```

If the denominator is zero or the required means are not finite, `D = NA`.

Calculate:

```text
D_24AS = normalized directionality of 24-nt antisense
D_24S  = normalized directionality of 24-nt sense
```

Primary directionality contrast:

```text
antisense_specific_directionality = D_24AS - D_24S
```

A positive value means 24-AS is more downstream-biased than the 24-S control track.

## 05.15 Antisense 23→24 composition endpoint

Using the anchor-window pooled mean densities:

```text
F24_AS_down
    = mean24AS_down / (mean23AS_down + mean24AS_down)

F24_AS_up
    = mean24AS_up / (mean23AS_up + mean24AS_up)
```

If a denominator is zero, the corresponding fraction is `NA`.

The composition shift is:

```text
delta_F24_AS = F24_AS_down - F24_AS_up
```

A positive value means that, among antisense 23+24-nt signal, the downstream neighbourhood is relatively more 24-nt dominated than the upstream neighbourhood.

This is a composition statistic. It does not necessarily imply an increase in total small-RNA abundance.

## 05.16 Historical circular-shift null

For each sample-virus-contig, weighting mode, anchor definition, and window:

- keep the 23-nt tracks and anchors fixed;
- shift the entire 24-AS track and 24-S track by the **same randomly selected circular shift**;
- recompute `antisense_specific_directionality` and `delta_F24_AS`;
- repeat 5000 times by default.

Using the same shift for 24-AS and 24-S preserves their mutual spatial relationship while disrupting their registration relative to the 23-nt anchors.

Allowed historical shifts are non-zero circular shifts satisfying:

```text
min(shift_bins, n_bins - shift_bins) > max_window_bins
```

where `max_window_bins = 500 / 10 = 50` under default settings.

Thus shifts within 500 nt of zero in either circular direction are excluded when possible.

If no such shifts exist for a short reference, the historical function falls back to **all non-zero circular shifts**. The canonical replication must report whenever this fallback occurs.

For each permutation index, different contigs may receive different randomly drawn allowed shifts; higher-level null statistics are then calculated by aggregating the same permutation index across contigs.

## 05.17 Historical empirical permutation P-value

For the one-sided hypothesis that the observed statistic is greater than expected under the shift null:

```text
p = (1 + number of null statistics >= observed statistic)
    / (1 + number of valid null statistics)
```

The tail direction is fixed in advance.

## 05.18 Historical pair-balanced aggregation

For a fixed:

```text
weighting × anchor definition × window
```

v1.4.1 takes the median of the finite sample-virus-contig statistics.

This is called **pair-balanced** in the archived package because each eligible row contributes comparably rather than being weighted by sequencing depth.

The historical bootstrap independently resamples those rows with replacement and recalculates the median.

This historical bootstrap does **not** preserve sample clustering and is retained only for exact reproduction/sensitivity.

## 05.19 Historical virus-balanced aggregation

For each biological analysis unit/virus:

1. calculate the median across its eligible sample-virus-contigs;
2. calculate the median across virus medians.

The historical hierarchical bootstrap resamples viruses and then observations within selected viruses.

This is a robustness view intended to reduce domination by a virus represented in many samples. With few viruses, its uncertainty must be interpreted cautiously.

## 05.20 Canonical sample-balanced aggregation

This is the primary cross-dataset aggregation for the new canonical analysis.

For every fixed weighting × anchor × window × endpoint:

1. calculate the endpoint for each eligible sample-virus-contig;
2. within each sample, take the median across its finite virus-contig endpoint values;
3. take the median across sample-level medians.

### Canonical confidence interval

Bootstrap the sample IDs with replacement. When a sample is selected, retain all of its eligible virus-contig observations together, recompute its within-sample median, and then recompute the across-sample median.

### Canonical permutation null

For each permutation index:

1. apply the v1.4.1 same-shift 24S/24AS null independently within each contig;
2. calculate the endpoint per contig;
3. collapse contig endpoints to a median within each sample;
4. take the median across samples.

The empirical P-value uses the same `(b+1)/(m+1)` one-sided formula.

This retains the v1.4.1 spatial null while preventing samples containing multiple eligible viruses from acting like multiple independent samples at the final aggregation level.

## 05.21 Multiple testing

### Historical replication

Reproduce v1.4.1 exactly: BH adjustment across the three windows `100/250/500 nt`, separately for each:

```text
weighting × anchor definition × endpoint × aggregation type
```

Thus each historical BH family contains three P-values.

### Canonical inference

For the primary sample-balanced analysis, define one confirmatory family per biological endpoint across all predefined combinations:

```text
2 weighting modes × 2 anchor definitions × 3 windows = 12 tests
```

The two endpoint families are:

- `delta_F24_AS`;
- `antisense_specific_directionality`.

Apply BH separately to those two 12-test families.

Pair-balanced, virus-balanced, leave-one-virus-out, cross-correlation, and any future parameter-sensitivity analyses are robustness/descriptive outputs rather than extra routes to a primary claim.

## 05.22 Leave-one-virus-out robustness

For each weighting × anchor × window, omit each virus in turn and recalculate the pair-level observed medians.

The canonical pipeline may additionally report sample-balanced leave-one-virus-out summaries, but these are sensitivity analyses rather than a new confirmatory family.

## 05.23 Descriptive cross-correlation

Cross-correlation is retained as spatial description, not primary transitivity evidence.

Historical implementation:

- transform anchor and 24-nt tracks with `log1p`;
- evaluate lags from −500 to +500 nt in 10-nt steps;
- positive lag means the 24-nt target track is displaced downstream relative to the 23-nt anchor track;
- calculate Pearson correlation at each lag where at least 10 overlapping bins exist and both vectors have non-zero variance.

A lag-asymmetry summary compares mean positive-lag correlation with mean negative-lag correlation for each 24-nt strand, then compares antisense minus sense.

Zero-lag co-localisation alone is not interpreted as transitivity.

## 05.24 Historical regression checkpoints

The historical replication should match the archived v1.4.1 results to numerical tolerance.

Important diagnostic checkpoints include:

```text
14 samples in coordinate diagnostics
19 eligible positive-sense sample-virus units
0 metadata conflicts
0 strand mismatches
```

For `unique_sequence × balanced23`, archived pair-balanced `delta_F24_AS` values are approximately:

```text
100 nt  = -0.000113
250 nt  = +0.002682
500 nt  = +0.003662
```

with archived pair-level BH-adjusted values approximately:

```text
100 nt  = 0.471706
250 nt  = 0.018596
500 nt  = 0.000600
```

The corresponding archived `antisense_specific_directionality` results are not significant and are near zero/negative rather than showing a convincing absolute downstream 24-AS wave.

These numbers are **regression checkpoints only**. They are never loaded as inputs to the new analysis.

## 05.25 Interpretation rule

The viral analysis may support statements such as:

- 23- and 24-nt spatial association;
- downstream/upstream asymmetry;
- an antisense-specific directional effect if present;
- a downstream shift in antisense 23/24 length composition;
- consistency with secondary/amplification-associated biology.

It must not by itself establish:

- that every 23-mer is primary;
- that every 24-mer is secondary;
- that RdRP directly synthesizes 24-mers;
- that 24-mers are Dicer-independent;
- a universal propagation distance;
- a specific Dicer/Ago/RdRP paralogue;
- host-mRNA transitivity.

## 05.26 Required outputs

```text
results/05_viral_transitivity/
    coordinate_qc.tsv
    eligible_positive_sense_units.tsv
    historical_v1.4.1_replication/
        transitivity_by_pair.tsv
        pair_balanced_results.tsv
        virus_balanced_results.tsv
        leave_one_virus_out.tsv
        cross_correlation.tsv
        regression_check.tsv
    canonical_transitivity_analysis/
        transitivity_by_pair.tsv
        transitivity_by_sample.tsv
        sample_balanced_results.tsv
        pair_balanced_sensitivity.tsv
        virus_balanced_sensitivity.tsv
        leave_one_virus_out.tsv
        multiple_testing_summary.tsv
        final_transitivity_summary.tsv
    figures/
```

---

## 6. Reproducibility requirements

All new outputs are written under `results/` in the canonical repository.

Configuration must record at minimum:

```text
target_lengths = [23, 24]
steprna_passenger_range = [15, 30]
steprna_sensitivity_range = [18, 28]
transitivity_bin_size_nt = 10
transitivity_windows_nt = [100, 250, 500]
transitivity_anchor_percentile = 90
transitivity_anchor_min_separation_nt = 50
transitivity_min_anchors = 3
transitivity_permutations = 5000
bootstrap_replicates = 5000
random_seed = 20260810
```

Each run records:

- pipeline Git commit;
- configuration file;
- software versions;
- legacy-core path;
- input identity/checksums where practical;
- run date;
- random seed.

Relevant software versions include Python, Snakemake, stepRNA, Bowtie2 used by stepRNA, pysam, NumPy, pandas, SciPy, and plotting/statistical packages actually used.

---

## 7. Required deterministic tests

### Stage 02

- 5p1/5p2/3p2/3p1 extraction;
- reverse-complement antisense expectation;
- strand-weighted combined expectation;
- zero expected frequency → `NA`.

### Stage 03–04

- correct File-A class;
- opposite-strand File-B selection;
- passenger-length filters;
- official stepRNA sign convention;
- parsing of official stepRNA outputs;
- separate passenger recovery and geometry denominators.

### Stage 05 coordinate and weighting logic

- alignment midpoint calculation;
- 10-nt bin assignment;
- QNAME abundance split across `k` loci sums to original abundance;
- duplicate QNAME alignments do not double-count a locus;
- duplicate physical sequences collapse correctly;
- one unique sequence across `k` loci sums to total weight 1;
- strand mismatch triggers canonical QC failure.

### Stage 05 anchor logic

- percentile is calculated over non-zero bins only;
- exact 90th-percentile candidate rule;
- deterministic tie-breaking;
- 50-nt separation rule;
- minimum three-anchor rule.

### Stage 05 window logic

- anchor bin excluded;
- downstream/upstream sign convention;
- boundary truncation;
- denominator equals number of valid anchor-bin observations;
- overlapping anchor windows reproduce historical repeated-bin weighting.

### Stage 05 endpoints

- exact normalized `D` formula;
- exact `D_24AS - D_24S` formula;
- exact `F24_AS` formula;
- zero denominator → `NA`;
- positive synthetic downstream 24 composition gives positive `delta_F24_AS`.

### Stage 05 permutation and aggregation

- same shift applied to 24-AS and 24-S within a contig;
- allowed-shift exclusion rule;
- short-reference fallback is reported;
- `(b+1)/(m+1)` P-value formula;
- historical pair median;
- historical virus median-of-medians;
- canonical sample median-of-medians;
- sample-cluster bootstrap keeps observations from a selected sample together;
- historical three-window BH family;
- canonical 12-test BH family;
- fixed seed reproduces identical outputs.

---

## 8. Explicitly superseded analyses

Older downstream results remain historical references only, including:

- legacy fixed-length 23/24 summaries;
- legacy custom Dicer-overhang summaries;
- v1.4.0 transitivity implementation;
- any v1.4.1 result loaded directly instead of recomputed from frozen inputs.

They may be used for numerical regression checks, never as silent analytical inputs.

---

## 9. Explicitly outside current scope

This build excludes:

- CHH host analysis;
- Pero host analysis;
- Ago2 host analysis;
- host transitivity;
- inferred host-trigger analysis;
- ViennaRNA vdCHIBIN accessibility scoring;
- vdCHIBIN thermodynamic asymmetry scoring;
- vdCHIBIN final ranking;
- construct architecture comparison;
- Nectar Designer integration.

---

## 10. Definition of success

The first canonical viral pipeline is complete when:

1. a fresh clone can point to the frozen validated core;
2. Stage 00 validates without modifying the core;
3. Stages 01–04 regenerate 23/24, enrichment, and official stepRNA/Dicer outputs from frozen inputs;
4. official stepRNA is the primary Dicer-overhang method;
5. sample clustering is respected in canonical cross-dataset uncertainty estimates;
6. historical Stage 05 reproduces the uploaded v1.4.1 results within numerical tolerance;
7. canonical Stage 05 recomputes the same biological endpoints with sample-balanced inference;
8. all major choices are configuration-controlled;
9. all critical calculations have deterministic tests;
10. all metrics are defined in `docs/METRIC_DICTIONARY.md`;
11. one Snakemake entry point regenerates the entire downstream viral analysis;
12. no manual movement or copying of intermediate outputs is required;
13. any deviation from historical results is surfaced explicitly rather than hidden.

---

## 11. Methodological references

The main methodological basis includes:

- Murcott B, Pawluk RJ, Protasio AV, Akinmusola RY, Lastik D, Hunt VL. 2022. *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. Frontiers in Bioinformatics 2:994871. DOI: 10.3389/fbinf.2022.994871.
- Benjamini Y, Hochberg Y. 1995. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. Journal of the Royal Statistical Society Series B 57:289–300.
- Phipson B, Smyth GK. 2010. *Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn*. Statistical Applications in Genetics and Molecular Biology 9:Article 39.
- Saravanan V, Berman GJ, Sober SJ. 2020. *Application of the hierarchical bootstrap to multi-level data in neuroscience*. Used here as general methodological support for respecting nested/clustered observations; the biological application in this project is different.
- Damayo J et al. 2026 preprint, *Primary and secondary antiviral RNAi responses throughout Varroa destructor life stages reveal the vertical transmission of viruses*. This is biological motivation for the historical 23→24 hypothesis. Mechanistic claims must remain limited to what the available data directly support.

---

## 12. Canonical workflow summary

```text
READ-ONLY VALIDATED LEGACY CORE
              │
              ▼
     00_validate_legacy_core
              │
              ▼
       01_viral_23_24
              │
              ▼
    02_terminal_enrichment
              │
              ▼
       03_official_steprna
              │
              ▼
       04_dicer_features
              │
              ▼
05_viral_transitivity_consistency
      ├── historical_v1.4.1_replication
      └── canonical_transitivity_analysis
```

This is the complete scope of the first canonical downstream Varroa viral small-RNA pipeline.
