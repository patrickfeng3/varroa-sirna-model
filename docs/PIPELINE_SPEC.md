# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.16  
**Status:** Stages 00–06 implemented and validated; Stage 07 single-nucleotide empirical analysis implemented and validated, with the v0.16 fixed-width regional-GC extension scientifically specified pending implementation; Stage 08 candidate biophysics pre-specified and deferred until Stage 07 completion  
**Scope:** Viral small-RNA analysis, generic transcript candidate preparation, Stage 07 single-nucleotide and fixed-width regional-GC empirical guide-sequence association analysis, and pre-specified Stage 08 candidate biophysics  
**Host transitivity:** Excluded  
**Candidate ranking:** Excluded through Stage 08

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
01 Viral length landscape and fixed 23/24-nt populations
        ↓
02 Terminal nucleotide enrichment
        ↓
03 Official stepRNA
        ↓
04 Duplex-geometry evidence and geometry-conditioned features
        ↓
05 Viral spatial/transitivity-consistency analysis
        ↓
06 Generic transcript candidate enumeration

07 Varroa empirical guide-sequence association landscape
   (viral small-RNA dataset; independent computational branch)
        ↓
08 Generic candidate biophysics:
   target accessibility + duplex-end asymmetry
```

A key reproducibility principle is that Stage 05 has two logically separate branches:

- **historical v1.4.1 specification reconstruction** — reconstructs the documented historical algorithm from frozen canonical inputs and compares observed effect sizes with archived checkpoints.
- **canonical_transitivity_analysis** — preserves the spatial endpoints but uses sample-balanced inference and mechanism-neutral interpretation.

The original v1.4.1 source/result package and exact historical RNG stream are unavailable in the currently audited project environment. Therefore the historical branch must not be described as an exact source-code or Monte Carlo replication. Its observed effect-size regression passed; its historical permutation P/BH stream was not exactly reproduced and is retained as a provenance limitation.

Historical checkpoints are regression targets only, never analytical inputs.

Stage 05 is an **analysis-only stage**. It does not create a vdCHIBIN candidate score, does not alter Stage 02–04 ranking features, and does not assign 23- or 24-nt populations to a biochemical pathway.

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

### 3.0 Canonical evidence-building principle

The new pipeline rebuilds the biological analysis in a defined sequence from validated frozen inputs. Historical results may guide expectations and later serve as regression/validation targets, but they are **not treated as conclusions already established by the new pipeline**.

For each stage:

1. define the analysis from the authoritative specification;
2. calculate the result from the validated upstream inputs;
3. interpret the newly generated result;
4. only then compare it with historical outputs.

If a canonical result differs materially from the historical analysis, the discrepancy must be investigated and reported rather than hidden or forced to match. Historical downstream result tables must never be used as inputs to recreate the canonical result they are supposed to validate.

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
references/consensus/<sample>.<analysis_unit>.final.fa
references/consensus/<sample>.<analysis_unit>.final.background_masked.fa
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

# 01 — Viral length landscape and fixed 23/24-nt populations

## Purpose

Stage 01 is the first biological analysis produced by the canonical pipeline.

It must first reconstruct the viral small-RNA length landscape from the validated read-level inputs and only then perform the dedicated 23/24-nt comparison.

Historical knowledge that 23- and 24-nt vsiRNAs are prominent is treated as an expectation to test, not as an assumption built into the calculation.

Stage 01 therefore has two linked parts:

- **01A — 15–35-nt viral small-RNA length landscape**
- **01B — dedicated 23/24-nt population and strand analysis**

No Dicer, terminal-enrichment, transitivity, accessibility, or candidate-design inference is performed in Stage 01.

---

## 01A — Reconstruct the 15–35-nt viral small-RNA length landscape

### Biological question

Across primary-eligible Varroa virus infections, which small-RNA lengths are prominent, and is that pattern supported both by accumulated abundance and by diversity of distinct RNA sequences?

### Why 15–35 nt?

The corrected upstream preprocessing retained the 15–35-nt small-RNA range. Stage 01 reconstructs the complete retained viral length spectrum rather than beginning by filtering to 23 and 24 nt.

Rows outside 15–35 nt, if encountered in the validated read-level tables, are reported in Stage 01 QC and are not included in the canonical length-spectrum denominator.

### Inputs

Validated:

```text
tables/<sample>/<sample>.read_level_features.tsv.gz
results/descriptive/eligibility.tsv
```

Metadata required for reporting are taken from the eligibility table, including:

- `sample`;
- `analysis_unit`;
- `biological_virus`;
- `polarity`;
- `primary_eligible`.

No remapping is performed.

### Canonical inclusion criteria

For the length-spectrum analysis, retain rows satisfying:

```text
sample × virus corresponds to a primary_eligible sample × analysis_unit pair
mapping_mode = exact
virus_assignment = assigned
strand ∈ {sense, antisense}
length ∈ [15, 35]
```

The read-level `virus` field is matched to eligibility `analysis_unit`.

Cross-virus ambiguous reads are excluded from virus-specific primary summaries.

### Weighting mode A — abundance

For each eligible row, abundance contribution is the numeric read-level `count` field.

Do not use the number of table rows as a substitute for abundance.

For each sample-virus unit and length `L`:

```text
length_count_abundance(L)
    = sum(count for eligible rows of length L)
```

### Weighting mode B — true unique sequence

For Stage 01, a distinct sequence is defined within:

```text
sample × analysis_unit × length × strand
```

Each distinct sequence contributes exactly 1 within that unit, regardless of its read-level abundance or how many read names carry it.

The same literal sequence may count again in another sample or another biological analysis unit because the biological unit is different.

For a length-spectrum total, unique sense and unique antisense sequence counts are summed after strand-specific deduplication.

### Length fraction

Raw counts are retained, but the primary comparable length-spectrum quantity is the within-pair fraction:

```text
length_fraction(L)
    = length_count(L)
      / sum(length_count(15), ..., length_count(35))
```

This is calculated separately for abundance and unique-sequence weighting.

If the 15–35-nt denominator is zero, all length fractions for that unit are `NA` and the unit is reported in QC.

### Length rank

Within each sample-virus unit and weighting mode, lengths are ranked from highest to lowest `length_fraction`.

Use **standard competition ranking with the minimum rank for ties**:

```text
1, 2, 2, 4, ...
```

Thus tied lengths receive the same best applicable rank and no arbitrary genomic or lexical tie-breaking is introduced.

A length is descriptively `top1` when `rank <= 1` and `top3` when `rank <= 3`.

Top-1/top-3 frequencies are descriptive robustness summaries, not formal hypothesis tests.

### Pair-level output

For every:

```text
sample × analysis_unit × weighting_mode × length
```

retain at least:

- sample;
- analysis unit;
- biological virus;
- polarity;
- weighting mode;
- length;
- length count;
- length fraction;
- length rank;
- top-1 indicator;
- top-3 indicator.

### Sample-level aggregation

Several viruses can occur in the same biological sample. For the canonical sample-balanced view, calculate within each sample and weighting mode:

- median `length_fraction` across eligible sample-virus units for each length;
- median `length_rank` across eligible sample-virus units for each length;
- number of contributing sample-virus units.

This preserves virus-level calculations while preventing samples containing more viruses from automatically contributing more independent weight to the cross-dataset summary.

### Across-dataset summary

For each length and weighting mode, report at minimum:

- sample-balanced median length fraction;
- sample-clustered 95% bootstrap confidence interval for that median;
- sample-balanced median rank;
- number of contributing biological samples;
- number of contributing sample-virus units;
- pair-level fraction of sample-virus units in which the length is top 1;
- pair-level fraction of sample-virus units in which the length is top 3.

The top-1/top-3 pair-level frequencies are explicitly descriptive and must not replace the sample-balanced primary summary.

### Statistical emphasis

Stage 01A is primarily descriptive. It does not require a P-value to decide which lengths are prominent.

Prominence is evaluated from effect-size summaries, distributions, ranks, reproducibility across units, and sample-clustered uncertainty.

Only after the canonical Stage 01A outputs have been generated may they be compared with the historical length-distribution results.

---

## 01B — Dedicated 23/24-nt population and strand analysis

### Biological question

Once the full length landscape has been reconstructed, how do the 23-nt and 24-nt viral small-RNA populations differ in abundance, sequence diversity, and strand bias?

### Inclusion criteria

Use the same Stage 01A canonical population definition, then restrict to:

```text
length ∈ {23, 24}
```

No additional biological eligibility filter is introduced.

### Four primary populations

For every eligible sample-virus unit, preserve separately:

```text
23 sense      (23S)
23 antisense  (23AS)
24 sense      (24S)
24 antisense  (24AS)
```

Calculate them under both abundance and unique-sequence weighting.

### Counts

For each weighting mode, retain:

```text
n23_sense
n23_antisense
n23_total
n24_sense
n24_antisense
n24_total
```

In abundance mode these are sums of read-level `count`.

In unique-sequence mode these are numbers of distinct sequences under the Stage 01 uniqueness definition.

### Strand fractions

For each sample-virus unit and weighting mode:

```text
antisense_fraction_23
    = 23AS / (23S + 23AS)

antisense_fraction_24
    = 24AS / (24S + 24AS)
```

If the relevant denominator is zero, report `NA`.

Do **not** call `antisense_fraction_24` `F24_AS`; the latter name is reserved in Stage 05 for the different quantity describing 24-nt composition within the antisense 23+24 population.

### Direct 24-vs-23 strand-bias effect size

Calculate:

```text
delta_antisense_fraction_24_minus_23
    = antisense_fraction_24
      - antisense_fraction_23
```

Interpretation:

```text
>0 = the 24-nt population is more antisense-biased than the 23-nt population
 0 = equal antisense fraction
<0 = the 24-nt population is less antisense-biased
```

This is a descriptive population effect size, not a mechanistic classification of individual RNAs.

### 23/24 composition

For descriptive context, also calculate:

```text
length23_fraction_among_23_24
    = n23_total / (n23_total + n24_total)

length24_fraction_among_23_24
    = n24_total / (n23_total + n24_total)
```

If `n23_total + n24_total = 0`, both are `NA`.

These quantities describe the relative representation of the two lengths and are distinct from Stage 05 `F24_AS`.

### Sample-level aggregation

Within each biological sample and weighting mode, take the median across eligible sample-virus units for continuous pair-level metrics such as:

- `antisense_fraction_23`;
- `antisense_fraction_24`;
- `delta_antisense_fraction_24_minus_23`;
- `length23_fraction_among_23_24`;
- `length24_fraction_among_23_24`.

Retain the number of contributing virus units.

Raw counts remain available at pair level and may be summarized descriptively, but the primary cross-dataset biological comparison should emphasize fractions/effect sizes rather than pooled raw read totals.

### Across-dataset aggregation

For each weighting mode, report sample-balanced point estimates and sample-clustered 95% bootstrap confidence intervals for the principal fraction/effect-size metrics.

Retain:

- number of contributing samples;
- number of contributing sample-virus units;
- pair-balanced descriptive summaries for comparison;
- no automatic formal significance test unless a later pre-specified biological question requires one.

### Interpretation limits

Stage 01 may establish observations such as:

- 23 and/or 24 nt are prominent in the newly reconstructed length landscape;
- one length class is more abundant or sequence-diverse than another;
- 24 nt is more or less antisense-biased than 23 nt;
- these patterns recur across samples and viruses.

Stage 01 must **not** conclude from length and strand bias alone that:

```text
23 nt = primary Dicer products
24 nt = secondary RdRP products
all 23-mers are primary
all 24-mers are secondary
```

Those mechanistic questions belong to later stages.

---

## Stage 01 QC/accounting

Write a compact accounting table documenting at minimum:

- number of samples represented in the eligibility table;
- number of primary-eligible sample-virus units;
- read-level rows examined;
- exact + assigned rows retained;
- abundance retained in the 15–35-nt range;
- distinct Stage 01 sequences retained;
- rows outside the canonical 15–35-nt range;
- 23-nt abundance and unique-sequence totals;
- 24-nt abundance and unique-sequence totals;
- zero-denominator/zero-signal units;
- any unexpected categorical value encountered.

Unexpected categories must be reported rather than silently coerced.

---

## Stage 01 outputs

```text
results/01_viral_23_24/
│
├── qc/
│   └── stage01_accounting.tsv
│
├── length_spectrum/
│   ├── length_distribution_by_pair.tsv
│   ├── length_distribution_by_sample.tsv
│   └── length_distribution_across_dataset.tsv
│
├── fixed_23_24/
│   ├── 23_24_counts_by_pair.tsv
│   ├── 23_24_fractions_by_pair.tsv
│   ├── 23_24_by_sample.tsv
│   └── 23_24_across_dataset.tsv
│
└── figures/
```

### Minimum figures

Generate only figures with a clear scientific purpose:

1. **15–35-nt length spectrum** — abundance and unique-sequence views, showing the newly reconstructed length landscape without assuming 23/24 dominance.
2. **23-vs-24 antisense-fraction comparison** — pair-level `antisense_fraction_23` versus `antisense_fraction_24`, with the equality line shown for interpretation.
3. **23/24 strand-bias distributions** — distributions of the two antisense fractions, separated clearly by weighting mode.

Avoid redundant cosmetic variants.

---

## Historical comparison rule

Historical length-distribution and fixed-23/24 results are **not inputs** to Stage 01.

After all canonical Stage 01 outputs are calculated and frozen for the run, historical outputs may be used as regression/interpretive checks.

A material discrepancy must be documented and investigated. The canonical calculation must never be altered solely to force agreement with the historical result.

---

# 02 — Length-matched terminal nucleotide enrichment

## Purpose

Stage 02 measures whether terminal nucleotides in the newly reconstructed 23-nt and 24-nt viral small-RNA populations occur more or less often than expected from the viral sequence space that was actually available in each infection.

Stage 02 does **not** assume that any terminal preference exists. It reconstructs observed and expected terminal composition independently, calculates enrichment within each biological sample-virus unit, and only then summarizes the pattern across samples.

The primary cross-dataset result is sample-balanced. A separate **pooled-abundance** result is retained as a secondary descriptive view of the total accumulated molecular pool; it must never replace the sample-balanced biological summary.

## Biological question

Do observed 23-nt and 24-nt *Varroa* viral small RNAs contain terminal nucleotide patterns more or less often than expected from the viral sequence actually available for processing?

Secondary descriptive question:

> Which terminal patterns dominate the total accumulated eligible viral small-RNA molecules when sequencing abundance is allowed to weight the complete dataset?

These are different questions and are reported separately.

---

## 02.1 Inputs

Use only validated frozen-core inputs:

```text
tables/<sample>/<sample>.read_level_features.tsv.gz
results/descriptive/eligibility.tsv
references/consensus/<sample>.<analysis_unit>.final.background_masked.fa
```

Stage 01 outputs may be used for regression/QC comparison of population totals, but Stage 02 must be calculable directly from the validated frozen inputs and must not require historical terminal-enrichment outputs.

No remapping is performed.

---

## 02.2 Canonical inclusion criteria

Retain observed read-level rows only when:

```text
sample × virus matches a primary_eligible sample × analysis_unit pair
mapping_mode = exact
virus_assignment = assigned
strand ∈ {sense, antisense}
length ∈ {23, 24}
```

Cross-virus ambiguous reads are excluded from virus-specific primary summaries.

Observed sequences are interpreted in the DNA-alphabet representation present in FASTQ/FASTA files:

```text
A, C, G, T
```

Biologically, `T` corresponds to uridine (`U`) in the RNA molecule. Canonical TSV outputs use `T` consistently to avoid silent T/U conversion inside the workflow.

Observed sequence strings must agree with their declared length. Unexpected non-ACGT bases in an included observed 23/24-nt sequence are reported as a data-integrity failure rather than silently coerced.

---

## 02.3 Terminal coordinates

For every physical RNA sequence in its own 5′→3′ orientation:

```text
5p1 = first nucleotide
5p2 = second nucleotide
3p2 = penultimate nucleotide
3p1 = final nucleotide
```

Example:

```text
5′ A C ........ T G 3′
   ↑ ↑          ↑ ↑
 5p1 5p2      3p2 3p1
```

The names `3p2` and `3p1` are therefore defined from the physical 3′ end, not from left-to-right reference coordinates.

### Observed antisense orientation

An observed antisense read-level sequence is already represented as the sequenced RNA in its own 5′→3′ orientation.

Therefore:

```text
observed antisense sequence → use directly
```

Do **not** reverse-complement observed antisense reads before extracting terminal nucleotides.

### Expected antisense orientation

Expected antisense RNAs are generated from reference-orientation viral windows:

```text
reference window
    ↓
reverse complement
    ↓
hypothetical antisense RNA in 5′→3′ orientation
```

This is where reverse complementation belongs.

---

## 02.4 Weighting modes for the observed population

Calculate all primary Stage 02 enrichment quantities under both Stage 01 weighting modes.

### Abundance weighting

Each included observed RNA contributes its numeric read-level `count`.

For one sample-virus × length × strand population:

```text
observed_total_weight
    = sum(count)
```

A read with `count = 1000` contributes 1000 times as much observed molecular abundance as a read with `count = 1`.

### Unique-sequence weighting

Use the same Stage 01 sequence identity:

```text
sample × analysis_unit × length × strand × sequence
```

Each distinct sequence contributes exactly 1 within that unit, regardless of read abundance or number of read names carrying it.

Sense and antisense sequences are deduplicated separately before any combined-strand calculation.

### Expected-background weighting

The matched viral background is a **sequence-opportunity background**. Every fully supported genomic start position contributes one candidate window. Repeated identical window sequences at different supported positions therefore contribute once per genomic position.

This same positional opportunity background is used for abundance and unique-sequence observed modes. The two modes differ in how the **observed** RNA population is weighted; they do not redefine the viral genome itself.

---

## 02.5 Observed terminal frequency

For every:

```text
sample
× analysis_unit
× length {23,24}
× strand {sense,antisense}
× weighting_mode {abundance,unique_sequence}
× terminal_position {5p1,5p2,3p2,3p1}
× nucleotide {A,C,G,T}
```

calculate:

```text
observed_fraction
    = observed weight carrying nucleotide b at position p
      / total observed weight for that population
```

If the observed population denominator is zero, `observed_fraction` is `NA` for all four nucleotides at that position and the zero-signal unit is reported in QC.

For every valid non-zero population and position:

```text
Σ observed_fraction(A,C,G,T) = 1
```

within numerical tolerance.

---

## 02.6 Matched expected viral background

For every primary-eligible sample-virus unit and target length `L ∈ {23,24}`, enumerate all length-`L` windows independently within each FASTA record of the corresponding sample-specific depth-masked consensus.

A window is **fully supported** only if every base in the window is one of:

```text
A, C, G, T
```

Any window containing `N` or another non-ACGT character is excluded from the expected background.

Windows must never cross FASTA-record/contig boundaries.

### Sense expectation

For every fully supported reference-orientation window:

```text
expected sense RNA = reference window
```

Extract `5p1`, `5p2`, `3p2`, and `3p1` directly.

### Antisense expectation

For the same fully supported reference window:

```text
expected antisense RNA = reverse_complement(reference window)
```

Then extract terminal positions from that antisense sequence in its own 5′→3′ orientation.

### Expected fraction

For nucleotide `b` and position `p`:

```text
expected_fraction(b,p)
    = number of fully supported candidate windows carrying b at p
      / total number of fully supported candidate windows
```

Sense and antisense expected fractions are calculated separately.

If no fully supported window exists for a required sample-virus × length background, the corresponding expected frequencies and enrichment values are `NA` and the unit is prominently reported in QC.

For every valid background and position:

```text
Σ expected_fraction(A,C,G,T) = 1
```

within numerical tolerance.

---

## 02.7 Combined-strand observed and expected frequencies

Retain three strand scopes:

```text
sense
antisense
combined
```

The combined observed population is the direct combination of the included sense and antisense observed populations under the same weighting mode.

The combined expected background must **not** be forced to a 50:50 sense/antisense mixture.

For each sample-virus × length × weighting mode, define:

```text
wS  = observed sense weight / (observed sense weight + observed antisense weight)
wAS = observed antisense weight / (observed sense weight + observed antisense weight)
```

when the denominator is non-zero.

Then:

```text
expected_fraction_combined(b,p)
    = wS  × expected_fraction_sense(b,p)
      + wAS × expected_fraction_antisense(b,p)
```

Thus the expected combined composition is matched to the actual strand mixture of that biological unit and weighting mode.

If the combined observed denominator is zero, `wS`, `wAS`, combined observed frequencies, combined expected frequencies, and combined enrichment are `NA`.

---

## 02.8 Pair-level enrichment ratio

For every valid sample-virus × length × strand scope × weighting mode × terminal position × nucleotide:

```text
enrichment_ratio
    = observed_fraction / expected_fraction
```

Interpretation:

```text
1   = observed as often as sequence availability predicts
>1  = enriched
<1  = depleted
0   = nucleotide was not observed although it was available in the background
```

If:

```text
expected_fraction = 0
```

then:

```text
enrichment_ratio = NA
```

No arbitrary pseudocount is added.

If the observed denominator is zero, the enrichment is also `NA` rather than zero.

Stage 02 enrichment is an empirical association and does not identify which molecular process produced it.

---

## 02.9 Canonical sample-balanced aggregation

The primary cross-dataset summary treats the biological sequencing sample as the top-level replication unit.

For each fixed:

```text
length
× strand_scope
× weighting_mode
× terminal_position
× nucleotide
```

perform:

1. calculate `enrichment_ratio` separately for every eligible sample-virus unit;
2. within each sample, take the median of finite enrichment ratios across its eligible viruses;
3. across samples, take the median of those sample-level medians.

This produces:

```text
sample_balanced_median_enrichment_ratio
```

The point estimate is the **median of matched pair-level ratios**, not a ratio of separately aggregated observed and expected frequencies.

### Canonical uncertainty

Use a sample-clustered percentile bootstrap:

1. resample sample IDs with replacement;
2. when a sample is selected, retain its eligible virus observations together;
3. recompute its within-sample median;
4. recompute the across-sample median;
5. report the 2.5th and 97.5th percentiles.

Record the seed, requested bootstrap replicates, valid replicates, and interval method in the output/provenance.

No read-level or virus-row bootstrap is substituted for the sample-clustered primary uncertainty estimate.

---

## 02.10 Pair-balanced historical/regression summary

Also retain:

```text
pair_median_enrichment_ratio
    = median of finite enrichment_ratio values
      across eligible sample-virus units
```

This is useful for historical regression because it corresponds most closely to the previous design-facing `median_enrichment_ratio` logic.

It is **not** the canonical primary cross-dataset inference because samples containing multiple eligible viruses contribute multiple pair observations.

If a historical-compatible field named exactly `median_enrichment_ratio` is exported, it must be explicitly documented as the pair-balanced legacy-compatible quantity; the canonical sample-balanced field must have a different unambiguous name.

---

## 02.11 Secondary pooled-abundance description

In addition to the primary sample-balanced abundance analysis, calculate one explicitly secondary molecular-pool view for **abundance weighting only**.

This asks:

> What terminal enrichment characterizes the total accumulated eligible molecules across the complete dataset when high-abundance/high-depth infections are allowed to contribute proportionally more molecular weight?

It does **not** estimate the typical biological sample.

For each fixed:

```text
length
× strand_scope
× terminal_position
× nucleotide
```

let, for sample-virus unit `u`:

```text
N_u = total observed abundance of the relevant length/strand scope
O_u = observed abundance carrying nucleotide b at position p
E_u = matched expected_fraction(b,p) for that unit
```

using only units for which the matched expected fraction is defined.

Calculate:

```text
pooled_abundance_observed_fraction
    = Σ O_u / Σ N_u
```

and the abundance-matched expected frequency:

```text
pooled_abundance_expected_fraction
    = Σ (N_u × E_u) / Σ N_u
```

Then:

```text
pooled_abundance_enrichment_ratio
    = pooled_abundance_observed_fraction
      / pooled_abundance_expected_fraction
```

For `combined`, `E_u` is the abundance-mode observed-strand-weighted combined expectation defined above.

This formulation preserves each unit's matched viral background while allowing its accumulated molecular abundance to determine its contribution to the pooled molecular view.

If the pooled abundance denominator is zero or the pooled expected fraction is zero, report `NA`.

### Interpretation and inferential status

`pooled_abundance_enrichment_ratio` is **secondary descriptive output only**.

It must:

- be labelled `pooled_abundance` or equivalent in every output;
- report total contributing observed abundance, number of contributing sample-virus units, and number of contributing samples;
- not replace `sample_balanced_median_enrichment_ratio` as the primary biological summary;
- not be given a read-level P-value or treated as though millions of reads were millions of independent biological replicates;
- not automatically become the later vdCHIBIN design reference without an explicit Stage 08/09 decision.

No formal confidence interval is required for this pooled descriptive quantity. If a clustered sensitivity interval is ever added later, it must resample biological samples, not individual reads.

---

## 02.12 23-versus-24 enrichment-landscape comparison

Only after the 23-nt and 24-nt enrichment landscapes have been independently reconstructed, compare them using matched terminal features.

Each length has:

```text
4 terminal positions × 4 nucleotides = 16 matched features
```

For the canonical comparison, calculate Spearman rank correlation between the 23- and 24-nt **sample-balanced median enrichment ratios** for the same 16 features.

At minimum report separately for:

```text
combined strand scope
antisense strand scope
```

and separately for abundance and unique-sequence weighting.

The antisense comparison is especially relevant to later guide-oriented design, but Stage 02 itself does not convert it into a candidate score.

A pooled-abundance 23-vs-24 correlation may be reported descriptively if useful, but it is not required for the canonical Stage 02 result.

Correlation measures similarity of the enrichment landscapes; it does not establish a shared enzyme or processing pathway.

---

## 02.13 QC/accounting

Write a compact Stage 02 accounting table reporting at minimum:

- samples represented;
- primary-eligible sample-virus units;
- observed read-level rows examined;
- exact + assigned eligible 23/24 rows retained;
- retained abundance and distinct-sequence totals by length and strand;
- observed zero-signal populations;
- observed length mismatches;
- unexpected observed bases/categories;
- number of background FASTA records examined;
- number of fully supported 23-nt and 24-nt windows per sample-virus unit;
- units with zero valid background windows;
- excluded background windows containing masked/ambiguous sequence;
- maximum deviation of valid observed A/C/G/T frequency sums from 1;
- maximum deviation of valid expected A/C/G/T frequency sums from 1;
- number of finite/non-finite pair-level enrichment values;
- number of samples contributing to each canonical dataset-level feature;
- total molecular abundance contributing to each pooled-abundance result.

Masked `N` bases in background consensuses are expected and cause overlapping windows to be excluded. Unexpected non-ACGT bases in included observed RNA sequences are not silently accepted.

---

## 02.14 Required outputs

```text
results/02_terminal_enrichment/
│
├── qc/
│   └── stage02_accounting.tsv
│
├── observed/
│   └── terminal_observed_by_pair.tsv
│
├── background/
│   └── terminal_expected_by_pair.tsv
│
├── enrichment/
│   ├── terminal_enrichment_by_pair.tsv
│   ├── terminal_enrichment_by_sample.tsv
│   ├── terminal_enrichment_across_dataset.tsv
│   └── terminal_enrichment_pooled_abundance.tsv
│
├── comparisons/
│   └── enrichment_23_vs_24.tsv
│
└── figures/
```

### `terminal_enrichment_across_dataset.tsv`

Must include at minimum:

```text
length
strand_scope
weighting_mode
terminal_position
nucleotide
sample_balanced_median_enrichment_ratio
ci_low
ci_high
n_samples
n_sample_virus_units
pair_median_enrichment_ratio
bootstrap_replicates_requested
bootstrap_replicates_valid
bootstrap_seed
ci_method
```

### `terminal_enrichment_pooled_abundance.tsv`

Must include at minimum:

```text
length
strand_scope
terminal_position
nucleotide
pooled_abundance_observed_fraction
pooled_abundance_expected_fraction
pooled_abundance_enrichment_ratio
pooled_observed_total_weight
n_samples
n_sample_virus_units
analysis_role = secondary_descriptive
```

---

## 02.15 Minimum figures

Figures are generated only after numerical Stage 02 outputs have passed review.

Minimum useful figures are:

1. **Terminal-enrichment landscape** — clear 23-nt and 24-nt terminal-feature visualization, with abundance and unique-sequence views separated rather than overplotted.
2. **23-vs-24 enrichment comparison** — matched 16-feature scatter for the canonical sample-balanced enrichment values, including the reported Spearman rho; antisense should be clearly available as a dedicated view.

The pooled-abundance result does not require a separate figure unless it materially changes interpretation relative to the sample-balanced abundance result.

Avoid redundant cosmetic variants.

---

## 02.16 Historical comparison rule

Historical fixed-length terminal-enrichment outputs are **not inputs** to the canonical Stage 02 calculation.

After the complete canonical Stage 02 outputs have been generated and interpreted, historical outputs may be read for regression comparison.

Numerical comparisons are made only when eligibility, weighting, denominator, strand scope, background definition, and aggregation are genuinely comparable.

If a historical quantity is pair-balanced while the canonical quantity is sample-balanced, label it **not directly comparable** rather than presenting a misleading numerical difference.

A material discrepancy in directly comparable pair-level quantities must be investigated and documented; the canonical calculation must not be changed solely to force historical agreement.

---

## 02.17 Interpretation limits

Stage 02 may establish observations such as:

- a terminal nucleotide is enriched or depleted relative to matched viral sequence availability;
- an enrichment is reproduced across biological samples;
- the pattern persists or changes between abundance and unique-sequence weighting;
- pooled molecular abundance emphasizes the same or a different terminal pattern;
- 23- and 24-nt terminal landscapes are similar or different.

Stage 02 must **not** by itself conclude that a preference was caused specifically by:

```text
Dicer cleavage
Argonaute loading
strand selection
RdRP synthesis
RNA stability
```

or that a terminal nucleotide guarantees RNAi efficacy.

It is an empirical sequence-association layer. Mechanistic Dicer analysis begins in Stage 03, and candidate scoring is deliberately deferred to later vdCHIBIN stages.

---


# 03 — Official stepRNA duplex-geometry reconstruction

## Purpose and biological question

Stage 03 asks whether the canonical 23-nt and/or 24-nt viral small-RNA populations contain recoverable complementary opposite-strand partners whose duplex-end geometry is consistent with Dicer/Dicer-like processing.

The primary question is:

> Among canonically defined viral small RNAs, what complementary passenger strands can be recovered, what 5′/3′ end-distance geometries do the reconstructed duplexes show, and what passenger lengths are represented?

Stage 03 is a **measurement/reconstruction stage**. It does not yet make the canonical cross-dataset claim that one length class has stronger Dicer evidence than another. Sample-aware biological aggregation and the direct 23-versus-24 Dicer comparison belong to Stage 04.

The primary published method is:

> Murcott B, Pawluk RJ, Protasio AV, Akinmusola RY, Lastik D, Hunt VL. 2022. *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. Frontiers in Bioinformatics 2:994871. DOI: 10.3389/fbinf.2022.994871.

---

## 03.1 Method priority

Use the official published **stepRNA** implementation as the primary duplex-geometry method.

The canonical environment pins:

```text
stepRNA == 1.0.6
```

unless a future explicit project decision changes the version and updates this specification.

The exact versions of Python, Bowtie2, pysam, Biopython, NumPy, and other runtime dependencies actually used must be recorded in Stage 03 provenance.

The project-native/historical Dicer reconstruction is **not** used to generate the primary Stage 03 result. It is retained for Stage 04 secondary validation and historical regression only.

There is no silent fallback to a custom stepRNA reimplementation if the official tool cannot be executed reproducibly.

---

## 03.2 Mandatory software preflight

Before running the biological dataset, Stage 03 must perform a cheap preflight that records:

- installed stepRNA version;
- Bowtie2 version;
- Python version;
- relevant Python package versions;
- successful invocation of the official stepRNA executable;
- successful construction/use of a Bowtie2 index by stepRNA;
- a synthetic or bundled-example check sufficient to confirm that the official output can be generated and parsed;
- a sign-convention check confirming that the parser interprets a known reference 5′ underhang as positive and a known reference 3′ overhang as negative.

If the preflight fails, Stage 03 stops with an explicit failure. It must not start a large batch of biological runs while repeatedly attempting undocumented workarounds.

---

## 03.3 Canonical biological scope

Use the same biological eligibility established upstream.

A sample-virus unit is included only when:

```text
primary_eligible == true
```

Observed reads used to construct Stage 03 inputs must additionally satisfy:

```text
mapping_mode == exact
virus_assignment == assigned
strand in {sense, antisense}
```

No remapping to the virus is performed in Stage 03.

The read-level `virus` field must match eligibility `analysis_unit`.

Unexpected or inconsistent categories are reported rather than silently coerced.

---

## 03.4 File A — focal/reference populations

Generate File A separately for each eligible sample-virus unit and each focal class:

```text
23-nt sense
23-nt antisense
24-nt sense
24-nt antisense
```

A File-A population with no eligible focal sequence is not sent to stepRNA; it is recorded as a zero-signal population in Stage 03 QC.

### Canonical File-A representation

The biological stepRNA run uses **collapsed distinct focal sequences**:

```text
one FASTA record per distinct
sample × analysis_unit × focal_length × focal_strand × sequence
```

Each FASTA record receives a stable opaque identifier.

The corresponding observed abundance from the read-level `count` field is stored in a separate manifest and is not encoded in a way that requires stepRNA to infer project-specific abundance semantics from the FASTA header.

This collapsed-input policy follows the published use of collapsed sRNA sequence data and avoids creating enormous abundance-expanded FASTA files.

---

## 03.5 File B — potential passenger pool

For each File-A run, File B contains potential passenger sequences from:

```text
the same sample
× the same analysis_unit
× the opposite mapped viral strand
× length 15–30 nt inclusive
```

Only exact, assigned reads are eligible.

File B is also collapsed to one FASTA record per distinct passenger sequence within that run, with stable opaque identifiers and a manifest retaining observed abundance.

The canonical primary passenger-length range is therefore:

```text
15–30 nt
```

This broad range is intentionally wider than the 23/24-nt focal classes so that a 23- or 24-nt focal RNA can recover a complementary passenger of a different length.

No requirement is imposed that the passenger itself be 23 or 24 nt.

---

## 03.6 Strand/orientation handling

File-A and File-B sequences are written in their **observed physical 5′→3′ sequence orientation**.

The pipeline must not manually reverse-complement the observed passenger FASTA before calling stepRNA.

Complementary orientation is determined by the official stepRNA/Bowtie2 alignment procedure.

Restricting File B to the opposite **mapped viral strand** defines the biological candidate pool; it does not change the physical sequence orientation supplied to stepRNA.

---

## 03.7 Official alignment behaviour

The canonical run preserves the official stepRNA alignment logic:

- Bowtie2 index built from File A;
- File B aligned to File A;
- local alignment;
- no mismatches in the complementary aligned segment under the official default method;
- terminal soft clipping retained so overhang/underhang distances can be reconstructed.

Do not replace this with genome-coordinate pairing, a custom local aligner, or mismatch-tolerant matching in the primary Stage 03 analysis.

Any future mismatch-tolerant analysis is exploratory and must have a separate name, configuration, and output path.

---

## 03.8 Official signed-distance convention

All parsed Stage 03 outputs preserve the official stepRNA sign convention relative to the File-A focal/reference RNA:

```text
negative distance = File-A reference overhang
0                 = blunt end
positive distance = File-A reference underhang
```

The 5′ and 3′ distances must always be stored in separate explicit fields:

```text
steprna_5p_distance
steprna_3p_distance
```

No sign conversion is applied merely to make a plot visually intuitive.

---

## 03.9 Passenger recovery is distinct from duplex geometry

Stage 03 reports passenger recovery separately from geometry.

### Unique-reference passenger recovery

```text
passenger_recovery_fraction_unique
    = number of distinct File-A focal sequences with ≥1 recovered passenger
      / number of distinct eligible File-A focal sequences
```

### Abundance-weighted passenger recovery

```text
passenger_recovery_fraction_abundance
    = sum(observed focal abundance for File-A sequences with ≥1 recovered passenger)
      / sum(observed focal abundance for all eligible File-A sequences)
```

These two quantities ask different questions:

- unique mode: how broadly passenger recovery is represented across focal sequence diversity;
- abundance mode: what fraction of the accumulated focal molecular population belongs to sequences for which a passenger is recoverable.

Low passenger recovery does **not** by itself demonstrate absence of Dicer processing. Passenger strands may be unstable, under-sampled, filtered, or otherwise unrecovered.

---

## 03.10 Marginal and same-duplex geometry are complementary Stage 03 outputs

Stage 03 must retain **both**:

1. the complete signed 5′ and 3′ **marginal end-distance spectra** produced from official stepRNA output; and
2. the complete **same-duplex joint geometry spectrum** reconstructed from the classified official stepRNA alignments.

These answer different questions and must never be conflated.

### Marginal end-distance spectrum

For each biological run retain, at minimum:

```text
sample
analysis_unit
biological_virus
focal_length
focal_strand
end                 # 5p or 3p
signed_distance
official_duplex_count
official_unique_reference_count
official stepRNA enrichment/log-ratio field(s)
official stepRNA Wald Z field(s)
```

A marginal distance of `0` means that one analysed end is flush/blunt relative to File A. It does **not** imply that the opposite end of the same duplex is also `0`.

### Same-duplex joint geometry spectrum

For every recovered focal/passenger duplex, retain the pair:

```text
(steprna_5p_distance, steprna_3p_distance)
```

using the two values derived from the **same classified alignment** through the already validated official-distance parser.

For each biological run and each observed `(d5,d3)` combination retain:

```text
sample
analysis_unit
biological_virus
focal_length
focal_strand
steprna_5p_distance
steprna_3p_distance
official_duplex_count
total_recovered_duplexes
joint_duplex_fraction
run_id
```

where:

```text
joint_duplex_fraction(d5,d3)
    = official_duplex_count(d5,d3)
      / total_recovered_duplexes
```

For every non-empty run:

```text
Σ joint_duplex_fraction(d5,d3) = 1
```

up to numerical precision.

The marginal and joint spectra are both required because a strong marginal `0` peak can be generated by many different same-duplex geometries, for example:

```text
(0,0), (0,-1), (0,+1), (0,-2), ...
```

at one end and:

```text
(0,0), (+1,0), (-1,0), (+2,0), ...
```

at the other.

Therefore:

> “distance 0 is prominent at individual ends”

must never be rewritten as:

> “most reconstructed duplexes are fully blunt.”

The geometry analysis remains descriptive of reconstructed sequencing relationships; it does not by itself identify the nuclease responsible for cleavage.

---

## 03.11 Native stepRNA counts versus project abundance weighting

The official stepRNA output files are preserved exactly.

In particular, the native stepRNA:

```text
*_overhang.csv
*_unique_overhang.csv
```

outputs must **not** be silently relabelled as the canonical project-wide `abundance_weighted` and `unique_sequence` metrics without validating what the official file counts represent.

For canonical Varroa abundance summaries, focal-sequence abundance is taken from the upstream read-level `count` field stored in the File-A manifest.

Stage 03 therefore keeps three concepts distinct:

1. **official duplex/alignment counts** produced by stepRNA;
2. **official unique-reference counts** produced by stepRNA;
3. **project abundance-weighted focal-reference support**, calculated from upstream focal-sequence abundance after geometry reconstruction.

This avoids accidentally multiplying or redefining abundance through passenger multiplicity.

Stage 04 will use the parsed reference/geometry support and focal abundance metadata for sample-aware population-level inference.

---

## 03.12 Passenger-length output

For every recovered focal/passenger relationship, retain the passenger length reported by stepRNA.

At pair level, produce a passenger-length summary by:

```text
sample
analysis_unit
focal_length
focal_strand
passenger_length
```

with official duplex/alignment counts and, where recoverable without ambiguity, unique focal-reference support counts.

Passenger-length distributions are descriptive Stage 03 outputs. They are not treated as proof that a particular nuclease generated the focal RNA.

---

## 03.13 Pre-specified Varroa 2-nt joint geometry

The historical Varroa analysis identified the following geometry as a feature of interest:

```text
5′ distance = +2
3′ distance = -2
```

under the official stepRNA sign convention.

This corresponds to:

```text
File-A 5′ underhang of 2 nt
File-A 3′ overhang of 2 nt
```

The canonical name is:

```text
varroa_2nt_joint_geometry
```

It is **pre-specified before inspecting the new Stage 03 result**.

It must not be called a universal “canonical Dicer geometry,” because Dicer/Dicer-like end geometries vary across organisms, pathways, and substrates.

### Joint-geometry support

The joint classification must use the 5′ and 3′ distances from the **same reconstructed focal/passenger duplex**.

Retain at least:

```text
n_recovered_duplexes
n_joint_geometry_duplexes
varroa_2nt_joint_duplex_fraction

n_focal_references
n_recovered_focal_references
n_focal_references_supporting_joint_geometry
varroa_2nt_reference_fraction_all
varroa_2nt_reference_fraction_recovered
```

A focal reference “supports” the joint geometry if at least one of its recovered passenger alignments has:

```text
5p = +2
3p = -2
```

The duplex-level and reference-level quantities must not be conflated.

For abundance weighting, also retain the fraction of focal molecular abundance represented by focal references that support the joint geometry. This is a **reference-support** metric, not a probability distribution over mutually exclusive distances; a focal reference can support more than one geometry through different recovered passengers.

---

## 03.14 No single Stage 03 “Dicer score”

Stage 03 does not collapse the evidence into one project-specific Dicer score.

The following remain separate:

- passenger recovery;
- full 5′ marginal distance spectrum;
- full 3′ marginal distance spectrum;
- full same-duplex `(d5,d3)` joint-geometry spectrum;
- official stepRNA enrichment/log-ratio output;
- official stepRNA Wald Z output;
- passenger-length distribution;
- pre-specified joint-geometry support.

The project-specific `Δ_Dicer` statistic belongs to Stage 04 as a secondary validation statistic.

---

## 03.15 Primary versus sensitivity passenger range

The required canonical Stage 03 run uses:

```text
15–30 nt
```

passengers.

A narrower:

```text
18–28 nt
```

analysis is a **named sensitivity analysis**, not part of the mandatory primary run.

It may be executed later if a robustness question makes it scientifically useful, but it must:

- use a separate output path;
- preserve all other settings;
- be labelled `passenger_18_28_sensitivity`;
- never replace the primary 15–30-nt result silently.

This prevents an optional sensitivity run from doubling the required Stage 03 computation without a defined need.

---

## 03.16 Required raw official outputs

For every successful official stepRNA run, preserve the native output directory or equivalent raw files needed for provenance, including the official summary files described by the published method:

- overhang counts;
- unique-reference overhang counts;
- passenger-number output;
- passenger-length output;
- overhang/underhang-type output;
- relevant classified BAM/alignment products needed to reconstruct joint 5′/3′ geometry.

Raw official outputs are immutable Stage 03 products and are not manually edited.

---

## 03.17 Canonical parsed outputs

Create stable project-level tables so Stage 04 does not depend on fragile stepRNA filenames.

Required parsed outputs:

```text
passenger_recovery_by_pair.tsv
overhang_spectrum_by_pair.tsv
passenger_length_by_pair.tsv
joint_geometry_spectrum_by_pair.tsv
joint_geometry_by_pair.tsv
joint_geometry_references.tsv.gz
```

`joint_geometry_references.tsv.gz` contains the focal reference identifiers/sequences needed to recover the pre-specified Dicer-supported subset in Stage 04 without rerunning stepRNA.

The parsed tables must retain enough identifiers to trace every summary back to:

```text
sample
analysis_unit
focal_length
focal_strand
stepRNA run ID
```

and, where reference-level support is reported, to the stable File-A focal identifier.

---

## 03.18 Required QC

Stage 03 QC must report at least:

- number of primary-eligible sample-virus units;
- number of possible File-A populations;
- number with non-zero focal input;
- distinct File-A sequences by focal length/strand;
- total focal abundance represented;
- distinct File-B passenger sequences by run;
- passenger-length range actually present;
- stepRNA version;
- Bowtie2 version;
- successful/failed stepRNA runs;
- zero-passenger runs;
- focal references with ≥1 passenger;
- maximum/minimum signed distances observed;
- passenger lengths observed;
- malformed/unparseable official output rows;
- runs missing an expected official output file;
- consistency of focal identifiers between manifests and parsed output;
- consistency of joint-geometry counts with the classified alignments used to derive them;
- exact per-run agreement between the sum of all same-duplex joint-geometry counts and `n_recovered_duplexes`;
- joint-geometry fraction sum equal to 1 for every non-empty run, within numerical precision;
- regression agreement between the `(+2,-2)` count recovered from the full joint spectrum and the pre-existing pre-specified `(+2,-2)` output.

A zero-passenger biological run is valid data and is not automatically a pipeline failure.

A software failure, missing required output, parser inconsistency, or identifier mismatch is a failure.

---

## 03.19 Outputs

```text
results/03_steprna/
│
├── qc/
│   ├── stage03_accounting.tsv
│   └── stage03_joint_geometry_spectrum_accounting.tsv
│
├── provenance/
│   ├── software_versions.tsv
│   └── run_manifest.tsv
│
├── inputs/
│   ├── input_manifest.tsv
│   ├── focal_reference_manifest.tsv.gz
│   └── passenger_manifest.tsv.gz
│
├── raw/
│   └── <one official stepRNA output directory per successful run>/
│
├── parsed/
│   ├── passenger_recovery_by_pair.tsv
│   ├── overhang_spectrum_by_pair.tsv
│   ├── passenger_length_by_pair.tsv
│   ├── joint_geometry_spectrum_by_pair.tsv
│   ├── joint_geometry_by_pair.tsv
│   └── joint_geometry_references.tsv.gz
│
└── sensitivity/
    └── passenger_18_28/        # absent unless explicitly run
```

Stage 03 does not need publication figures before the calculations and parser have been validated.

---

## 03.20 Interpretation limits

A reproducibly enriched duplex-end or same-duplex geometry is evidence **consistent with structured small-RNA processing**, and may be consistent with Dicer/Dicer-like processing, but geometry alone is not sufficient to assign a specific nuclease or distinguish primary Dicer products from secondary/RdRP-derived products.

Stage 03 does not directly observe cleavage in vivo and does not by itself:

- identify a specific Varroa Dicer/Dicer-like protein;
- prove that every RNA in a length class was Dicer-generated;
- prove that an RNA lacking a recovered passenger was not Dicer-generated;
- establish that 23 nt is “primary” and 24 nt is “secondary”;
- establish RdRP-dependent amplification;
- define a candidate-window scoring metric.

Those higher-level biological comparisons and Dicer-conditioned sequence analyses begin in Stage 04.

---


# 04 — Duplex-geometry evidence aggregation and geometry-conditioned sequence features

## 04.1 Purpose and analysis provenance

Stage 03 used the official stepRNA implementation to reconstruct complementary focal/passenger relationships and their signed 5′/3′ end geometry. Stage 04 does **not** rerun stepRNA. It asks two higher-level questions:

1. **Population evidence:** which stepRNA geometry features are reproducible across biological samples, and how do 23- and 24-nt populations compare?
2. **Sequence specificity:** do focal RNAs supporting the pre-specified Varroa `(+2,-2)` joint geometry have terminal sequence features that differ from the general Stage 02 population or from other passenger-recovered focal RNAs?

Stage 04 is deliberately named around **duplex geometry**, not around a definitive “Dicer score.” Reconstructed geometry can be consistent with Dicer/Dicer-like processing, but the sequencing data do not directly observe cleavage or identify the responsible nuclease.

### Relationship to the observed Stage 03 landscape

Stage 03 and the completed Stage 04 post-processing establish two complementary observations:

1. **Marginal view:** distance `0` is a prominent and frequently top-ranked 5′ or 3′ end-distance, particularly for antisense focal populations.
2. **Same-duplex view:** fully blunt `(0,0)` duplexes are nevertheless a minority because a marginal `0` at one end is commonly paired with a non-zero distance at the other end.

The validated sample-balanced `(0,0)` fractions are approximately:

```text
23S   0.0215
23AS  0.0931
24S   0.0227
24AS  0.0311
```

and all four focal classes have:

```text
0/54 runs with >50% of reconstructed duplexes exactly (0,0)
```

The pre-specified `(+2,-2)` geometry is also a minority:

```text
23S   0.0153
23AS  0.0200
24S   0.0198
24AS  0.0058
```

The marginal 3′ `-2` component remains reproducibly represented, especially in antisense populations, but the corresponding 5′ `+2` component is substantially weaker.

These observations **must not be used to redefine the pre-specified joint geometry after the fact**.

Therefore:

- the complete marginal 5′/3′ spectra remain required;
- the complete same-duplex `(d5,d3)` spectrum is also required;
- `(0,0)` is explicitly distinguished from marginal distance `0`;
- the `(+2,-2)` joint geometry remains a pre-specified secondary feature;
- neither `(0,0)` nor `(+2,-2)` is treated as a universal mechanistic signature;
- special emphasis on marginal 3′ `-2` remains post-Stage-03 exploratory/descriptive unless independently justified later.

---

## 04.2 Inputs

Use existing canonical outputs only. Do not remap reads and do not rerun stepRNA.

Primary Stage 03 inputs:

```text
results/03_steprna/parsed/passenger_recovery_by_pair.tsv
results/03_steprna/parsed/overhang_spectrum_by_pair.tsv
results/03_steprna/parsed/passenger_length_by_pair.tsv
results/03_steprna/parsed/joint_geometry_spectrum_by_pair.tsv
results/03_steprna/parsed/joint_geometry_by_pair.tsv
results/03_steprna/parsed/joint_geometry_references.tsv.gz
results/03_steprna/inputs/focal_reference_manifest.tsv.gz
```

Stage 02 inputs for matched terminal background and general enrichment:

```text
results/02_terminal_enrichment/background/terminal_expected_by_pair.tsv
results/02_terminal_enrichment/enrichment/terminal_enrichment_by_pair.tsv
results/02_terminal_enrichment/enrichment/terminal_enrichment_across_dataset.tsv
```

Eligibility/sample metadata remain defined by the validated frozen-core eligibility table already used upstream.

The frozen legacy core remains read-only.

---

## 04.3 Analysis units and focal-strand rule

Official stepRNA runs remain stratified as:

```text
23S
23AS
24S
24AS
```

Do **not** merge sense and antisense focal runs into one canonical geometry statistic. The focal-reference orientation changes the interpretation of signed distances, and passenger availability differs strongly by focal strand.

For later vdCHIBIN candidate development, `antisense` is the design-relevant focal strand because the eventual guide is antisense to the target mRNA. Sense-strand results remain an important biological comparison and must still be reported.

---

## 04.4 Canonical sample-aware aggregation of the full stepRNA spectrum

### Primary effect quantity

For each fixed:

```text
focal_length ∈ {23,24}
focal_strand ∈ {sense,antisense}
end ∈ {5p,3p}
signed_distance
official_view ∈ {duplex, unique_reference}
```

use the corresponding **official stepRNA log-ratio** from Stage 03 as the primary per-run effect quantity.

The native stepRNA Wald Z-score is retained as a run-level descriptive statistic but is **not combined into a new population-level Z-score or P-value**.

### Within-sample aggregation

For each fixed feature above:

```text
sample_steprna_log_ratio_median
    = median of finite pair-level official stepRNA log-ratios
      across eligible viruses within one sequencing sample
```

Retain the number of contributing sample-virus units.

### Across-sample aggregation

```text
sample_balanced_steprna_log_ratio
    = median of sample_steprna_log_ratio_median
      across contributing biological samples
```

Uncertainty uses the canonical sample-clustered percentile bootstrap:

- resample biological samples with replacement;
- preserve all virus observations belonging to a sampled library;
- recompute the within-sample median;
- recompute the across-sample median;
- use the configured fixed seed and replicate count;
- report the 95% percentile interval and number of valid replicates.

This is an estimation procedure. A pointwise 95% interval for every signed distance must not be presented as a simultaneous family-wise significance statement.

### Same-duplex joint-geometry aggregation

For every fixed:

```text
focal_length ∈ {23,24}
focal_strand ∈ {sense,antisense}
steprna_5p_distance = d5
steprna_3p_distance = d3
```

begin from the Stage 03 pair-level:

```text
joint_duplex_fraction(d5,d3)
```

and aggregate canonically:

```text
joint fraction within sample-virus
→ median across viruses within biological sample
→ median across biological samples
→ sample-clustered percentile bootstrap 95% CI
```

The canonical cross-dataset quantity is:

```text
sample_balanced_joint_duplex_fraction(d5,d3)
```

The full `(d5,d3)` table is retained. At minimum, explicitly report:

```text
(0,0)
(+2,-2)
```

along with:

- the most common joint geometry in each sample-virus run;
- the number/fraction of runs in which `(0,0)` is the most common joint geometry;
- the number/fraction of runs in which `joint_duplex_fraction(0,0) > 0.5`;
- a pooled duplex fraction as a **secondary descriptive** quantity only.

A geometry being the single most common category does not imply that it constitutes a majority of duplexes.

### Dominant/enriched distance summaries

Within each focal class/end/official view, report:

- distance with the highest sample-balanced median log-ratio;
- sample-balanced log-ratio for distance `0`;
- sample-balanced log-ratio for the pre-specified marginal component (`+2` at 5′ or `-2` at 3′);
- full signed spectrum in the output table.

Do not discard other distances merely because one distance is largest.

---

## 04.5 Passenger recovery remains a separate observability metric

For each focal class, sample-balance the Stage 03 quantities:

```text
passenger_recovery_fraction_unique
passenger_recovery_fraction_abundance
```

using the same:

```text
median across viruses within sample
→ median across samples
→ sample-clustered bootstrap 95% CI
```

Passenger recovery is **not** treated as a Dicer-activity score.

This separation is essential because a focal RNA can only be geometrically reconstructed if a complementary passenger is present in the sequenced library. Differences in sense/antisense population abundance can therefore alter passenger recoverability independently of cleavage mechanism.

---

## 04.6 Pre-specified Varroa joint-geometry population summaries

The pre-specified feature remains exactly:

```text
steprna_5p_distance = +2
steprna_3p_distance = -2
```

using the two distances from the **same reconstructed focal/passenger duplex**.

For each focal class, sample-balance all Stage 03 joint-support quantities, including:

```text
varroa_2nt_joint_duplex_fraction
varroa_2nt_reference_fraction_all
varroa_2nt_reference_fraction_recovered
varroa_2nt_reference_fraction_abundance_all
varroa_2nt_reference_fraction_abundance_recovered
```

### Primary joint-support views

For interpretation of geometry conditional on observability, emphasize:

```text
varroa_2nt_reference_fraction_recovered
varroa_2nt_reference_fraction_abundance_recovered
```

because their denominator is restricted to focal references for which at least one passenger was actually recoverable.

`*_fraction_all` quantities remain useful but combine passenger recoverability with geometry support and must be labelled accordingly.

`varroa_2nt_joint_duplex_fraction` remains a distinct duplex-level quantity; it must never be substituted for a focal-reference support fraction.

---

## 04.7 Direct paired 23-versus-24 comparisons

23-versus-24 comparisons are made **within matched sample-virus units** before biological aggregation.

For a fixed focal strand and fixed metric `M`:

```text
paired_delta_24_minus_23(M)
    = M(24 nt) - M(23 nt)
```

The canonical comparison is then:

```text
same sample-virus paired difference
→ median across viruses within sample
→ median across samples
→ sample-clustered bootstrap 95% CI
```

Apply this to the pre-specified joint-support metrics and, where directly meaningful, to the official stepRNA log-ratio for a fixed end and signed distance.

Required focal-strand comparisons are separate:

```text
23S versus 24S
23AS versus 24AS
```

Do not compare `23S` directly with `24AS` as if they were the same focal-strand analysis.

### Design-relevant emphasis

The `23AS versus 24AS` comparison is the primary later-design view because the eventual vdCHIBIN guide strand is antisense to the target transcript. The sense comparison remains part of the biological interpretation.

---

## 04.8 Secondary pair-balanced and virus-balanced summaries

The canonical inference is sample-balanced.

Retain two explicitly secondary robustness views where estimable:

### Pair-balanced

```text
median across sample-virus units
```

This is useful for regression against historical pair-balanced outputs but gives samples with more eligible viruses more influence.

### Virus-balanced

For a fixed metric:

```text
median across samples within biological_virus
→ median across biological viruses
```

A virus-balanced bootstrap resamples biological-virus identities and then contributing observations within selected viruses when implemented.

This view asks whether a conclusion is broadly shared across virus taxa rather than driven by a repeatedly observed virus. It is a sensitivity analysis, not the primary estimator.

All three aggregation labels must remain explicit in output tables.

---

## 04.9 Historical custom `Δ_Dicer` is regression-only

The historical project reported a custom `Δ_Dicer` statistic from the pre-canonical pipeline. It is **not** an official stepRNA metric and is not part of the primary canonical Stage 04 inference.

Do not reconstruct it approximately from the new stepRNA outputs.

Historical regression is allowed only if the **exact archived v1.4.0 implementation and configuration** can be located and executed/read reproducibly. In that case:

- preserve its original target geometry/distance definition;
- preserve its original comparison-distance set;
- preserve its original weighting and permutation scheme;
- preserve its original random seed/replicate count if recorded;
- use Phipson-Smyth-style non-zero Monte Carlo permutation P-values if that is what the archived implementation specifies;
- write outputs under a clearly labelled `historical_regression/` path;
- compare the historical question with the new official stepRNA result without treating numerical agreement as required.

If the exact historical implementation/configuration cannot be recovered, record:

```text
historical_delta_dicer_status = not_reproduced_exact_definition_unavailable
```

and do **not** invent a replacement comparison-distance set.

The historical result may be scientifically reconcilable with Stage 03 even if `(+2,-2)` is not the dominant raw geometry: enrichment above a custom null/comparison spectrum and absolute prevalence are different questions.

---

# 04B — Geometry-conditioned terminal sequence features

## 04.10 Question

Do focal RNAs that support the **pre-specified** Varroa `(+2,-2)` geometry carry terminal nucleotide features beyond those already present in the general Stage 02 population?

Because `(+2,-2)` was pre-specified before Stage 03 results were inspected, it remains the primary geometry-conditioned subset. Do not replace it with blunt geometry or the marginal 3′ `-2` component after seeing the Stage 03 spectrum.

The term **geometry-conditioned** is preferred to **Dicer-conditioned** in canonical metric names because the sequence subset is defined computationally by recovered duplex geometry, not by direct biochemical assignment to Dicer.

---

## 04.11 Geometry-supporting and passenger-recovered focal subsets

Within each:

```text
sample
× analysis_unit
× focal_length {23,24}
× focal_strand {sense,antisense}
× weighting_mode {unique_sequence,abundance}
```

define:

### All focal subset

Every eligible focal reference in the Stage 03 focal manifest.

### Passenger-recovered subset

Every focal reference with at least one recovered passenger, regardless of geometry.

### Joint-geometry subset

Every focal reference with at least one reconstructed passenger duplex satisfying:

```text
5p = +2
3p = -2
```

The joint-geometry subset must be a subset of the passenger-recovered subset. Any violation is a Stage 04 QC failure.

### Weighting

Unique-sequence mode:

```text
one unit per distinct focal reference sequence
```

Abundance mode:

```text
weight each focal reference by focal_abundance
```

Do not weight by the number of recovered passenger alignments.

---

## 04.12 Terminal positions

Use exactly the Stage 02 physical RNA terminal definitions:

```text
5p1 = sequence[0]
5p2 = sequence[1]
3p2 = sequence[-2]
3p1 = sequence[-1]
```

Canonical alphabet remains:

```text
A C G T
```

with `T` corresponding biologically to RNA `U`.

Do not reverse-complement an observed focal sequence before terminal extraction.

---

## 04.13 Reuse the Stage 02 matched viral background

Do not independently recreate a new expected terminal background in Stage 04.

For a matching:

```text
sample
× analysis_unit
× focal_length
× focal_strand
× terminal_position
× nucleotide
```

reuse the canonical Stage 02 `expected_fraction` from:

```text
terminal_expected_by_pair.tsv
```

This guarantees that general and geometry-conditioned enrichments are compared against the same sample-, virus-, length-, strand-, and position-matched viral sequence opportunity.

If a required Stage 02 expected value is absent or undefined, the Stage 04 enrichment/contrast is `NA` and the missing match is reported in QC.

---

## 04.14 Absolute geometry-conditioned enrichment

For terminal feature `f` within one sample-virus unit:

```text
joint_observed_fraction(f)
    = weighted frequency of f among joint-geometry-supporting focal references
```

```text
E_joint_absolute(f)
    = joint_observed_fraction(f)
      / matched Stage 02 expected_fraction(f)
```

Interpretation:

> How enriched is feature `f` among focal sequences supporting the pre-specified joint geometry relative to matched viral sequence opportunity?

This is an **absolute geometry-conditioned enrichment**, not yet evidence that the feature is specifically associated with geometry rather than with passenger recoverability or with the overall Stage 02 population.

If the joint-supporting subset is empty, report `NA` rather than zero.

If the expected fraction is zero or undefined, report `NA`.

No pseudocount is introduced.

---

## 04.15 General-population contrast

Let:

```text
E_all(f)
```

be the matching Stage 02 pair-level enrichment for the same sample-virus × length × strand × weighting mode × terminal feature.

Define:

```text
joint_vs_all_log2_contrast(f)
    = log2(E_joint_absolute(f) / E_all(f))
```

Interpretation:

```text
0  = joint-geometry subset resembles the overall observed small-RNA enrichment
>0 = feature is relatively more enriched in the joint-geometry subset
<0 = feature is relatively less enriched in the joint-geometry subset
```

If either enrichment is non-positive or undefined, report `NA`; do not add a pseudocount solely to force a finite logarithm.

---

## 04.16 Passenger-recovered control contrast

Passenger recovery is an ascertainment step: geometry can only be classified for focal sequences whose complementary passenger was recovered. Therefore Stage 04 adds a control defined **before inspecting geometry-conditioned terminal-nucleotide results**.

For feature `f`:

```text
recovered_observed_fraction(f)
    = weighted frequency of f among all passenger-recovered focal references
```

```text
E_recovered_absolute(f)
    = recovered_observed_fraction(f)
      / matched Stage 02 expected_fraction(f)
```

Then:

```text
joint_vs_recovered_log2_contrast(f)
    = log2(E_joint_absolute(f) / E_recovered_absolute(f))
```

Because the same matched viral background appears in numerator and denominator, this is algebraically equivalent (when all quantities are defined) to:

```text
log2(joint_observed_fraction / recovered_observed_fraction)
```

Interpretation:

> Among focal sequences for which passenger recovery was possible in the observed library, is terminal feature `f` relatively more or less common in those supporting the pre-specified joint geometry?

This is the preferred **geometry-specific sequence contrast** because it reduces, but does not eliminate, confounding by passenger recoverability.

No pseudocount is added. Zero or undefined required quantities produce `NA`.

---

## 04.17 Sample-aware aggregation of geometry-conditioned sequence features

For each fixed:

```text
focal_length
focal_strand
weighting_mode
terminal_position
nucleotide
metric
```

calculate the metric first within each sample-virus unit.

Canonical aggregation then uses:

```text
median across viruses within biological sample
→ median across samples
→ sample-clustered percentile bootstrap 95% CI
```

Report:

- sample-balanced median;
- CI lower/upper;
- number of contributing samples;
- number of contributing sample-virus units;
- number of undefined pair-level values;
- pair-balanced median as a secondary descriptive comparator.

For later candidate development, the primary geometry-conditioned sequence view is:

```text
focal_strand = antisense
```

Sense-strand sequence results are retained but not used as the default future guide-scoring reference.

---

## 04.18 Redundancy and concordance analysis

Stage 04 does **not** assume that geometry-conditioned information deserves a new candidate-scoring dimension.

For the 16 terminal features (`4 positions × 4 nucleotides`), report descriptive Spearman correlations separately by focal length and weighting mode.

### General-versus-joint enrichment correlation

```text
rho_joint_vs_general
    = Spearman correlation across the 16 matched features between:
      sample-balanced E_joint_absolute
      and
      sample-balanced Stage 02 E_all
```

A high positive correlation means the absolute joint-geometry enrichment landscape largely tracks the general terminal-enrichment landscape.

### Abundance-versus-unique contrast concordance

```text
rho_joint_contrast_abundance_vs_unique
    = Spearman correlation across the 16 matched features between
      sample-balanced joint_vs_recovered_log2_contrast
      in abundance and unique-sequence modes
```

This assesses whether any geometry-specific sequence pattern is broadly shared between molecular accumulation and sequence diversity.

These correlations are descriptive; no automatic mechanistic conclusion follows from a high or low rho.

---

## 04.19 No automatic `D` score in Stage 04

Stage 04 does **not** create:

```text
Dicer_score
Dicer_compatibility_score
D
```

or any candidate-window weight.

Instead it produces effect estimates, uncertainty, recoverability controls, and redundancy summaries.

For the currently validated dataset, Stage 04 does **not** justify carrying a separate geometry/Dicer-derived sequence metric into vdCHIBIN ranking. The `(+2,-2)` geometry is a minority same-duplex feature, its geometry-conditioned terminal sequence effects are not consistently strong across weighting modes, and they overlap with the general Stage 02 terminal-enrichment landscape.

A future analysis may revisit this only if independent evidence establishes a sufficiently estimable, reproducible, biologically interpretable, and non-redundant geometry-derived feature. Any such change must be documented prospectively and must not be retroactively encoded into Stage 04.

---

## 04.20 Exploratory marginal 3′ `-2` result

Because Stage 03 showed notable support for the marginal 3′ `-2` distance, Stage 04 may report its sample-balanced official stepRNA log-ratio as part of the **already-required full spectrum**.

However:

- `3p = -2` is not promoted to the primary geometry-conditioned sequence subset;
- no new candidate metric is built from it in Stage 04;
- any dedicated sequence analysis conditioned on `3p = -2` would be post-hoc/exploratory and requires a separately named future analysis.

---

## 04.21 QC

Stage 04 QC must report at least:

- Stage 03 runs represented;
- samples and sample-virus units represented;
- all four focal classes represented;
- duplicate/missing fixed `sample × analysis_unit × focal_length × focal_strand` keys;
- complete same-duplex joint-spectrum accounting for every recovered run;
- exact agreement between summed same-duplex counts and total recovered duplexes;
- per-run same-duplex fraction sum equal to 1 within numerical precision;
- exact regression agreement between full-spectrum `(+2,-2)` counts and the pre-existing `(+2,-2)` joint output;
- missing signed-distance rows required for matched 23/24 comparisons;
- missing/undefined official log-ratios;
- joint-support references absent from the focal manifest;
- joint-support references not present in the passenger-recovered subset;
- focal-abundance mismatches between Stage 03 manifests and Stage 04 calculations;
- missing Stage 02 expected-background matches;
- missing Stage 02 general-enrichment matches;
- terminal-frequency sums across A/C/G/T for non-empty subsets;
- empty joint-support subsets;
- empty passenger-recovered subsets;
- undefined log2 contrasts;
- number of samples contributing to each canonical summary;
- bootstrap requested/valid replicate counts;
- historical-regression status.

A legitimately empty geometry-supporting subset is not a pipeline failure; its derived ratios are `NA` and the condition is reported.

Identifier mismatches, impossible subset relationships, inconsistent abundances, or missing required upstream keys are failures.

---

## 04.22 Outputs

```text
results/04_duplex_geometry/
│
├── qc/
│   ├── stage04_accounting.tsv
│   └── stage04_joint_geometry_spectrum_accounting.tsv
│
├── population/
│   ├── full_spectrum_by_sample.tsv
│   ├── full_spectrum_across_dataset.tsv
│   ├── passenger_recovery_across_dataset.tsv
│   ├── joint_geometry_by_sample.tsv
│   ├── joint_geometry_across_dataset.tsv
│   ├── joint_geometry_spectrum_by_sample.tsv
│   ├── joint_geometry_spectrum_across_dataset.tsv
│   ├── joint_geometry_mode_by_pair.tsv
│   └── joint_geometry_spectrum_summary.tsv
│
├── comparisons/
│   └── paired_23_vs_24.tsv
│
├── sequence_features/
│   ├── geometry_terminal_by_pair.tsv
│   ├── geometry_terminal_by_sample.tsv
│   ├── geometry_terminal_across_dataset.tsv
│   ├── geometry_specific_contrasts.tsv
│   └── redundancy.tsv
│
└── historical_regression/
    └── <only if exact archived v1.4.0 definition is recoverable>
```

No publication figures are required for the initial Stage 04 computational implementation. Figures can be generated after the numerical results have been reviewed.

---

## 04.23 Interpretation limits

Stage 04 may support statements such as:

- a particular signed end-distance feature is reproducibly enriched across biological samples;
- marginal distance `0` is or is not prominent at individual duplex ends;
- fully blunt `(0,0)` duplexes are a minority or substantial same-duplex feature;
- the pre-specified `(+2,-2)` geometry is a minority or substantial same-duplex feature;
- 23 and 24 nt differ or resemble one another for a specified geometry metric;
- geometry-supporting focal RNAs do or do not show terminal sequence preferences beyond the general and passenger-recovered populations.

Stage 04 must **not** by itself claim:

```text
all 23-mers are Dicer products
all 24-mers are secondary/RdRP products
marginal distance 0 means most duplexes are fully blunt
blunt duplexes prove non-Dicer processing
(+2,-2) proves a specific Varroa Dicer enzyme
the observed geometry distribution is sufficient to classify Dicer versus RdRP origin
3′ -2 alone is a validated candidate-design rule
geometry-conditioned enrichment is automatically an independent D score
```

The safest mechanistic language remains **consistent with**, **associated with**, or **supports**, with alternative explanations retained.

---


# 05 — Viral spatial/transitivity-consistency analysis

## 05.1 Biological question

Are local 23-nt viral small-RNA hotspots associated with a reproducible downstream change in the 24-nt antisense population, either in absolute spatial directionality or in antisense 23:24 length composition?

The analysis asks whether such a spatial relationship is **consistent with** amplification/transitivity-associated biology. It does **not** assume in advance that:

```text
23 nt = primary/Dicer product
24 nt = secondary/RdRP product
```

Stages 03–04 showed that duplex geometry is insufficient to make that pathway assignment, so Stage 05 treats 23 nt and 24 nt first as empirically distinct length populations.

This is an observational analysis of natural viral infections. Viral replication, replication intermediates, subgenomic RNA production, local sequence/mappability, RNA stability, and library ascertainment can themselves create spatial structure. A positive Stage 05 result therefore cannot by itself prove RdRP-mediated transitivity or identify the nuclease/polymerase responsible.

### Analysis-only rule

Stage 05 is performed first as biological/spatial analysis only.

It does **not**:

- create a per-window vdCHIBIN transitivity score;
- change candidate-ranking weights;
- label a candidate as “primary-like” or “secondary-like”;
- promote a spatial endpoint into Stage 08/09 design scoring;
- infer a universal propagation distance.

Any later use of Stage 05 in construct-level design interpretation would require a separate, prospectively documented decision after the Stage 05 results are reviewed.

## 05.2 Historical v1.4.1 provenance status

The historical strengthened-transitivity method is retained as a **specification reconstruction**, not an exact source-code replication.

During canonical Stage 05 implementation, the named original v1.4.1 source/result packages could not be located in either:

- the canonical repository; or
- the frozen read-only legacy core.

Accordingly record:

```text
historical_source_package_status = unavailable
historical_rng_stream_status = unavailable
historical_raw_p_checkpoint_status = unavailable
historical_effect_size_regression = PASS
historical_permutation_regression = NOT_EXACTLY_REPRODUCED
```

The historical observed effect-size checkpoints were reproduced within the documented approximate tolerance after coordinate logic was corrected.

The historical permutation/BH values were **not** forced to match. The implementation must not alter seed, iteration order, shift rules, endpoints, anchor definitions, or weighting merely to reproduce archived Monte Carlo values when the exact historical RNG provenance is unavailable.

For backwards-compatible output paths, a directory named:

```text
historical_v1.4.1_replication/
```

may be retained, but its provenance/manifest must explicitly state that the contents are a **specification reconstruction**.
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

## 05.11 23-nt spatial anchor scores

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

These scores identify reproducible local 23-nt hotspots for spatial analysis.

They do **not** imply that the hotspot is biologically primary, Dicer-derived, or upstream of 24-nt production. Historical v1.4.1 terminology referring to “primary-like” anchors is retained only inside the historical replication/provenance where needed for exact reproducibility.

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

Strand-controlled directionality contrast:

```text
antisense_specific_directionality = D_24AS - D_24S
```

A positive value means 24-AS is more downstream-biased than the 24-S control track.

This is a spatial endpoint. It is not, by itself, a primary-versus-secondary or Dicer-versus-RdRP classifier.

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

The use of a circular shift is a **statistical randomization device**. It does not assert that the analysed viral reference is biologically circular. The shift preserves within-track spatial autocorrelation while disrupting registration of the 24-nt tracks relative to the fixed 23-nt anchor locations.

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

For the primary sample-balanced observational analysis, define one pre-specified inferential family per biological endpoint across all predefined combinations:

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

## 05.24A Validated canonical Stage 05 result

Canonical Stage 05 completed successfully with:

```text
targeted tests                = 21/21 PASS
dry-run                       = Stage 05 job only
runtime                       = 193.757 s
QC                            = 11 PASS, 1 WARN, 0 FAIL, 2 INFO
samples                       = 14
eligible +ssRNA units         = 19
biological viruses            = 3
metadata conflicts            = 0
strand mismatches             = 0
weight-normalization deviation= 0
```

The three contributing viruses are:

```text
BMLV
VDV-5
VDV-9
```

The single QC warning records two unit/weighting/anchor combinations with fewer than the required number of retained anchors. This is an estimability warning, not a coordinate/strand/normalization failure.

### Primary endpoint 1 — antisense 23:24 composition

The sample-balanced `delta_F24_AS` result shows the strongest and most reproducible positive spatial association at 250–500 nt.

Key canonical estimates:

| Weighting | Anchor | Window | Estimate | 95% bootstrap CI | Raw permutation P | BH P |
|---|---|---:|---:|---:|---:|---:|
| abundance | balanced23 | 250 | +0.006185 | [-0.000957, +0.015907] | 0.024995 | 0.037493 |
| abundance | balanced23 | 500 | +0.016160 | [+0.007511, +0.025209] | 0.002000 | 0.007199 |
| abundance | combined23 | 250 | +0.005753 | [+0.001319, +0.011406] | 0.002400 | 0.007199 |
| abundance | combined23 | 500 | +0.013266 | [+0.006500, +0.019799] | 0.002999 | 0.007199 |
| unique_sequence | balanced23 | 250 | +0.003698 | [+0.002126, +0.006149] | 0.006999 | 0.011998 |
| unique_sequence | balanced23 | 500 | +0.003344 | [+0.000876, +0.006022] | 0.000400 | 0.004799 |
| unique_sequence | combined23 | 250 | +0.001928 | [-0.001638, +0.003904] | 0.005199 | 0.010398 |
| unique_sequence | combined23 | 500 | +0.002153 | [+0.000332, +0.005789] | 0.001000 | 0.005999 |

All eight pre-specified 250/500-nt combinations are positive and BH-significant in the spatial permutation analysis.

The 100-nt analyses do not show a consistent BH-significant positive shift.

Interpretation:

> downstream regions around 23-nt spatial hotspots contain a modestly larger **fraction** of 24-nt antisense molecules within the combined 23AS+24AS population, particularly over 250–500 nt.

The effect is a change in composition, not a direct percentage increase in total small-RNA abundance.

### Primary endpoint 2 — antisense-specific absolute directionality

No pre-specified combination shows evidence for **positive** antisense-specific 24-nt downstream directionality after BH correction:

```text
BH-adjusted P = 1.0 for all 12 combinations
```

Several point estimates are close to zero or negative. Some bootstrap intervals for unique-sequence analyses are entirely negative; because the pre-specified permutation alternative tests for a positive downstream effect, these observations do not support the hypothesised positive directionality endpoint.

Therefore the canonical result is:

> evidence for a downstream **23:24 antisense compositional shift**, but no evidence for an absolute antisense-specific downstream wave of 24-nt signal relative to the 24S control.

### Leave-one-virus-out robustness

For all eight BH-significant 250/500-nt `delta_F24_AS` results, leave-one-virus-out estimates remained positive.

Observed leave-one-virus-out ranges:

| Weighting | Anchor | Window | LOO estimate range |
|---|---|---:|---:|
| abundance | balanced23 | 250 | +0.004906 to +0.006669 |
| abundance | balanced23 | 500 | +0.012925 to +0.016160 |
| abundance | combined23 | 250 | +0.004996 to +0.010006 |
| abundance | combined23 | 500 | +0.012620 to +0.016646 |
| unique_sequence | balanced23 | 250 | +0.002397 to +0.004277 |
| unique_sequence | balanced23 | 500 | +0.001491 to +0.003947 |
| unique_sequence | combined23 | 250 | +0.000167 to +0.001928 |
| unique_sequence | combined23 | 500 | +0.001388 to +0.002775 |

This shows that the positive compositional shift is not solely dependent on any one of BMLV, VDV-5, or VDV-9.

### Historical regression status

Observed historical effect sizes match the archived approximate checkpoints.

Exact historical permutation/BH reproduction is unresolved because the original historical source/RNG stream and raw P checkpoints are unavailable.

The archived versus reconstructed BH values for the historical unique-sequence × balanced23 `delta_F24_AS` branch were:

| Window | Archived BH | Reconstructed BH |
|---|---:|---:|
| 100 nt | 0.471706 | 0.366927 |
| 250 nt | 0.018596 | 0.016497 |
| 500 nt | 0.000600 | 0.002400 |

These values must be retained as a provenance comparison only. The canonical sample-balanced inference is the inferential result used for biological interpretation.

---

## 05.25 Interpretation rule

The viral analysis may support statements such as:

- 23- and 24-nt spatial association;
- downstream/upstream asymmetry around predefined 23-nt hotspots;
- an antisense-specific 24-nt directional effect if present;
- a downstream shift in antisense 23/24 length composition;
- a spatial pattern **consistent with** amplification/transitivity-associated biology.

The two primary endpoints must remain distinct:

- positive `antisense_specific_directionality` supports an **absolute antisense-specific downstream bias** of the 24-nt population relative to the 24S control;
- positive `delta_F24_AS` supports a **relative downstream shift in antisense 23:24 composition toward 24 nt**.

A positive `delta_F24_AS` with little or no positive `antisense_specific_directionality` therefore means:

> the downstream antisense population is relatively more 24-nt rich,

not:

> there is a demonstrated downstream wave of newly generated 24-nt siRNAs.

Stage 05 must not by itself establish:

- that a 23-nt hotspot is a primary/Dicer cleavage site;
- that every 23-mer is primary;
- that every 24-mer is secondary;
- that a 23→24 spatial association proves biochemical precursor→product order;
- that RdRP directly synthesizes 24-mers;
- that 24-mers are Dicer-independent;
- that spatial directionality proves RNAi transitivity rather than viral replication-associated structure;
- a universal propagation distance;
- a specific Dicer/Ago/RdRP paralogue;
- host-mRNA transitivity.

### Carry-forward rule

Stage 05 produces **analysis outputs only**.

No Stage 05 quantity, including:

```text
antisense_specific_directionality
delta_F24_AS
crosscorr_23_to_24
lag_asymmetry_AS_minus_S
```

is automatically converted into a per-window vdCHIBIN ranking metric.

For the currently validated dataset, Stage 05 remains **analysis-only**. Its positive `delta_F24_AS` result is a regional population-level spatial observation and is not carried into Stage 06–09 as an intrinsic score assigned independently to each hypothetical 24-nt vdCHIBIN candidate window.

Any later construct-level use would require a separately justified and prospectively documented rule.

## 05.26 Required outputs

```text
results/05_viral_transitivity/
    coordinate_qc.tsv
    eligible_positive_sense_units.tsv
    historical_v1.4.1_replication/   # specification reconstruction; original source/RNG unavailable
        transitivity_by_pair.tsv
        pair_balanced_results.tsv
        virus_balanced_results.tsv
        leave_one_virus_out.tsv
        cross_correlation.tsv
        regression_check.tsv
    canonical_transitivity_analysis/   # mechanism-neutral canonical spatial analysis
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

# 06 — Generic transcript target preparation and exhaustive candidate enumeration

## 06.1 Purpose

Stage 06 is the first design-facing stage. Its implementation must be **target-agnostic**.

The core Stage 06 engine must work for any supplied mature/spliced mRNA transcript sequence, not only Vd-CHIBIN.

Vd-CHIBIN / `XM_022792159.1` is the **first validated target instance and regression fixture**, not a hard-coded special case.

Stage 06 answers only:

> For each registered transcript target and requested candidate length, what complete target intervals exist, what is the exact sense/mRNA sequence, what is the corresponding 5′→3′ antisense guide sequence, and what optional transcript-coordinate annotation does each interval overlap?

Stage 06 does **not**:

- run ViennaRNA;
- calculate thermodynamic asymmetry;
- join Stage 02 empirical terminal enrichment;
- use Stage 03/04 duplex geometry;
- use Stage 05 transitivity;
- rank candidates;
- assign candidates to Dicer/RdRP pathways;
- select positive or negative controls;
- compare constructs.

---

## 06.2 Biological coordinate object: the mature transcript

The computational target is a **transcript sequence**, not a genomic DNA interval.

For eukaryotic targets, use the mature/spliced transcript sequence in 5′→3′ orientation. Introns are therefore absent from the Stage 06 coordinate system.

Each transcript isoform is a separate target because alternative splicing changes transcript sequence and transcript-coordinate positions.

The generic Stage 06 engine must not infer transcript sequence by concatenating genomic coordinates during routine candidate enumeration. If a genomic annotation is the original source, transcript extraction/validation belongs in a separate reference-preparation step.

This separation prevents genomic strand, exon order, introns, and alternative isoforms from silently contaminating the candidate-coordinate system.

---

## 06.3 Canonical target registry

The canonical Stage 06 input is a tracked manifest:

```text
resources/targets/target_manifest.tsv
```

One row represents one transcript target.

Required columns:

```text
target_id
transcript_id
display_name
organism
molecule_type
fasta_path
fasta_record_id
annotation_path
expected_length_nt
sequence_sha256_uppercase_dna
candidate_lengths_nt
source_database
source_accession_version
```

### `target_id`

Project-stable, filesystem-safe identifier used in candidate IDs.

It must be unique within the manifest and must not be inferred from gene symbol alone.

### `transcript_id`

Identifier for the exact transcript sequence.

Where an external accession is used, preserve the accession **with version** whenever one exists.

Examples may include:

```text
NM_...
XM_...
ENST...
user-defined transcript ID
synthetic transcript ID
```

The engine must not require a RefSeq prefix.

### `fasta_path` and `fasta_record_id`

The FASTA may be a single-record or multi-record file.

`fasta_record_id` identifies the exact record to use.

The normalized sequence must be interpreted as the transcript 5′→3′ sequence.

### `sequence_sha256_uppercase_dna`

SHA-256 is calculated over the normalized uppercase `A/C/G/T` sequence characters only, excluding FASTA headers and whitespace.

The hash locks the exact transcript sequence independently of filename.

### `candidate_lengths_nt`

Comma-separated positive integers, for example:

```text
23,24
```

The generic engine must enumerate whatever lengths are listed for that target.

The current Vd-CHIBIN project instance uses `23,24`, but the Stage 06 code must not hard-code those values.

### `annotation_path`

Optional transcript-coordinate region annotation.

Use:

```text
NA
```

if no region annotation is available.

Candidate enumeration must still work when annotation is unavailable.

### Source fields

`source_database` and `source_accession_version` are provenance fields.

They may be `NA` for user-defined, de novo assembled, or synthetic transcript targets.

---

## 06.4 Generic transcript-coordinate annotation

Stage 06 uses a simple transcript-coordinate annotation table rather than parsing arbitrary genomic GFF3 directly.

If present, the annotation TSV must contain:

```text
transcript_id
region_label
start_1based
end_1based
coordinate_system
```

Coordinates are:

```text
1-based inclusive transcript coordinates
```

`region_label` is not restricted to Vd-CHIBIN-specific labels.

Common protein-coding mRNA labels may include:

```text
5_prime_UTR
CDS
3_prime_UTR
```

but user-defined labels are allowed.

The core enumerator must **not assume that all transcripts contain all three labels**.

Annotation may be:

- complete;
- partial;
- absent.

For the canonical simple-region table, annotated intervals for one transcript must be ordered and non-overlapping.

Gaps are permitted and are treated as `unannotated`.

If the annotation file is `NA`, annotation-derived candidate fields are `NA`.

This design allows the same engine to handle:

- complete RefSeq/GenBank/Ensembl mRNAs;
- predicted transcripts;
- partial transcripts;
- alternative isoforms;
- de novo assembled transcripts;
- synthetic mRNAs/cDNAs.

---

## 06.5 Generic validation

For every target row:

1. locate the FASTA record identified by `fasta_record_id`;
2. remove whitespace and normalize nucleotide case;
3. permit DNA `T` or RNA `U` input, but normalize internally to uppercase DNA alphabet `A/C/G/T`;
4. reject sequences containing unresolved non-ACGT characters unless a future specification explicitly defines ambiguity handling;
5. verify `expected_length_nt`;
6. verify `sequence_sha256_uppercase_dna`;
7. verify all requested candidate lengths are positive integers and do not exceed transcript length;
8. if annotation is present, verify:
   - matching `transcript_id`;
   - coordinates lie within 1..transcript length;
   - start ≤ end;
   - rows are ordered/non-overlapping after sorting;
   - no duplicate region interval rows.

A sequence/hash mismatch is a FAIL.

The core code must not contain:

```text
XM_022792159.1
710
Vd-CHIBIN
330
665
```

as algorithmic constants.

Those values belong only to the current target registry/reference fixture and target-specific regression tests.

---

## 06.6 Exhaustive generic candidate enumeration

For a transcript of length `L` and requested candidate length `w`:

```text
start_1based = 1 ... L - w + 1
end_1based   = start_1based + w - 1
n_candidates = L - w + 1
```

Every complete interval is retained.

No accessibility, sequence-composition, thermodynamic, empirical, geometry, or transitivity filter is applied.

For a target requesting multiple lengths, enumerate each length independently.

For multiple target-manifest rows, enumerate each target independently and concatenate the resulting rows into the canonical Stage 06 table.

---

## 06.7 Generic candidate identifier

Candidate identifier:

```text
TARGET_ID__LENGTHnt__START_END
```

with zero-padded transcript coordinates.

Examples for the current target:

```text
Vd_CHIBIN__23nt__0001_0023
Vd_CHIBIN__24nt__0001_0024
```

Candidate identity is defined by:

```text
target_id
transcript_id
candidate length
transcript coordinates
```

Sequence identity alone is not sufficient because identical k-mers may occur at multiple positions or targets.

---

## 06.8 Sequence orientation

For every candidate:

```text
target_sequence_dna
    = exact mature-transcript/sense slice, 5′→3′

target_sequence_rna
    = target_sequence_dna with T→U, 5′→3′

antisense_guide_sequence_rna
    = reverse complement of target_sequence_rna, 5′→3′
```

The generic guide rule is:

> the guide sequence is antisense to the supplied target transcript.

Therefore:

```text
guide 5′ nt ↔ complement of target interval final nt
guide 3′ nt ↔ complement of target interval first nt
```

This orientation rule is independent of the transcript's genomic strand because Stage 06 operates in transcript coordinates on the supplied 5′→3′ mRNA sequence.

Target accessibility in Stage 07 must use the target/sense transcript orientation, not the reverse complement.

---

## 06.9 Generic annotation of candidates

Required output fields:

```text
start_region
end_region
overlap_regions
crosses_annotation_boundary
annotation_status
```

### If annotation is complete or partial

`start_region` and `end_region` contain the region label covering the relevant coordinate.

If the coordinate is not covered by any annotation row:

```text
unannotated
```

`overlap_regions` is the ordered semicolon-separated list of all region labels and any unannotated segment touched by the candidate.

`crosses_annotation_boundary = TRUE` when the candidate spans more than one region/unannotated segment.

### If annotation is unavailable

Use:

```text
annotation_status = unavailable
start_region = NA
end_region = NA
overlap_regions = NA
crosses_annotation_boundary = NA
```

Candidate enumeration remains valid.

The Stage 06 engine must not discard a candidate because it crosses a region boundary or lacks annotation.

---

## 06.10 Canonical output schema

Primary output:

```text
results/06_targets/target_candidates.tsv
```

Required columns:

```text
target_id
transcript_id
display_name
organism
candidate_id
candidate_length_nt
start_1based
end_1based
target_sequence_dna
target_sequence_rna
antisense_guide_sequence_rna
annotation_status
start_region
end_region
overlap_regions
crosses_annotation_boundary
```

Optional exact filtered exports may be produced per target or per length, but the combined table is the source of truth.

Do not add ranking columns in Stage 06.

---

## 06.11 Generic QC

For every target × candidate-length stratum verify:

```text
observed candidate count = L - w + 1
```

Also verify:

- first interval = `1..w`;
- last interval = `(L-w+1)..L`;
- every interval has length `w`;
- every legal start occurs exactly once;
- no illegal start occurs;
- candidate IDs are globally unique;
- target sequence equals the exact transcript slice;
- RNA conversion is exact;
- antisense guide is the exact reverse complement;
- no candidate is filtered.

Across the complete output verify:

```text
total observed rows
=
Σtargets Σlengths (L - w + 1)
```

Target-specific annotation counts may be tested as **regression fixtures**, not generic algorithm assumptions.

---

## 06.12 Current Vd-CHIBIN regression fixture

The current first target instance is:

```text
target_id      = Vd_CHIBIN
transcript_id  = XM_022792159.1
display_name   = Vd-CHIBIN
organism       = Varroa destructor
molecule_type  = mRNA
length         = 710
candidate lengths = 23,24
```

Locked SHA-256:

```text
4a0d25aa05b269a118ed1b952dca63ccd1c0a7978fc42295faf3bf650e43ea42
```

Transcript-coordinate annotation:

```text
5_prime_UTR  1–329
CDS          330–665
3_prime_UTR  666–710
```

Expected Vd-CHIBIN enumeration:

```text
23 nt = 688 candidates
24 nt = 687 candidates
combined = 1,375
```

Expected Vd-CHIBIN boundary counts:

```text
23 nt:
5_prime_UTR→CDS = 22
CDS→3_prime_UTR = 22

24 nt:
5_prime_UTR→CDS = 23
CDS→3_prime_UTR = 23
```

These values validate the current target fixture only.

They must not appear as generic Stage 06 algorithmic constants.

---

## 06.13 Generic software architecture requirement

The Stage 06 Python implementation should expose target-agnostic functions such as:

```text
load_target_manifest(...)
load_transcript_sequence(...)
load_transcript_regions(...)
validate_target(...)
enumerate_candidates(...)
annotate_candidate(...)
reverse_complement_rna(...)
```

The command-line entry point should accept at minimum:

```text
--target-manifest
--output-root
```

The core enumerator must not import project-specific Vd-CHIBIN constants.

Vd-CHIBIN-specific expected values belong in test fixtures or target metadata.

This separation is required so that adding another mRNA normally means:

1. add/register its transcript FASTA;
2. optionally add its transcript-coordinate annotation;
3. add one row to `target_manifest.tsv`;
4. rerun Stage 06.

No Python code change should normally be required.

---

## 06.14 Outputs

```text
resources/targets/
└── target_manifest.tsv

results/06_targets/
│
├── target_reference_summary.tsv
├── target_candidates.tsv
│
├── qc/
│   └── stage06_accounting.tsv
│
└── provenance/
    └── stage06_manifest.tsv
```

No figures are required.

No upstream viral analysis is rerun.

---

## 06.15 Interpretation and carry-forward

A Stage 06 candidate row means only:

> this exact interval exists in this exact locked transcript sequence and has this exact antisense guide sequence and annotation context.

It does not imply efficacy or mechanism.

Stage 07 and later candidate analyses must preserve:

```text
target_id
transcript_id
candidate_length_nt
```

as explicit strata/identifiers.

This makes the downstream framework reusable across genes, transcript isoforms, organisms, and future target panels.

---


# 07 — Varroa empirical guide-sequence association landscape

## 07.1 Purpose

Stage 07 extends the Stage 02 terminal-nucleotide analysis across the **entire physical 5′→3′ sequence** of viral 23-nt and 24-nt small RNAs.

Its purpose is to ask:

> Which nucleotide-position features are reproducibly over- or under-represented among observed Varroa viral small RNAs relative to the matched viral sequence opportunity, and which features are associated with disproportionate accumulation among the observed sequences?

Stage 07 contains three logically separate components:

1. **literature-specified validation features**
   - adenine at antisense position 10 (`A10`);
   - continuous GC fraction across antisense positions 9–14 (`GC9_14`);

2. **unbiased single-position discovery**
   - A/C/G/T at every physical position of 23-nt and 24-nt RNAs;

3. **fixed-width regional GC discovery**
   - continuous GC fraction in every contiguous 6-nt window;
   - 23 nt: `GC1_6 ... GC18_23`;
   - 24 nt: `GC1_6 ... GC19_24`;
   - each window is reported in both 5′-based and equivalent 3′-relative coordinates.

The regional width is fixed at **6 nt before canonical implementation** because the published `GC9_14` feature spans six nucleotides. Stage 07 v0.16 does not optimize window width.

Stage 07 performs **no motif/k-mer discovery**, no moving motif scan, and no variable-width regional search.

The stage is an empirical small-RNA **representation/accumulation association analysis**, not an efficacy experiment.

It must not describe a Stage 07 association as:

```text
proven AGO loading preference
proven Dicer preference
proven RdRP preference
proven RNAi efficacy feature
expected Vd-CHIBIN knockdown
```

because the underlying libraries are total small-RNA sequencing rather than AGO-IP/RISC-purified libraries or prospective knockdown experiments.

---

## 07.2 Academic rationale

The external hypotheses come from Cedden et al. (2025), who experimentally compared insect siRNA sequence features in *Tribolium castaneum* and reported associations of higher insecticidal efficacy with, among other variables, adenine at guide position 10 and higher GC content across guide positions 9–14.

Stage 07 does **not** assume that these insect-derived efficacy associations transfer to *Varroa destructor*.

Instead, it asks whether the same sequence features are associated with representation or accumulation in the existing Varroa viral small-RNA population.

Unbiased positional nucleotide analysis also has direct methodological precedent. High-throughput Argonaute-loading experiments have compared nucleotide frequencies at each randomized guide position between selected/loaded and input populations to reveal sequence preferences without specifying each nucleotide-position combination in advance.

The matched-background principle is especially important here: apparent positional preferences must be interpreted relative to the nucleotide composition of the sequence opportunities from which the small RNAs could have arisen.

Cedden et al. tested 21-nt siRNAs in *Tribolium castaneum* and found that higher GC content specifically across antisense positions 9–14, rather than whole-guide GC, was associated with higher efficacy. The same study explicitly cautioned that parameters may require modification in more evolutionarily distant arthropods such as mites. Stage 07 therefore treats `GC9_14` as a literature-guided validation feature while using a **fixed 6-nt sliding regional scan** as a project-specific way to ask whether the relevant local-GC region is shifted in the Varroa small-RNA population.

The regional scan is deliberately narrow. It does not establish that any six-nucleotide boundary is mechanistically privileged, and adjacent windows are expected to be strongly correlated because they share five of six positions. A run of neighboring significant windows is therefore interpreted as evidence for a **broad regional composition pattern**, not as several independent biological discoveries.

Single-position and regional-GC analyses are also not statistically or biologically independent: regional GC is mathematically determined by the constituent G/C frequencies. Later feature integration must therefore assess redundancy before any scoring or weighting decision.

---

## 07.3 Pilot-analysis disclosure

A non-canonical exploratory pilot was run before this formal specification.

The pilot:

- used the 21 read-level libraries;
- retained exact, assigned 23/24-nt reads;
- compared abundance-weighted frequencies with frequencies among distinct observed sequences;
- did **not** apply canonical `primary_eligible` filtering;
- did **not** compare against the matched viral-window background;
- did **not** calculate formal uncertainty or multiple-testing-adjusted inference.

The pilot showed:

```text
A10 abundance/unique ratio:
23AS ≈ 0.943
24AS ≈ 0.955

GC9–14 abundance-minus-unique:
23AS ≈ -0.0035
24AS ≈ -0.0104
```

and nominated several exploratory position-specific signals, most notably:

```text
23 nt: A at position 21 = A at L-2
24 nt: A at position 22 = A at L-2
```

These pilot observations are **not canonical evidence**.

In particular:

- `A10` and `GC9_14` remain literature-specified validation features because they were motivated from the dsRIP study before the pilot;
- `A(L-2)` and all other pilot-nominated positions remain exploratory discoveries and receive no privileged statistical treatment;
- the formal Stage 07 analysis must be interpreted with the knowledge that this dataset has already been inspected once.

Accordingly, Stage 07 should emphasize effect sizes, confidence intervals, replication across samples/lengths/weightings, and conservative discovery correction rather than presenting formal P-values as if the dataset were completely untouched.

### 07.3.1 Post-v0.15 regional-GC exploratory disclosure

After the v0.15 single-nucleotide Stage 07 implementation had passed its canonical checks, a second **non-canonical, read-only** exploratory analysis derived fixed 6-nt regional GC values from the existing `positional_by_pair.tsv`.

This exploratory scan did not rerun Stage 07 or read the frozen core. It nominated a broad antisense accumulation-associated GC-depletion/AU-enrichment region on the 3′ side of both lengths, with the strongest shared 3′-aligned windows approximately:

```text
3p5–10:
23AS accumulation ΔGC ≈ -0.01623
24AS accumulation ΔGC ≈ -0.01949

3p6–11:
23AS accumulation ΔGC ≈ -0.01552
24AS accumulation ΔGC ≈ -0.01945
```

The corresponding sense effects were substantially smaller in the exploratory output.

The same exploration also confirmed that the previously observed single-nucleotide `A at 3p3` association cannot by itself explain the strongest regional result, because the strongest `3p5–10` and `3p6–11` windows do not include `3p3`.

These observations are **pilot-nominated**, not confirmatory. In the formal v0.16 regional analysis:

- `3p5–10`, `3p6–11`, and every other non-`GC9_14` six-nucleotide window receive identical exploratory treatment;
- no window is granted a special threshold, direction, or multiplicity exemption because it looked interesting in the pilot;
- the published `GC9_14` feature remains the only literature-specified regional GC hypothesis;
- no alternative window width is searched.

---

## 07.4 Inputs

Use only the same validated frozen-core biological inputs required by Stage 02:

```text
tables/<sample>/<sample>.read_level_features.tsv.gz
results/descriptive/eligibility.tsv
references/consensus/<sample>.<analysis_unit>.final.background_masked.fa
```

Stage 02 outputs may be read **only for regression/QC comparison**.

Stage 06 target-candidate outputs are not computational inputs to Stage 07.

No remapping is performed.

No stepRNA run is performed.

No ViennaRNA run is performed.

No live database access is required.

---

## 07.5 Canonical observed population

Use the same primary population as Stage 02.

Retain observed read-level rows only when:

```text
sample × virus matches a primary_eligible sample × analysis_unit pair
mapping_mode = exact
virus_assignment = assigned
strand ∈ {sense, antisense}
length ∈ {23, 24}
```

Canonical expected scope from the validated frozen core:

```text
20 primary samples
54 primary-eligible sample-virus units
```

If the observed canonical counts differ from the validated Stage 01/02 population, Stage 07 must stop rather than silently redefining eligibility.

Canonical sequence alphabet in TSV outputs remains:

```text
A, C, G, T
```

where `T` is the DNA-alphabet representation of RNA uridine.

---

## 07.6 Physical sequence orientation and position numbering

Every observed small-RNA sequence is interpreted in its own sequenced physical 5′→3′ orientation.

For a sequence of length `L`:

```text
position_5p = 1 ... L
```

where:

```text
position 1 = physical 5′ nucleotide
position L = physical 3′ nucleotide
```

Also record:

```text
position_from_3p = L - position_5p + 1
```

for interpretation.

Thus:

```text
23-mer position 21 = 3p3
24-mer position 22 = 3p3
```

Observed antisense sequences are used directly.

Do **not** reverse-complement observed antisense reads.

Expected antisense sequences are generated by reverse-complementing the matched reference-orientation viral windows, exactly as in Stage 02.

This distinction is mandatory.

---

## 07.7 Observed weighting modes

Retain the two canonical Stage 02 weighting modes.

### Unique-sequence weighting

Within:

```text
sample × virus × length × strand
```

each distinct sequence contributes weight `1`.

### Abundance weighting

Each included read-level observation contributes its canonical numeric `count`, using the same abundance accounting as Stage 02.

Stage 07 must not invent a new abundance-collapsing rule.

The Stage 02 population/weighting implementation should be reused or regression-tested rather than reinterpreted.

---

## 07.8 Matched viral sequence-opportunity background

For every primary-eligible:

```text
sample × virus × length
```

enumerate all fully supported windows independently within each FASTA record of the sample-specific depth-masked viral consensus.

A supported window contains only:

```text
A, C, G, T
```

and never crosses a FASTA-record boundary.

For sense expectation:

```text
expected sense RNA = reference-orientation window
```

For antisense expectation:

```text
expected antisense RNA = reverse_complement(reference-orientation window)
```

Every supported genomic start contributes one opportunity, even when identical sequence strings occur at multiple positions.

This is the same sequence-opportunity background as Stage 02.

It controls for the fact that a nucleotide can only be observed at a guide position if the viral substrate provides that nucleotide in the corresponding sequence opportunity.

It does not control for all biological or technical biases.

---

# 07A — Full positional representation landscape

## 07.9 Positional observed and expected fractions

For every:

```text
sample
× virus
× length {23,24}
× strand {antisense,sense}
× weighting_mode {unique_sequence,abundance}
× position {1...L}
× nucleotide {A,C,G,T}
```

calculate:

```text
observed_fraction
=
observed weight carrying nucleotide b at position p
/
total observed weight
```

and:

```text
expected_fraction
=
number of supported background windows carrying nucleotide b at p
/
number of supported background windows
```

For every valid population and position:

```text
Σ_b observed_fraction(b) = 1
Σ_b expected_fraction(b) = 1
```

within numerical tolerance.

---

## 07.10 Representation enrichment ratio

Define:

```text
representation_enrichment
=
observed_fraction / expected_fraction
```

Interpretation:

```text
= 1  observed at the frequency predicted by viral sequence opportunity
> 1  over-represented
< 1  under-represented
```

If:

```text
expected_fraction = 0
```

then enrichment is `NA`.

No pseudocount is added.

Also calculate the zero-safe paired difference:

```text
representation_delta_fraction
=
observed_fraction - expected_fraction
```

The enrichment ratio is the main interpretable effect-size ratio.

The delta fraction is the main zero-safe quantity for paired inferential testing.

---

## 07.11 Accumulation contrast

For a fixed:

```text
sample × virus × length × strand × position × nucleotide
```

define:

```text
unique_fraction
=
observed_fraction under unique-sequence weighting

abundance_fraction
=
observed_fraction under abundance weighting
```

Then:

```text
accumulation_ratio
=
abundance_fraction / unique_fraction
```

when `unique_fraction > 0`.

Equivalent, when both representation enrichments are finite:

```text
accumulation_ratio
=
representation_enrichment_abundance
/
representation_enrichment_unique
```

because the same matched expected fraction appears in both denominators.

Interpretation:

```text
> 1
sequences carrying this nucleotide-position feature
contribute disproportionately high accumulated abundance
relative to how common the feature is among distinct observed sequences

≈ 1
little abundance-versus-diversity shift

< 1
feature is disproportionately associated with lower accumulated abundance
```

Also calculate:

```text
accumulation_delta_fraction
=
abundance_fraction - unique_fraction
```

for zero-safe paired inference.

If `accumulation_ratio > 0`, optionally report:

```text
log2_accumulation_ratio
=
log2(accumulation_ratio)
```

No pseudocount is used.

Stage 07 must not translate an accumulation ratio into a fold-change in RNAi efficacy.

---

# 07B — Literature-specified validation features

## 07.12 A10

For antisense 23-nt and 24-nt populations, report the complete Stage 07 positional metrics for:

```text
position = 10
nucleotide = A
```

under both unique-sequence and abundance weighting, together with the abundance-versus-unique accumulation contrast.

The relevant external hypothesis is:

> A at antisense position 10 was associated with higher siRNA efficacy in the Tribolium dsRIP study.

The Stage 07 question is narrower:

> Is A10 over-represented and/or disproportionately accumulated in Varroa viral antisense 23/24-nt small-RNA populations relative to the matched viral sequence opportunity?

A negative or null result must be retained.

No `+10` scoring bonus is imported from dsRIP.

---

## 07.13 Continuous GC9–14

For every observed and expected sequence, calculate:

```text
GC9_14_fraction
=
number of G/C nucleotides at positions 9,10,11,12,13,14
/
6
```

For each sample-virus × length × strand calculate:

```text
observed_GC9_14_mean_unique
observed_GC9_14_mean_abundance
expected_GC9_14_mean
```

Then:

```text
GC9_14_delta_unique_vs_expected
=
observed_GC9_14_mean_unique
-
expected_GC9_14_mean
```

```text
GC9_14_delta_abundance_vs_expected
=
observed_GC9_14_mean_abundance
-
expected_GC9_14_mean
```

and:

```text
GC9_14_accumulation_delta
=
observed_GC9_14_mean_abundance
-
observed_GC9_14_mean_unique
```

The primary GC9–14 feature is **continuous**.

Do not threshold at 50% GC in the canonical analysis.

The dsRIP web tool's later scoring threshold is not treated as the underlying biological hypothesis.

---

## 07.14 Validation-family inference

The literature-specified antisense validation family consists of:

```text
A10:
2 lengths ×
{unique representation, abundance representation, accumulation}

GC9_14:
2 lengths ×
{unique representation, abundance representation, accumulation}
```

for:

```text
12 pre-specified validation tests
```

Use two-sided tests.

The formal analysis must not force the effect into the direction reported in Tribolium.

Apply Benjamini–Hochberg FDR correction across these 12 tests.

Because the same dataset was already viewed in the pilot, these tests are described as **formal literature-guided validation within the existing dataset**, not as independent replication on untouched data.

Primary reporting emphasizes:

```text
effect size
95% sample-clustered bootstrap CI
raw P
BH-adjusted P
direction
```

---

# 07C — Unbiased single-position discovery

## 07.15 Discovery search space

Calculate every nucleotide at every physical position.

No nucleotide-position combination is excluded from the output.

For new discovery claims, define the internal-position search space:

```text
positions 3 ... L-2
```

because:

```text
positions 1, 2, L-1, L
```

were already explicitly analysed in Stage 02 and are reserved primarily for regression/previous-evidence interpretation.

Thus the internal discovery families contain:

```text
23 nt:
19 internal positions × 4 bases = 76 hypotheses per endpoint

24 nt:
20 internal positions × 4 bases = 80 hypotheses per endpoint
```

A10 is present in this unbiased landscape even though it is also reported separately as a literature-specified feature.

Pilot-nominated features such as `A at L-2` receive no special correction or pass criterion.

---

## 07.16 Primary discovery scope

Primary biological discovery scope:

```text
strand = antisense
```

because antisense viral small RNAs are in the orientation complementary to viral target RNA and are therefore the most directly relevant population for later guide-oriented design hypotheses.

However:

```text
antisense read
≠
proven AGO-loaded guide
```

The complete sense analysis is retained as a **secondary comparator** using the same positional definitions and matched sense background.

A discovery that appears similarly in sense and antisense populations is interpreted more cautiously than an antisense-dominant signal.

---

## 07.17 Discovery endpoints and multiple testing

For each length separately, analyze three endpoint families:

```text
1. unique-sequence representation delta
2. abundance-weighted representation delta
3. abundance-versus-unique accumulation delta
```

For each base-position hypothesis, obtain a sample-level paired test as defined below.

Because nucleotide-position hypotheses are strongly dependent:

- the four nucleotide fractions at one position are compositionally constrained;
- neighboring sequence positions can also be correlated;

use:

```text
Benjamini–Yekutieli FDR
```

as the **primary discovery correction** within each:

```text
length × endpoint × strand
```

family.

Also report conventional Benjamini–Hochberg adjusted P-values as a sensitivity/comparability column.

A canonical discovery claim requires at minimum:

```text
BY-adjusted P < 0.05
```

plus a finite, directionally interpretable effect size.

Effect-size magnitude and cross-sample consistency remain more important than the adjusted P-value alone.

---

# 07D — Sample-aware inference

## 07.18 Pair → sample → dataset aggregation

The biological sequencing sample remains the top-level replication unit.

For every fixed metric:

1. calculate the metric separately in every primary-eligible sample-virus unit;
2. within each sample, take the median across its eligible viruses;
3. across samples, take the median of sample-level medians.

For representation enrichment:

```text
dataset effect
=
median across samples
of within-sample median pair-level enrichment ratios
```

For representation delta, accumulation ratio/delta, and GC9–14 deltas, use the same hierarchy.

Do not pool all reads, viruses, or sample-virus units into one pseudo-replicate.

---

## 07.19 Sample-clustered uncertainty

For every canonical dataset-level effect size use a sample-clustered percentile bootstrap:

1. resample sample IDs with replacement;
2. retain all eligible viruses belonging to a resampled sample;
3. recompute within-sample medians;
4. recompute the dataset median;
5. use the 2.5th and 97.5th percentiles.

Canonical bootstrap count:

```text
5000 replicates
```

Record:

```text
seed
requested replicates
valid replicates
interval method
```

No read-level bootstrap is substituted for sample-level uncertainty.

---

## 07.20 Paired inferential test

For inferential P-values, use the sample-level paired **difference** rather than the enrichment ratio.

Representation:

```text
sample_delta
=
within-sample median(
    observed_fraction - expected_fraction
)
```

Accumulation:

```text
sample_delta
=
within-sample median(
    abundance_fraction - unique_fraction
)
```

GC9–14 uses the corresponding observed-minus-expected or abundance-minus-unique sample-level delta.

Primary paired test:

```text
two-sided exact sign test
```

across non-zero sample-level deltas.

This tests whether the direction of the paired difference is consistently positive or negative across independent samples while making minimal distributional assumptions.

Record:

```text
n_samples_total
n_samples_nonzero
n_positive
n_negative
raw_p
```

If too few non-zero sample differences are available for a meaningful test, return `NA` and report the estimability limitation rather than substituting another test silently.

---

# 07E — Stage 02 regression and quality control

## 07.21 Exact terminal regression

The full Stage 07 positional engine must reproduce Stage 02 at the four previously analysed terminal positions.

Mapping:

```text
position 1   ↔ 5p1
position 2   ↔ 5p2
position L-1 ↔ 3p2
position L   ↔ 3p1
```

Regression must cover:

```text
23 nt and 24 nt
sense and antisense
unique-sequence and abundance weighting
all A/C/G/T nucleotides
pair-level observed fractions
pair-level expected fractions
pair-level enrichment ratios
sample-balanced aggregate enrichment
```

Expected numerical tolerance:

```text
<= 1e-12
```

for quantities that are algebraically identical.

Any systematic failure is treated as an orientation/background/weighting bug and blocks canonical Stage 07 completion.

---

## 07.22 Accounting QC

Stage 07 must verify:

- primary eligibility matches Stage 02;
- expected primary sample count = 20;
- expected primary-eligible sample-virus units = 54;
- mapping mode is exact;
- virus assignment is assigned;
- only 23/24 nt are included;
- observed sequence length agrees with declared length;
- no unexpected nucleotide alphabet is silently coerced;
- sense/antisense orientation is preserved;
- expected windows never cross FASTA-record boundaries;
- unsupported/N-containing background windows are excluded;
- A/C/G/T observed and expected fractions each sum to 1 at every valid position;
- no pseudocount is used in enrichment calculations;
- all 23 positions and all 24 positions are present where estimable;
- exactly 18 six-nucleotide windows are generated for 23 nt and 19 for 24 nt;
- 5′-based and 3′-relative regional coordinates map to the same physical window;
- direct per-sequence regional GC equals the constituent-position G+C derivation within numerical tolerance;
- `GC9_14` regional calculations reproduce the existing canonical GC9–14 calculations;
- `GC9_14` is not double-counted in the exploratory regional multiple-testing family;
- no window width other than 6 nt is generated;
- Stage 02 terminal regression passes.

---

# 07F — Fixed-width regional GC landscape

## 07.23 Rationale and search space

The canonical regional analysis uses one pre-specified width:

```text
window_width_nt = 6
```

For a guide/read of length `L`, enumerate all contiguous six-nucleotide windows:

```text
start_5p = 1 ... L-5
end_5p   = start_5p + 5
```

Therefore:

```text
23 nt -> 18 windows
24 nt -> 19 windows
```

The same physical window must also be reported in 3′-relative coordinates:

```text
near_3p = L - end_5p + 1
far_3p  = L - start_5p + 1
```

Example:

```text
23 nt positions 14–19 -> 3p5–10
24 nt positions 15–20 -> 3p5–10
```

The 5′ and 3′ coordinate labels describe the **same hypothesis**. They must never be treated as separate tests or separately corrected discoveries.

The canonical v0.16 search does not scan widths other than 6 nt.

---

## 07.24 Per-sequence regional GC fraction

For one observed or expected sequence and one six-nucleotide window:

```text
regional_gc6_fraction
=
number of G/C nucleotides in the six positions
/
6
```

Allowed values:

```text
0, 1/6, 2/6, 3/6, 4/6, 5/6, 1
```

For every:

```text
sample × virus × length × strand × six-nt window
```

calculate:

```text
observed_gc6_mean_unique
observed_gc6_mean_abundance
expected_gc6_mean
```

where:

- `observed_gc6_mean_unique` gives each distinct observed sequence weight 1;
- `observed_gc6_mean_abundance` uses the canonical Stage 02/07 abundance weights;
- `expected_gc6_mean` gives each fully supported matched viral sequence opportunity weight 1.

For fixed-width windows, the mean regional GC fraction is algebraically equal to the mean across the six constituent positional `P(G)+P(C)` frequencies. An implementation may exploit this equivalence only if deterministic tests confirm exact agreement with direct sequence-level calculation.

---

## 07.25 Regional representation and accumulation effects

Define:

```text
regional_gc6_delta_unique_vs_expected
=
observed_gc6_mean_unique
-
expected_gc6_mean
```

```text
regional_gc6_delta_abundance_vs_expected
=
observed_gc6_mean_abundance
-
expected_gc6_mean
```

```text
regional_gc6_accumulation_delta
=
observed_gc6_mean_abundance
-
observed_gc6_mean_unique
```

Interpretation:

### Positive unique/abundance representation delta

Observed small RNAs are more GC-rich in that region than matched viral sequence opportunity predicts.

### Negative unique/abundance representation delta

Observed small RNAs are more GC-poor/AU-rich in that region than matched viral sequence opportunity predicts.

### Positive accumulation delta

Among sequences that are observed, higher-copy products are more GC-rich in that region than the diversity of observed sequences.

### Negative accumulation delta

Among sequences that are observed, higher-copy products are more GC-poor/AU-rich in that region than the diversity of observed sequences.

These are **percentage-point composition differences**, not fold-changes in RNA abundance or RNAi efficacy.

No GC threshold is imposed.

No regional score or bonus is calculated.

---

## 07.26 Regional aggregation and uncertainty

Use exactly the Stage 07 sample-aware hierarchy:

```text
sample-virus regional metric
        ↓
median eligible viruses within sample
        ↓
sample-level regional metric
        ↓
median across samples
```

For every canonical dataset-level regional effect:

- use the same 5000-replicate sample-clustered percentile bootstrap;
- use the same fixed documented seed/provenance requirements;
- use the same two-sided exact sign test on non-zero sample-level deltas.

The sample remains the top-level biological replication unit.

---

## 07.27 Regional multiple-testing families

`GC9_14` is already a literature-specified validation feature in Section 07B.

Therefore:

- `GC9_14` remains present in the complete regional tables;
- its canonical inferential adjustment comes from the **12-test literature-validation BH family**;
- it is **not counted again** as a novel regional discovery.

All other six-nucleotide windows are exploratory.

For each:

```text
length × endpoint × strand
```

regional discovery family, apply:

```text
primary: Benjamini–Yekutieli FDR
sensitivity: Benjamini–Hochberg FDR
```

where endpoint is one of:

```text
unique_representation
abundance_representation
accumulation
```

Family sizes excluding the literature-specified `GC9_14` window are:

```text
23 nt -> 17 exploratory windows per endpoint × strand
24 nt -> 18 exploratory windows per endpoint × strand
```

A novel regional window is statistically supported only when:

```text
regional BY-adjusted P < 0.05
```

with a finite effect and estimable sample-level inference.

Because neighboring windows overlap heavily, several adjacent supported windows must not be counted as independent biological evidence. Interpretation should describe the **shared broad region** and report the neighboring-window pattern.

The pilot-nominated `3p5–10` and `3p6–11` windows receive no special treatment.

---

## 07.28 Antisense primary analysis and sense comparator

As for the single-nucleotide analysis:

```text
primary biological scope = antisense
secondary comparator = sense
```

The complete regional scan is calculated for both strands.

A regional effect that is much stronger or more consistent in antisense than in sense may be described as **antisense-associated**.

It must not be described as a proven AGO-loading feature because total small-RNA libraries do not isolate AGO-bound molecules.

---

## 07.29 Relationship between single-nucleotide and regional evidence

The regional-GC landscape is a complementary summary of the single-position nucleotide landscape, not an independent source of sequence data.

For example:

```text
low GC in 3p5–10
```

can arise from depletion of G/C and/or enrichment of A/T at one or several constituent positions.

Therefore Stage 07 must preserve both views:

```text
single nucleotide × position
regional continuous GC
```

but later integration must explicitly assess:

- whether a regional effect is distributed across several positions or dominated by one position;
- whether a named single-nucleotide feature and a regional-GC feature are strongly correlated;
- whether combining them would double-count the same underlying empirical tendency.

Stage 07 itself does not resolve redundancy by assigning weights.

---

## 07.30 Relationship to duplex-end thermodynamic asymmetry

Regional GC composition is not equivalent to duplex-end thermodynamic stability.

The planned Stage 08 thermodynamic asymmetry feature uses nearest-neighbour stack energies at the guide and passenger 5′ duplex ends.

A six-nucleotide regional GC window:

- discards sequence order;
- may lie partly or wholly inward from the terminal stacks;
- is only a coarse proxy for local duplex stability.

A 3′-side GC-depletion signal could therefore be independent of, partially correlated with, or directionally opposed to the classical preference for a relatively less-stable guide 5′ end than passenger 5′ end.

Stage 07 must not force agreement between these evidence types.

Their relationship is evaluated only after Stage 08 features exist.

---

# 07G — Outputs

## 07.31 Canonical output tree

```text
results/07_empirical_sequence/
│
├── positional_by_pair.tsv
├── positional_by_sample.tsv
├── positional_summary.tsv
│
├── gc9_14_by_pair.tsv
├── gc9_14_by_sample.tsv
├── gc9_14_summary.tsv
│
├── regional_gc6_by_pair.tsv
├── regional_gc6_by_sample.tsv
├── regional_gc6_summary.tsv
├── regional_gc6_discovery.tsv
│
├── literature_validation.tsv
├── discovery_summary.tsv
├── sense_comparator.tsv
│
├── qc/
│   ├── stage07_accounting.tsv
│   └── stage02_terminal_regression.tsv
│
└── provenance/
    └── stage07_manifest.tsv
```

Canonical inference lives in TSV outputs.

Plots may be produced for interpretation, but no plot is a source of truth.

---

## 07.32 Required positional-summary fields

At minimum:

```text
strand
length
position_5p
position_from_3p
nucleotide
weighting_mode

n_samples
sample_balanced_observed_fraction
sample_balanced_expected_fraction
sample_balanced_representation_enrichment
sample_balanced_representation_delta_fraction

bootstrap_ci_low
bootstrap_ci_high

sign_test_n_nonzero
sign_test_n_positive
sign_test_n_negative
raw_p
bh_p
by_p
```

The accumulation summary additionally records:

```text
sample_balanced_unique_fraction
sample_balanced_abundance_fraction
sample_balanced_accumulation_ratio
sample_balanced_accumulation_delta_fraction
log2_accumulation_ratio
```

where estimable.

The regional-GC summary additionally records at minimum:

```text
strand
length
start_5p
end_5p
near_3p
far_3p
region_5p
region_3p
endpoint
n_samples

sample_balanced_observed_gc6
sample_balanced_expected_gc6
sample_balanced_regional_gc6_delta

bootstrap_ci_low
bootstrap_ci_high

sign_test_n_nonzero
sign_test_n_positive
sign_test_n_negative
raw_p

evidence_class
validation_bh_p
regional_bh_p
regional_by_p
```

`evidence_class` distinguishes at least:

```text
literature_validation_gc9_14
exploratory_regional_gc6
```

The `GC9_14` row must not be double-counted as an exploratory regional discovery.

---

# 07H — Interpretation limits

## 07.33 What Stage 07 can support

A robust result can support wording such as:

> nucleotide X at physical position p is reproducibly over-represented among observed Varroa viral antisense 23/24-nt small RNAs relative to matched viral sequence opportunity.

or:

> observed sequences carrying nucleotide X at position p contribute disproportionately high abundance relative to their prevalence among distinct observed sequences.

For regional GC it may support wording such as:

> a broad 3′-relative region is reproducibly GC-depleted/AU-enriched among high-abundance antisense products relative to the diversity of observed sequences.

When several overlapping six-nucleotide windows are supported, the preferred interpretation is a **regional pattern**, not a claim that the numerically strongest window boundary is uniquely mechanistic.

It may also support:

> this association is shared across 23/24 nt,

or:

> this association is stronger in antisense than in the sense comparator,

when the corresponding analyses justify those statements.

---

## 07.34 What Stage 07 cannot establish

Stage 07 does not directly distinguish among:

- Dicer cleavage preference;
- Argonaute loading preference;
- small-RNA stability;
- source-RNA abundance;
- viral replication/transcription structure;
- other processing biology;
- technical library-preparation bias.

Small-RNA library construction can itself introduce sequence-dependent representation bias, particularly through RNA-ligation steps and interactions between adapters and small-RNA sequence/structure.

The matched viral-window background corrects for **viral sequence opportunity**, but it does not remove sequence-dependent library-preparation bias.

Therefore even a highly reproducible Stage 07 positional association is an **empirical small-RNA representation/accumulation feature**, not automatically a molecular selection rule.

Direct AGO/RISC-IP sequencing or prospective wet-lab knockdown experiments would provide stronger evidence for loading or efficacy.

---

## 07.35 Relationship to candidate design

Stage 07 creates **organism-specific empirical sequence evidence** that may later be evaluated on Stage 06 target candidates.

It does not assign a candidate score in Stage 07.

A discovered feature is not automatically converted into points or a fixed weight.

Later integration must distinguish:

```text
literature-supported + Varroa-supported feature

Varroa-discovered association only

literature-supported feature not reproduced in Varroa
```

and must preserve the difference between:

```text
representation/accumulation association
```

and:

```text
experimentally validated knockdown efficacy
```

No motif/k-mer feature is introduced by this stage.

Both single-nucleotide and regional-GC outputs are eligible to be **carried forward as empirical candidate features for later evaluation**, but this does not imply inclusion in the final score.

Examples include:

```text
single-position indicator:
A at guide 3p3

continuous regional composition:
guide GC fraction in a selected 3′-relative six-nucleotide region
```

Any later choice of a named empirical feature must document:

- its Stage 07 evidence class;
- whether it was literature-specified or discovered in this dataset;
- redundancy with neighboring single-position/regional features;
- relationship to Stage 08 thermodynamic/structural features;
- whether it has direct efficacy validation.

No predictive train/test model is part of canonical Stage 07 v1.

A regularized multivariable predictive model may be considered later only if the transparent single-position analysis identifies sufficiently stable signal and a clear prediction target is defined.

---

# 08 — Generic candidate biophysics: target accessibility and duplex-end thermodynamic asymmetry

## 08.1 Purpose

Stage 08 calculates two **raw biophysical candidate features** for every Stage 06 transcript interval:

1. target-site accessibility on the supplied mature transcript; and
2. relative thermodynamic stability of the two 5′ ends of the perfectly complementary guide/passenger duplex.

The implementation must remain **target-agnostic**.

Vd-CHIBIN is the first regression target, not an algorithmic special case.

Stage 08 does **not**:

- join Stage 02 nucleotide-enrichment values;
- use Stage 03/04 duplex geometry;
- use Stage 05 transitivity;
- convert raw features to percentiles;
- calculate A, T, N, S or another aggregate score;
- apply an asymmetry pass/fail gate;
- rank candidates;
- select target regions;
- implement Stage 09.

The historical weighted A/T/N/S scanner is a regression reference only and is not automatically inherited by the canonical pipeline.

---

## 08.2 Inputs

Canonical inputs:

```text
resources/targets/target_manifest.tsv
results/06_targets/target_candidates.tsv
```

plus each registered transcript FASTA referenced by the target manifest.

No Stage 00–05 result table is a computational dependency.

No live database fetch is permitted during canonical Stage 08 execution.

Stage 08 must preserve at least:

```text
target_id
transcript_id
candidate_id
candidate_length_nt
start_1based
end_1based
target_sequence_rna
antisense_guide_sequence_rna
```

from Stage 06.

---

# 08A — Target-site accessibility

## 08.3 Biological interpretation

Accessibility is calculated on the **target transcript / sense mRNA**, not on the antisense guide.

The primary quantity is:

```text
P_unpaired(candidate)
```

defined as the equilibrium probability, under the RNAplfold local partition-function model, that the **entire candidate interval** is simultaneously unpaired.

Higher values mean that the complete target interval is predicted to be more available for intermolecular recognition.

This is a computational secondary-structure prediction, not a direct measurement of in-vivo RNA structure.

Published experimental and computational studies support target accessibility as a determinant of siRNA/RISC target recognition, while also showing that accessibility is only one contributor to silencing efficacy.

---

## 08.4 RNAplfold model

Use ViennaRNA `RNAplfold`.

Canonical software requirement:

```text
ViennaRNA = 2.7.2
```

or an explicitly updated version recorded in provenance and revalidated against the Stage 08 regression tests before becoming canonical.

RNAplfold computes local partition functions and unpaired probabilities using sliding windows.

For a candidate ending at transcript coordinate `e` with length `w`, extract the `_lunp` value corresponding to:

```text
[e-w+1 ... e]
```

at interval length `w`.

The extraction must be tested against hand-checked synthetic output so that end-coordinate versus start-coordinate indexing cannot be reversed.

---

## 08.5 Canonical accessibility parameter profile

The current project profile preserves the previously audited three RNAplfold scales:

Primary:

```text
label = W150_L100_main
W = 150
L = 100
```

Sensitivity:

```text
label = W100_L80_sensitivity
W = 100
L = 80
```

and:

```text
label = W200_L150_sensitivity
W = 200
L = 150
```

Canonical model temperature:

```text
T = 37.0 °C
```

For each target:

```text
u = maximum candidate_length_nt requested for that target
```

Thus the current Vd-CHIBIN 23/24 analysis uses:

```text
u = 24
```

and both 23-nt and 24-nt accessibility values are extracted from the same RNAplfold run for a given W/L parameter set.

The three W/L choices are **project modelling scales**, not universal biological constants. Their purpose is to test whether candidate accessibility is robust to reasonable changes in local folding context.

Stage 09 must not treat agreement across these settings as independent biological replication.

---

## 08.6 Generic handling of transcript length

The Stage 08 engine must not assume a 710-nt transcript.

For requested RNAplfold window `W_req`, maximum base-pair span `L_req`, and transcript length `N`:

```text
W_eff = min(W_req, N)
L_eff = min(L_req, W_eff)
```

Record both requested and effective values.

Require:

```text
u <= W_eff
```

for every enabled parameter set.

If this is not satisfied, fail clearly rather than silently truncating the requested unpaired interval.

---

## 08.7 Accessibility outputs

For each candidate record:

```text
accessibility_p_W150_L100
accessibility_p_W100_L80
accessibility_p_W200_L150
```

and descriptive robustness summaries:

```text
accessibility_p_min
accessibility_p_max
accessibility_p_range
```

where:

```text
accessibility_p_min
    = minimum probability across enabled parameter sets

accessibility_p_max
    = maximum probability across enabled parameter sets

accessibility_p_range
    = accessibility_p_max - accessibility_p_min
```

These remain raw/descriptive values.

Do **not** convert them to percentiles or combine them with fixed weights in Stage 08.

Important cross-length rule:

> the probability that an entire 24-nt interval is unpaired is a stricter event than the corresponding 23-nt event.

Therefore raw 23-nt and 24-nt full-interval probabilities must not automatically be treated as directly comparable efficacy scores.

Any cross-length normalization belongs in Stage 09.

---

# 08B — Duplex-end thermodynamic asymmetry

## 08.8 Biological rationale

For a perfectly complementary small-RNA duplex, the two strands can differ in the local stability of their 5′ ends.

Classical RNAi studies showed that the strand whose 5′ end is less stably paired is often preferentially selected as guide during RISC assembly. Drosophila studies provide direct mechanistic evidence that RISC-loading machinery can sense terminal duplex asymmetry.

However:

```text
favourable terminal asymmetry
≠
guaranteed guide loading
≠
guaranteed knockdown
```

The canonical project has not established a Varroa-specific quantitative loading rule.

Stage 08 therefore calculates **raw asymmetry features** and directional signs only.

---

## 08.9 Duplex definition

For each Stage 06 candidate:

```text
guide
    = antisense_guide_sequence_rna, 5′→3′

passenger
    = target_sequence_rna, 5′→3′
```

These sequences are exact reverse complements and define a perfectly paired RNA/RNA duplex for the purpose of the terminal-stability calculation.

No Dicer overhang is added.

No mismatch is introduced.

No candidate-specific secondary structure is folded for this calculation.

The metric is a local duplex-end stability surrogate, not the free energy of the entire small-RNA duplex.

---

## 08.10 Nearest-neighbour model

Use canonical Watson–Crick RNA/RNA nearest-neighbour stacking free energies from the ViennaRNA default Turner-2004 parameter model at 37 °C.

The implementation must obtain or validate the stacking parameters from the pinned ViennaRNA parameter set; it must not silently substitute unrelated DNA, RNA/DNA or alternative fitted parameters.

The provenance record must include:

```text
ViennaRNA version
parameter-set identity
temperature
```

The nearest-neighbour model is appropriate because RNA duplex stability depends on adjacent base-pair stacks rather than independent base-pair counts.

---

## 08.11 Main terminal asymmetry definition

Primary terminal length:

```text
4 paired nucleotides
```

This contains:

```text
3 nearest-neighbour stacks
```

For the guide 5′ end:

```text
guide_5p_stack_dg_4bp_kcal_mol
    = sum of the 3 Turner nearest-neighbour stack ΔG°37 terms
      across the first 4 paired nucleotides at the guide 5′ end
```

For the passenger 5′ end:

```text
passenger_5p_stack_dg_4bp_kcal_mol
    = analogous sum at the passenger 5′ end
```

Define:

```text
asymmetry_ddg_4bp_kcal_mol
    =
    guide_5p_stack_dg_4bp_kcal_mol
    -
    passenger_5p_stack_dg_4bp_kcal_mol
```

Because more-negative duplex free energy is more stable:

```text
asymmetry_ddg_4bp > 0
```

means:

> the desired antisense guide 5′ end is less stably paired than the passenger 5′ end.

This direction is consistent with the classical strand-selection tendency.

Stage 08 must not call this a pass/fail rule.

---

## 08.12 Five-nucleotide sensitivity feature

Also calculate a pre-specified terminal-length sensitivity feature using:

```text
5 paired nucleotides
= 4 nearest-neighbour stacks
```

with:

```text
guide_5p_stack_dg_5bp_kcal_mol
passenger_5p_stack_dg_5bp_kcal_mol
asymmetry_ddg_5bp_kcal_mol
```

and the same sign convention.

This sensitivity is retained because published Drosophila siRNA analyses have used terminal free-energy windows spanning approximately the first five nucleotides, and because the previous project workflow used a 5-bp sensitivity calculation.

The 4-bp and 5-bp features are correlated model summaries, not independent experiments.

---

## 08.13 Stack-only terminology

The Stage 08 quantities above are deliberately named:

```text
stack_dg
```

because they sum nearest-neighbour stacking terms.

They are **not** labelled as complete duplex formation free energies.

Unless explicitly added in a future specification, Stage 08 does not add:

- duplex initiation terms;
- 3′ DNA/RNA overhang terms;
- concentration-dependent melting terms;
- whole-duplex RNAcofold energies;
- target-opening penalties to the asymmetry calculation.

This explicit naming prevents the historical stack-only metric from being over-interpreted as a complete physical ΔG of duplex formation.

---

## 08.14 Asymmetry robustness annotations

Descriptive fields:

```text
asymmetry_direction_4bp
asymmetry_direction_5bp
asymmetry_direction_consistent
```

where direction is:

```text
guide_5p_less_stable     if ΔΔG > 0
equal_within_numeric_tol if |ΔΔG| <= numeric tolerance
guide_5p_more_stable     if ΔΔG < 0
```

and `asymmetry_direction_consistent` indicates whether the 4-bp and 5-bp calculations support the same non-zero direction.

This is annotation only.

No candidate is filtered in Stage 08.

---

# 08C — Generic outputs, QC and interpretation

## 08.15 Canonical outputs

```text
results/08_candidate_biophysics/
│
├── candidate_accessibility.tsv
├── candidate_thermodynamic_asymmetry.tsv
├── candidate_biophysics.tsv
│
├── qc/
│   └── stage08_accounting.tsv
│
└── provenance/
    ├── stage08_manifest.tsv
    └── rnaplfold_runs.tsv
```

`candidate_biophysics.tsv` is the canonical Stage 08 joined table.

It must contain exactly one row for every Stage 06 candidate.

---

## 08.16 Required canonical joined-table fields

Carry forward:

```text
target_id
transcript_id
candidate_id
candidate_length_nt
start_1based
end_1based
target_sequence_rna
antisense_guide_sequence_rna
```

Accessibility:

```text
accessibility_p_W150_L100
accessibility_p_W100_L80
accessibility_p_W200_L150
accessibility_p_min
accessibility_p_max
accessibility_p_range
```

Thermodynamic asymmetry:

```text
guide_5p_stack_dg_4bp_kcal_mol
passenger_5p_stack_dg_4bp_kcal_mol
asymmetry_ddg_4bp_kcal_mol

guide_5p_stack_dg_5bp_kcal_mol
passenger_5p_stack_dg_5bp_kcal_mol
asymmetry_ddg_5bp_kcal_mol

asymmetry_direction_4bp
asymmetry_direction_5bp
asymmetry_direction_consistent
```

No Stage 08 column may be called:

```text
A
T
N
S
score
rank
pass
fail
primary_score
secondary_score
```

---

## 08.17 Generic QC

### Accounting

For every target × candidate-length stratum:

```text
Stage 08 row count = Stage 06 row count
```

and globally:

```text
Stage 08 canonical rows = Stage 06 canonical rows
```

No candidate may disappear because of low accessibility or unfavourable asymmetry.

### Accessibility

Require:

- all probabilities finite;
- `0 <= P_unpaired <= 1`;
- exact Stage 06 coordinate match;
- candidate length equals extracted RNAplfold interval length;
- extraction uses candidate end coordinate and requested interval length correctly;
- sensitivity summary arithmetic is exact;
- one RNAplfold run per target × parameter set, not per candidate;
- requested/effective W/L/u values recorded.

### Thermodynamic asymmetry

Require:

- guide remains exact reverse complement of target;
- 4-bp calculation uses exactly 3 terminal stacks;
- 5-bp calculation uses exactly 4 terminal stacks;
- guide 5′ and passenger 5′ ends are not swapped;
- stored `ΔΔG = guide - passenger`;
- sign annotation matches stored ΔΔG;
- hand-calculated regression examples agree with implementation;
- no infeasible-energy sentinel such as `100000` is accepted.

### Genericity

Tests must include at least one non-Vd-CHIBIN synthetic transcript and non-23/24 candidate length for the reusable helper functions.

Target-specific regression tests may additionally use Vd-CHIBIN.

---

## 08.18 Vd-CHIBIN regression requirements

For the current first target:

```text
target_id = Vd_CHIBIN
transcript_id = XM_022792159.1
Stage 06 rows = 1,375
```

Expected Stage 08 rows:

```text
23 nt = 688
24 nt = 687
total = 1,375
```

Accessibility regression should compare the canonical W150/L100, W100/L80 and W200/L150 raw values against the previously audited Vd-CHIBIN accessibility outputs when those historical files are supplied as a regression fixture.

Thermodynamic regression should compare the new 4-bp/5-bp stack-sum calculation against the previously corrected Vd-CHIBIN thermodynamic table when that historical table is supplied as a regression fixture.

Historical files are **regression checks only**, not Stage 08 computational inputs.

If the historical regression fixture is absent, canonical Stage 08 must still run from Stage 06 + the registered transcript FASTA.

---

## 08.19 Interpretation limitations

Accessibility limitations:

- RNAplfold predicts a thermodynamic secondary-structure ensemble;
- it does not directly measure in-vivo structure;
- RNA-binding proteins, translation, modifications, cellular localization and non-equilibrium folding can alter real accessibility;
- W/L choices define a model of local structural context.

Thermodynamic-asymmetry limitations:

- the feature is based on a perfectly complementary duplex;
- the stack-only calculation is a local stability surrogate;
- actual small-RNA termini, overhangs, 5′ phosphorylation, Argonaute identity and loading cofactors can affect strand selection;
- classical asymmetry rules are supported strongly in model systems such as Drosophila, but no quantitative Varroa-specific loading function is assumed.

Therefore Stage 08 results should be described as:

```text
predicted target accessibility
and
predicted duplex-end thermodynamic asymmetry
```

not as direct measurements of RNAi efficacy.

---

## 08.20 Carry-forward to Stage 09

Stage 09 may integrate Stage 08 raw features with empirical Stage 02 terminal-nucleotide evidence.

Before defining Stage 09, the project must decide explicitly:

- within-target normalization;
- within-length versus cross-length comparison;
- whether accessibility robustness should use main/min/range or another rule;
- whether asymmetry sign should be a gate, a continuous feature, or both;
- whether any weights are justified.

The historical:

```text
A = 0.70 main + 0.30 robustness
T = 0.70 main + 0.30 robustness
S = 0.60A + 0.30T + 0.10N
```

must **not** be imported automatically.

No arbitrary ranking weights are defined in Stage 08.

---

## 6. Reproducibility requirements

All new outputs are written under `results/` in the canonical repository.

Configuration must record at minimum:

```text
target_lengths = [23, 24]
stage06_target_manifest = "resources/targets/target_manifest.tsv"
stage06_coordinate_system = "1-based inclusive transcript coordinates"

stage08_viennarna_version = "2.7.2"
stage08_temperature_c = 37.0
stage08_accessibility_parameter_sets = [
    ["W150_L100_main", 150, 100],
    ["W100_L80_sensitivity", 100, 80],
    ["W200_L150_sensitivity", 200, 150]
]
stage08_asymmetry_terminal_lengths_nt = [4, 5]
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

For Stage 05, provenance must additionally record:

```text
historical_source_package_status
historical_rng_stream_status
historical_raw_p_checkpoint_status
historical_effect_size_regression
historical_permutation_regression
```

Relevant software versions include Python, Snakemake, stepRNA, Bowtie2 used by stepRNA, pysam, NumPy, pandas, SciPy, and plotting/statistical packages actually used.

---

## 7. Required deterministic tests

### Stage 02

- exact 5p1/5p2/3p2/3p1 extraction;
- observed antisense sequences are not reverse-complemented;
- expected antisense windows are reverse-complemented correctly;
- 23-nt and 24-nt window enumeration is inclusive and never crosses FASTA-record boundaries;
- windows containing `N`/non-ACGT background sequence are excluded;
- abundance mode uses read-level `count`, not row count;
- unique-sequence mode uses the Stage 01 strand-specific sequence identity;
- observed and expected A/C/G/T frequencies each sum to 1 when defined;
- observed zero denominator → `NA`;
- expected zero frequency → enrichment `NA`;
- observed zero with positive expected frequency → enrichment 0;
- strand-weighted combined expectation uses the observed mixture for the matching weighting mode;
- sample-level median is taken across pair-level enrichment ratios, not as a ratio of separately aggregated frequencies;
- sample-balanced median and sample-clustered bootstrap are reproducible with fixed seed;
- pooled-abundance observed and abundance-matched expected fractions use the documented molecular weights;
- pooled-abundance output remains separate from the canonical sample-balanced result.

### Stage 03–04

- correct File-A class;
- opposite-strand File-B selection;
- passenger-length filters;
- official stepRNA sign convention;
- parsing of official stepRNA outputs;
- separate passenger recovery and geometry denominators.

### Stage 07 empirical guide-sequence landscape

Population/background tests:

- exact reuse of Stage 02 primary eligibility;
- exact/assigned 23/24-nt population only;
- observed antisense sequence used directly in physical 5′→3′ orientation;
- expected antisense generated only by reverse-complementing matched viral windows;
- complete A/C/G/T positional accounting for all 23/24 positions;
- expected windows exclude N/ambiguous bases and never cross FASTA-record boundaries;
- unique-sequence and abundance weightings reproduce Stage 02 definitions.

Metric tests:

- exact positional observed fraction;
- exact matched expected fraction;
- representation enrichment and delta;
- abundance-versus-unique accumulation ratio and delta;
- no pseudocount behavior at zero/undefined denominators;
- exact continuous GC9–14 calculation;
- exact A10 extraction;
- position-from-3′ mapping (`L-p+1`);
- exact 6-nt regional GC fraction;
- exactly 18 regional windows for 23 nt and 19 for 24 nt;
- exact 5′↔3′ regional coordinate conversion;
- exact equality of direct sequence-level regional GC and the mean of six constituent positional G+C frequencies;
- exact `GC9_14` regression between the literature feature and the regional engine;
- no alternative window-width generation.

Inference tests:

- pair → sample → dataset median hierarchy;
- 5000-sample clustered bootstrap with fixed seed in regression tests;
- exact sign-test accounting for positive/negative/non-zero sample deltas;
- BH validation-family correction;
- BY primary single-position discovery-family correction;
- separate single-position length/endpoint/strand correction families;
- BY primary regional-GC discovery-family correction;
- separate regional length/endpoint/strand correction families;
- `GC9_14` excluded from the regional exploratory family while retained in the full regional output.

Regression:

- position 1 = 5p1;
- position 2 = 5p2;
- position L-1 = 3p2;
- position L = 3p1;
- exact Stage 02 terminal observed/expected/enrichment regression within 1e-12.

### Stage 08 generic candidate biophysics

Accessibility tests:

- `_lunp` end-coordinate/interval-length indexing on a hand-checkable synthetic profile;
- probability range 0–1;
- one RNAplfold run per target × W/L parameter set;
- exact Stage 06 candidate accounting;
- correct extraction for arbitrary candidate lengths;
- DNA/RNA transcript normalization inherited from the locked target;
- W/L effective-value handling for transcripts shorter than requested W;
- failure when requested unpaired length exceeds effective W;
- exact min/max/range summaries.

Thermodynamic-asymmetry tests:

- exact reverse-complement orientation;
- guide 5′ end and passenger 5′ end identified correctly;
- 4-bp feature = exactly 3 Turner stack terms;
- 5-bp feature = exactly 4 Turner stack terms;
- `ΔΔG = guide - passenger`;
- positive sign means guide 5′ is less stable;
- hand-calculated canonical Watson–Crick examples;
- rejection of missing/non-finite/sentinel energies;
- non-Vd-CHIBIN synthetic candidate regression.

Vd-CHIBIN accounting regression:

- 688 23-nt rows;
- 687 24-nt rows;
- 1,375 joined rows;
- no candidate loss.

### Stage 06 generic target/candidate preparation

Generic tests:

- manifest parsing for one and multiple targets;
- single-record and multi-record FASTA selection by `fasta_record_id`;
- DNA and RNA input normalization to uppercase DNA alphabet;
- expected-length and SHA-256 validation;
- rejection of unresolved ambiguous bases;
- arbitrary legal candidate lengths;
- generic `L-w+1` enumeration;
- globally unique target-aware candidate IDs;
- exact target slicing;
- exact antisense reverse-complement orientation;
- annotation-present, partial-annotation and annotation-unavailable behavior;
- `unannotated` gap handling;
- no filtering of cross-boundary candidates;
- no hard-coded Vd-CHIBIN constants in core enumeration logic.

Vd-CHIBIN regression fixture tests:

- target ID `Vd_CHIBIN`;
- transcript `XM_022792159.1`;
- length 710;
- exact locked SHA-256;
- 688 23-nt candidates;
- 687 24-nt candidates;
- 1,375 combined candidates;
- known UTR/CDS boundary accounting for both lengths.


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
- Stage 09 feature integration and normalization;
- Stage 10 candidate ranking/selection;
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
7. canonical Stage 05 recomputes the same spatial endpoints with sample-balanced inference;
8. all major choices are configuration-controlled;
9. all critical calculations have deterministic tests;
10. all metrics are defined in `docs/METRIC_DICTIONARY.md`;
11. one Snakemake entry point regenerates the entire downstream viral analysis;
12. no manual movement or copying of intermediate outputs is required;
13. any deviation from historical results is surfaced explicitly rather than hidden;
14. Stage 06 deterministically enumerates every requested candidate length for every registered transcript target without code changes per gene; the current Vd-CHIBIN fixture reproduces 688 23-nt + 687 24-nt = 1,375 rows without ranking or upstream viral analysis;
15. Stage 07 reproduces the Stage 02 terminal landscape exactly while extending matched-background sequence association analysis to every 23/24-nt position and every fixed 6-nt regional-GC window, with literature-guided A10/GC9–14 validation, conservative dependent-test multiple-testing control, explicit pilot-discovery disclosure, and no efficacy/ranking claim;
16. Stage 08 preserves every Stage 06 row while calculating only raw target-accessibility and duplex-end-asymmetry features, with no candidate filtering, scoring, ranking, or Stage 00–05 rerun.

---

## 11. Methodological references

The main methodological basis includes:

- Murcott B, Pawluk RJ, Protasio AV, Akinmusola RY, Lastik D, Hunt VL. 2022. *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. Frontiers in Bioinformatics 2:994871. DOI: 10.3389/fbinf.2022.994871.
- Benjamini Y, Hochberg Y. 1995. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. Journal of the Royal Statistical Society Series B 57:289–300.
- Phipson B, Smyth GK. 2010. *Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn*. Statistical Applications in Genetics and Molecular Biology 9:Article 39.
- Saravanan V, Berman GJ, Sober SJ. 2020. *Application of the hierarchical bootstrap to multi-level data in neuroscience*. Used here as general methodological support for respecting nested/clustered observations; the biological application in this project is different.
- Damayo J et al. 2026 preprint, *Primary and secondary antiviral RNAi responses throughout Varroa destructor life stages reveal the vertical transmission of viruses*. This is biological motivation for the historical 23→24 hypothesis. Mechanistic claims must remain limited to what the available data directly support.
- Cedden D, Güney G, Rostás M, Bucher G. 2025. *Optimizing dsRNA sequences for RNAi in pest control and research with the dsRIP web platform*. BMC Biology 23:114. DOI: 10.1186/s12915-025-02219-6. Used here as the external source for the A10 and continuous GC9–14 hypotheses; its insect efficacy results are not treated as Varroa validation.
- Goh E, Okamura K. 2019. *Hidden sequence specificity in loading of single-stranded RNAs onto Drosophila Argonautes*. Nucleic Acids Research 47:3101–3116. DOI: 10.1093/nar/gky1300. Used as methodological precedent that positional/sequence composition can influence Argonaute-associated RNA populations; the present total-small-RNA data do not provide an equivalent AGO-loading assay.
- Jayaprakash AD, Jabado O, Brown BD, Sachidanandam R. 2011. *Identification and remediation of biases in the activity of RNA ligases in small-RNA deep sequencing*. Nucleic Acids Research 39:e141. DOI: 10.1093/nar/gkr693. Supports the explicit library-preparation-bias limitation.
- Benjamini Y, Yekutieli D. 2001. *The control of the false discovery rate in multiple testing under dependency*. Annals of Statistics 29:1165–1188. DOI: 10.1214/aos/1013699998. Supports conservative FDR control for the strongly dependent overlapping-window discovery families.

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
