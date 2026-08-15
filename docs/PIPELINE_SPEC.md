# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.12  
**Status:** Stages 00–05 implemented and validated; Stage 06 scientifically specified before implementation  
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
01 Viral length landscape and fixed 23/24-nt populations
        ↓
02 Terminal nucleotide enrichment
        ↓
03 Official stepRNA
        ↓
04 Duplex-geometry evidence and geometry-conditioned features
        ↓
05 Viral spatial/transitivity-consistency analysis
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

# 06 — Vd-CHIBIN target preparation and exhaustive 23/24-nt candidate enumeration

## 06.1 Purpose

Stage 06 establishes the **canonical vdCHIBIN candidate universe** before any accessibility, thermodynamic, empirical-enrichment, ranking, or construct-level scoring is applied.

It answers only:

> What exact 23-nt and 24-nt target intervals exist in the locked Vd-CHIBIN transcript, what sequences do they contain, what is the corresponding antisense guide orientation, and which transcript annotations do they overlap?

Stage 06 is deliberately sequence/coordinate preparation only.

It does **not**:

- run ViennaRNA;
- calculate strand asymmetry;
- join Stage 02 enrichment values;
- use Stage 03/04 duplex geometry;
- use Stage 05 transitivity;
- apply an asymmetry gate;
- calculate A/T/N/S;
- rank candidates;
- select positive or negative controls;
- compare 23-nt and 24-nt candidates by a single aggregate score;
- compare 24/48/96 constructs.

---

## 06.2 Locked target reference

Canonical target:

```text
project target: Vd-CHIBIN / vbchitin
RefSeq accession-version: XM_022792159.1
transcript length: 710 nt
```

Canonical transcript annotation:

```text
5′ UTR: 1–329
CDS:    330–665
3′ UTR: 666–710
```

All coordinates are **1-based inclusive**.

The canonical tracked reference bundle is:

```text
resources/vdchibin/XM_022792159.1.fasta
resources/vdchibin/XM_022792159.1.annotation.tsv
resources/vdchibin/XM_022792159.1.reference_manifest.tsv
```

The normalized uppercase DNA sequence must be exactly 710 nt, contain only `A/C/G/T`, and have:

```text
SHA-256 = 4a0d25aa05b269a118ed1b952dca63ccd1c0a7978fc42295faf3bf650e43ea42
```

where the hash is calculated over the 710 uppercase sequence characters only, with no FASTA header or newline characters.

The small tracked reference bundle is the canonical computational input. A live database fetch may be used only as an explicit reference-validation step; Stage 06 must not silently replace the locked accession-version sequence during routine analysis.

Historical custom Panel-B annotations are not part of the canonical Stage 06 annotation set.

---

## 06.3 Parallel candidate lengths

Canonical Stage 06 candidate lengths are:

```text
23 nt
24 nt
```

Both lengths are enumerated from the **same locked transcript**, using the same coordinate system, annotation rules, and guide-orientation rules.

This parallel treatment is important because Stages 01–05 established that both 23-nt and 24-nt viral small-RNA populations are biologically prominent while **not** establishing that one length uniquely represents primary/Dicer processing and the other uniquely represents secondary/RdRP processing.

Therefore:

```text
23 nt ≠ automatically “primary”
24 nt ≠ automatically “secondary”
```

The two candidate universes remain separate analysis strata through later stages unless a future specification explicitly defines a justified cross-length comparison.

Stage 06 does not prefer either length.

---

## 06.4 Exhaustive window enumeration

For transcript length:

```text
L = 710
```

and candidate length `w`, enumerate every complete interval with:

```text
start_1based = 1 ... L - w + 1
end_1based   = start_1based + w - 1
```

### 23 nt

```text
n_candidates_23 = 710 - 23 + 1 = 688
```

### 24 nt

```text
n_candidates_24 = 710 - 24 + 1 = 687
```

Total Stage 06 candidate rows:

```text
688 + 687 = 1,375
```

No candidate may extend beyond either transcript end.

No spacing, accessibility, sequence-composition, empirical-enrichment, geometry, or transitivity filter is applied at this stage.

Enumeration must be deterministic and exhaustive within each length.

---

## 06.5 Canonical candidate identifier

Each candidate receives a deterministic identifier:

```text
XM_022792159.1__LENGTHnt__START_END
```

with zero-padded 1-based inclusive coordinates.

Examples:

```text
XM_022792159.1__23nt__0001_0023
XM_022792159.1__23nt__0688_0710

XM_022792159.1__24nt__0001_0024
XM_022792159.1__24nt__0687_0710
```

Candidate identity is defined by:

```text
accession-version
candidate length
target interval coordinates
```

Sequence identity alone is not sufficient because repeated sequences at different transcript positions remain distinct target intervals.

---

## 06.6 Sequence orientation

For each candidate interval:

```text
target_sequence_dna
    = exact transcript/sense sequence, 5′→3′, DNA alphabet

target_sequence_rna
    = target_sequence_dna with T→U, 5′→3′

antisense_guide_sequence_rna
    = reverse complement of target_sequence_rna, 5′→3′
```

The desired guide orientation for later RNAi design is:

```text
antisense to the Vd-CHIBIN mRNA
```

The guide 5′ end therefore corresponds to the **target interval's 3′ end**, not its 5′ start coordinate.

Stage 06 must test this explicitly for both 23-nt and 24-nt candidates.

Do not reverse-complement the target before Stage 07 target-accessibility analysis; accessibility is evaluated on the mRNA/sense target interval.

---

## 06.7 Transcript-region annotation

Stage 06 records both **start-region grouping** and **actual interval overlap** so that the two concepts cannot be confused later.

For each candidate define:

```text
start_region
end_region
overlap_regions
crosses_annotation_boundary
```

`start_region` is the annotation containing `start_1based`.

`end_region` is the annotation containing `end_1based`.

`overlap_regions` is the ordered semicolon-separated list of all transcript annotations intersected by the full candidate interval.

Cross-boundary candidates remain valid Stage 06 candidates. They are **annotated, not discarded**.

Examples:

```text
23 nt candidate 308–330
start_region = 5_prime_UTR
end_region   = CDS
overlap_regions = 5_prime_UTR;CDS
crosses_annotation_boundary = TRUE
```

```text
24 nt candidate 307–330
start_region = 5_prime_UTR
end_region   = CDS
overlap_regions = 5_prime_UTR;CDS
crosses_annotation_boundary = TRUE
```

```text
23 nt candidate 666–688
start_region = 3_prime_UTR
end_region   = 3_prime_UTR
overlap_regions = 3_prime_UTR
crosses_annotation_boundary = FALSE
```

```text
24 nt candidate 666–689
start_region = 3_prime_UTR
end_region   = 3_prime_UTR
overlap_regions = 3_prime_UTR
crosses_annotation_boundary = FALSE
```

Any later region-specific selection rule must state explicitly whether it uses start-region grouping, complete containment, or overlap.

---

## 06.8 Expected deterministic counts

### 23-nt candidates

Total:

```text
688
```

By start coordinate region:

```text
5_prime_UTR starts = 329
CDS starts         = 336
3_prime_UTR starts = 23
total              = 688
```

Fully contained within one annotation:

```text
5_prime_UTR only = 307
CDS only         = 314
3_prime_UTR only = 23
total fully within one region = 644
```

Cross-boundary:

```text
5_prime_UTR → CDS = 22
CDS → 3_prime_UTR = 22
total cross-boundary = 44
```

### 24-nt candidates

Total:

```text
687
```

By start coordinate region:

```text
5_prime_UTR starts = 329
CDS starts         = 336
3_prime_UTR starts = 22
total              = 687
```

Fully contained within one annotation:

```text
5_prime_UTR only = 306
CDS only         = 313
3_prime_UTR only = 22
total fully within one region = 641
```

Cross-boundary:

```text
5_prime_UTR → CDS = 23
CDS → 3_prime_UTR = 23
total cross-boundary = 46
```

### Combined accounting

```text
23-nt candidates = 688
24-nt candidates = 687
total rows        = 1,375
```

These counts are deterministic QC expectations, not biological filters.

---

## 06.9 Canonical output schema

Primary candidate table:

```text
results/06_vdchibin_candidates/vdchibin_candidates.tsv
```

Required columns:

```text
candidate_id
accession
candidate_length_nt
start_1based
end_1based
target_sequence_dna
target_sequence_rna
antisense_guide_sequence_rna
start_region
end_region
overlap_regions
crosses_annotation_boundary
```

The primary table contains **both** 23-nt and 24-nt candidates.

Optional convenience exports may be generated:

```text
vdchibin_23nt_candidates.tsv
vdchibin_24nt_candidates.tsv
```

but they must be exact filtered views of the canonical combined table and must not become separate sources of truth.

Optional implementation/provenance columns may be added only if they do not alter candidate identity or silently introduce ranking features.

Do **not** add columns called:

```text
score
rank
primary_score
secondary_score
Dicer_score
transitivity_score
```

in Stage 06.

---

## 06.10 QC and deterministic tests

Stage 06 must verify at least:

### Reference

- exact accession-version string;
- sequence length = 710;
- uppercase normalized sequence contains only A/C/G/T;
- sequence SHA-256 matches the locked manifest;
- annotation intervals are contiguous, non-overlapping, ordered, and cover exactly 1–710;
- CDS = 330–665;
- UTR boundaries match the locked annotation.

### Enumeration — 23 nt

- exactly 688 candidates;
- first interval = 1–23;
- final interval = 688–710;
- every interval length = 23;
- starts increase by exactly 1;
- no duplicate candidate IDs;
- no missing start coordinate from 1–688;
- target sequence equals exact transcript slice.

### Enumeration — 24 nt

- exactly 687 candidates;
- first interval = 1–24;
- final interval = 687–710;
- every interval length = 24;
- starts increase by exactly 1;
- no duplicate candidate IDs;
- no missing start coordinate from 1–687;
- target sequence equals exact transcript slice.

### Combined table

- exactly 1,375 rows;
- candidate lengths are exactly `{23,24}`;
- exactly 688 rows have length 23;
- exactly 687 rows have length 24;
- no duplicate candidate IDs across the combined table.

### Orientation

For both lengths:

- `target_sequence_rna = target_sequence_dna.replace(T,U)`;
- antisense guide equals exact reverse complement of target RNA;
- guide length equals `candidate_length_nt`;
- guide 5′ nucleotide is complementary to the target interval's final nucleotide;
- guide 3′ nucleotide is complementary to the target interval's first nucleotide.

### Annotation — 23 nt

- start-region counts = 329 / 336 / 23;
- fully-contained counts = 307 / 314 / 23;
- exactly 44 cross-boundary candidates;
- 308–330 is labelled `5_prime_UTR;CDS`;
- 644–666 is labelled `CDS;3_prime_UTR`;
- 666–688 is fully `3_prime_UTR`.

### Annotation — 24 nt

- start-region counts = 329 / 336 / 22;
- fully-contained counts = 306 / 313 / 22;
- exactly 46 cross-boundary candidates;
- 307–330 is labelled `5_prime_UTR;CDS`;
- 643–666 is labelled `CDS;3_prime_UTR`;
- 666–689 is fully `3_prime_UTR`.

Any mismatch in reference identity, sequence hash, enumeration, or orientation is a **FAIL**.

---

## 06.11 Outputs

```text
results/06_vdchibin_candidates/
│
├── vdchibin_reference_summary.tsv
├── vdchibin_candidates.tsv
├── vdchibin_23nt_candidates.tsv          # optional exact filtered view
├── vdchibin_24nt_candidates.tsv          # optional exact filtered view
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

## 06.12 Interpretation and carry-forward

A Stage 06 row means only:

> this exact 23-nt or 24-nt interval exists in the locked Vd-CHIBIN transcript and has this exact antisense guide sequence and annotation context.

It does **not** mean that the candidate is efficient, accessible, favourably asymmetric, enriched, Dicer-compatible, transitivity-compatible, or selected.

All 1,375 candidates proceed equally into Stage 07 unless a future explicit sequence-integrity exclusion is added prospectively.

Stage 07 must retain 23 nt and 24 nt as separate analysis strata while calculating target accessibility and thermodynamic strand asymmetry using these exact candidate identities.

Any eventual comparison between 23-nt and 24-nt designs must distinguish:

```text
within-length ranking
```

from:

```text
between-length biological/design comparison
```

and must not assume raw metric values are directly comparable across lengths unless the metric definition explicitly supports that comparison.

---


## 6. Reproducibility requirements

All new outputs are written under `results/` in the canonical repository.

Configuration must record at minimum:

```text
target_lengths = [23, 24]
vdchibin_accession = "XM_022792159.1"
vdchibin_candidate_lengths_nt = [23, 24]
vdchibin_coordinate_system = "1-based inclusive"
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

### Stage 06 target/candidate preparation

- exact accession-version and 710-nt reference identity;
- locked sequence SHA-256;
- contiguous annotation coverage of 1–710;
- exact 688-window enumeration for 23 nt;
- exact 687-window enumeration for 24 nt;
- exact combined total of 1,375 candidates;
- first/last window coordinates for both lengths;
- target-sequence slicing;
- RNA T→U conversion;
- exact antisense reverse-complement orientation for both lengths;
- guide 5′/3′ endpoint correspondence;
- deterministic length-aware candidate IDs;
- start-region counts for both lengths;
- full-containment counts for both lengths;
- 44 cross-boundary 23-nt windows;
- 46 cross-boundary 24-nt windows;
- exact 23-nt boundary examples 308–330, 644–666 and 666–688;
- exact 24-nt boundary examples 307–330, 643–666 and 666–689.

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
- Stage 07 ViennaRNA vdCHIBIN accessibility/thermodynamic scoring;
- Stage 08 feature integration;
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
14. Stage 06 deterministically regenerates the locked parallel 23/24-nt Vd-CHIBIN candidate universe (688 + 687 = 1,375 rows) without running any ranking or upstream viral analysis.

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
