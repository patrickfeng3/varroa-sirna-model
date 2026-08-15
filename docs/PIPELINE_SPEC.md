# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.6  
**Status:** Stages 00–02 implemented and validated; Stage 03 scientifically specified before implementation  
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

## 03.10 Full signed geometry spectrum is the primary Stage 03 geometry output

Stage 03 must retain the complete signed 5′ and 3′ distance spectra produced by official stepRNA rather than reducing the analysis to a single 2-nt category.

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

where the exact official field names from stepRNA 1.0.6 are preserved in raw outputs and mapped transparently into the parsed schema.

The full spectrum is primary because Dicer-associated end geometry is pathway- and organism-dependent; Stage 03 must first describe the Varroa viral data rather than define Dicer as “2 nt or nothing.”

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
- full 5′ distance spectrum;
- full 3′ distance spectrum;
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
- consistency of joint-geometry counts with the classified alignments used to derive them.

A zero-passenger biological run is valid data and is not automatically a pipeline failure.

A software failure, missing required output, parser inconsistency, or identifier mismatch is a failure.

---

## 03.19 Outputs

```text
results/03_steprna/
│
├── qc/
│   └── stage03_accounting.tsv
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
│   ├── joint_geometry_by_pair.tsv
│   └── joint_geometry_references.tsv.gz
│
└── sensitivity/
    └── passenger_18_28/        # absent unless explicitly run
```

Stage 03 does not need publication figures before the calculations and parser have been validated.

---

## 03.20 Interpretation limits

A reproducibly enriched duplex-end geometry is evidence **consistent with Dicer/Dicer-like processing**.

Stage 03 does not directly observe cleavage in vivo and does not by itself:

- identify a specific Varroa Dicer/Dicer-like protein;
- prove that every RNA in a length class was Dicer-generated;
- prove that an RNA lacking a recovered passenger was not Dicer-generated;
- establish that 23 nt is “primary” and 24 nt is “secondary”;
- establish RdRP-dependent amplification;
- define a candidate-window scoring metric.

Those higher-level biological comparisons and Dicer-conditioned sequence analyses begin in Stage 04.

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
