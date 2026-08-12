# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.2
**Status:** Scientific specification prior to canonical implementation
**Scope:** Validated viral small-RNA analysis through viral spatial/transitivity-consistency analysis
**Host transitivity:** Excluded
**vdCHIBIN ranking:** Excluded from this build

---

# 1. Purpose

This repository contains the canonical, reproducible downstream analysis of *Varroa destructor* viral small-RNA sequencing data.

The expensive preprocessing, virus discovery, consensus generation and strict mapping stages have already been completed and independently audited.

The existing legacy project is therefore treated as a **read-only validated data core**.

The canonical repository will regenerate the current biological analyses without unnecessarily repeating the ~60 GB upstream processing.

Conceptual workflow:

```text
Validated legacy core
        ↓
00 Legacy-core validation
        ↓
01 Fixed 23/24-nt populations
        ↓
02 Terminal nucleotide enrichment
        ↓
03 Official stepRNA
        ↓
04 Dicer-feature analysis
        ↓
05 Viral spatial/transitivity-consistency analysis
```

---

# 2. Frozen legacy data core

The external dataset location is supplied locally through:

```text
config/paths.local.yaml
```

Example:

```text
/Users/patrickmod/Desktop/varroa_all_samples_pipeline_v1.0.0
```

This machine-specific path must not be committed to GitHub.

The legacy directory is treated as **read-only**.

## Validated reusable layers

The completed audit supports reuse of:

* 21 corrected processed small-RNA FASTQs
* 21 corrected preprocessing audit records
* 273 virus-discovery BAMs
* approved viral reference metadata
* selected sample-virus manifest
* sample-specific final viral consensuses
* sample-specific depth-masked background consensuses
* 21 competitive exact SAM files
* 21 competitive one-mismatch SAM files
* 21 read-level feature tables
* eligibility and mapping-summary tables

The audit reported:

```text
21/21 raw FASTQs valid
21/21 corrected processed FASTQs valid
273/273 discovery mappings present
42/42 strict SAMs structurally valid
21/21 read-level tables valid
0 warnings
0 failures
```

Therefore these expensive layers are frozen unless a later validation identifies a specific problem.

---

# 3. General statistical principles

These rules apply throughout the canonical pipeline.

## Biological clustering

Millions of reads from one library are not millions of biological replicates.

Likewise, several virus observations originating from the same sequencing library share:

* the same mite/sample context
* library preparation
* sequencing depth
* technical biases

Therefore sample-virus observations from the same sample must not automatically be treated as completely independent replicates.

## Primary uncertainty framework

Where cross-library uncertainty is required, the canonical pipeline will use **sample-clustered resampling**.

Conceptually:

```text
sample
 ├── virus A
 ├── virus B
 └── virus C
```

The sample is treated as the primary cluster.

For statistics calculated at sample-virus level:

1. calculate the required statistic within each sample-virus unit;
2. preserve those within-sample relationships;
3. resample samples with replacement;
4. recompute the overall summary.

If additional within-sample resampling is scientifically justified, it must be explicitly documented.

## Sensitivity analyses

Where useful, additionally report:

* sample-virus pair-balanced summaries
* biological-virus-balanced summaries
* leave-one-virus-out results
* abundance-weighted analyses
* true unique-sequence analyses

These are complementary views, not replacements for the sample-aware primary analysis.

## Existing mixed models

The corrected legacy pipeline already contains mixed-effects models incorporating sample, biological virus and sample-virus structure.

These may be retained as **secondary inferential/sensitivity analyses** where appropriate.

They do not replace the simpler descriptive and resampling summaries required for the fixed 23/24 analyses.

## Multiple testing

Whenever a family of related hypotheses is tested, the family must be defined before analysis.

Benjamini-Hochberg FDR correction is used for the predefined related tests unless a different correction is explicitly justified.

Both raw and adjusted P-values must be retained.

## Randomness

All bootstrap and permutation procedures must use:

* explicit random seed
* recorded number of iterations
* reproducible software implementation

---

# 00 — Validate legacy core

## Purpose

Confirm that the external frozen core still contains all inputs required by the canonical downstream workflow.

This stage performs no remapping.

## Required inputs

At minimum:

```text
results/descriptive/eligibility.tsv

tables/<sample>/<sample>.read_level_features.tsv.gz

references/consensus/
    <sample>.<virus>.final.fa
    <sample>.<virus>.final.background_masked.fa

alignments/
    <sample>.all_viruses.exact.sam
```

plus relevant:

* sample manifests
* viral metadata
* reference locks
* workflow completion records

## Required checks

Confirm:

* expected 21 libraries
* expected corrected preprocessing provenance
* expected read-level schemas
* exact SAM presence and structural validity
* required sample-virus consensuses
* depth-masked backgrounds
* consistent sample names
* consistent virus names
* no accidental write path inside `legacy_core`

## Outputs

```text
results/00_validation/
    legacy_core_validation.tsv
    legacy_core_validation.md
```

## Failure behaviour

Any missing or inconsistent dependency required by a downstream stage causes an explicit failure.

There must be no silent fallback to an older output.

---

# 01 — Fixed 23-nt and 24-nt viral small-RNA populations

## Biological question

What are the abundance, strand distribution and reproducibility of the 23-nt and 24-nt viral small-RNA populations?

## Inputs

Validated:

```text
read_level_features.tsv.gz
eligibility.tsv
viral metadata
```

No reads are remapped.

## Primary inclusion criteria

Use:

```text
mapping_mode = exact
virus_assignment = assigned
strand = sense or antisense
primary_eligible = true
length = exactly 23 or 24 nt
```

Cross-virus ambiguous reads are excluded from virus-specific primary summaries.

## Required dimensions

Retain:

* sample
* virus
* biological virus
* viral polarity
* strand
* length
* sequence
* abundance

## Weighting modes

### Abundance mode

Repeated observations of the same sequence retain their observed abundance.

Question answered:

> What dominates the accumulated sequenced small-RNA population?

### Unique-sequence mode

Within the specified biological analysis unit, each distinct sequence contributes once.

Question answered:

> Is the pattern distributed across many distinct sequence species?

These modes must remain completely separate.

## Required summaries

For 23 nt and 24 nt separately:

* total abundance
* distinct sequence count
* sense abundance
* antisense abundance
* sense fraction
* antisense fraction
* pair-level 23:24 relationship
* cross-sample summary

## Primary uncertainty

Use sample-clustered resampling for across-dataset confidence intervals.

Pair-balanced summaries remain useful descriptive/sensitivity outputs.

## Outputs

```text
results/01_viral_23_24/
```

including:

```text
23_24_counts_by_pair.tsv
23_24_fractions_by_pair.tsv
23_24_strand_bias_by_pair.tsv
23_24_across_samples.tsv
23_24_across_pairs.tsv
```

and clear descriptive figures.

---

# 02 — Length-matched terminal nucleotide enrichment

## Biological question

Do observed Varroa 23- and 24-nt viral small RNAs contain terminal nucleotide patterns more or less frequently than expected from the viral sequence actually available for processing?

## Inputs

* exact eligible 23/24 read-level data
* sample-specific depth-masked viral consensuses
* eligibility table

The corrected 23/24 implementation already established this matched-background principle.

## Positions

For each physical sequenced RNA:

```text
5p1 = first nucleotide
5p2 = second nucleotide
3p2 = penultimate nucleotide
3p1 = final nucleotide
```

Terminal bases must refer to the physical sequenced RNA orientation.

## Observed frequency

For each:

```text
sample × virus × length × strand × weighting
```

calculate observed frequencies of:

```text
A
C
G
U
```

at each terminal position.

## Expected frequency

Enumerate every fully supported window of the **same length** in the sample-specific depth-masked viral consensus.

### Sense

Use reference orientation.

### Antisense

Use the reverse-complement orientation.

### Combined-strand result

Do **not** average sense and antisense expectations 50:50.

Use the observed sense/antisense mixture for that biological unit to weight the two expected backgrounds.

Conceptually:

```text
expected_combined
=
w_sense × expected_sense
+
w_antisense × expected_antisense
```

where the weights correspond to the observed strand mixture for the relevant analysis unit.

## Enrichment ratio

```text
enrichment_ratio
=
observed_fraction
/
expected_fraction
```

Interpretation:

```text
1.0   = observed as often as sequence availability predicts
>1.0  = enriched
<1.0  = depleted
```

## Primary across-dataset summary

For each:

```text
length
strand scope
terminal position
nucleotide
weighting mode
```

report:

* median enrichment
* sample-aware bootstrap CI
* number of contributing samples
* number of contributing sample-virus units

The design-facing statistic remains:

```text
median_enrichment_ratio
```

## 23-vs-24 comparison

Calculate correlations between the matched 23- and 24-nt enrichment landscapes.

At minimum:

* overall
* antisense-specific

Use Spearman correlation as the main rank-association summary.

## Outputs

```text
results/02_terminal_enrichment/
```

including:

```text
fixed_length_positional_nucleotides_by_pair.tsv
fixed_length_positional_nucleotides_across_samples.tsv
fixed_length_positional_nucleotides_across_pairs.tsv

ALL_VIRUSES_23nt_positional_nucleotide_ratios.tsv
ALL_VIRUSES_24nt_positional_nucleotide_ratios.tsv

ALL_VIRUSES_23nt_strand_specific_positional_nucleotide_ratios.tsv
ALL_VIRUSES_24nt_strand_specific_positional_nucleotide_ratios.tsv

23_vs_24_enrichment_correlations.tsv
```

and associated figures.

## Interpretation limit

This measures empirical enrichment among sequenced Varroa vsiRNAs.

It does not independently distinguish effects caused by:

* Dicer cleavage
* Ago loading
* strand selection
* RNA stability
* library bias
* another biological process

---

# 03 — Official stepRNA Dicer-overhang analysis

## Biological question

Do the 23- and/or 24-nt populations show guide/passenger duplex-end geometry consistent with Dicer-like processing?

## Primary software

Use the **official published stepRNA implementation** as the primary overhang detector.

The previous pipeline-native reconstruction may be retained only as:

* an independent validation
* a diagnostic comparison
* a source for project-specific aggregation

It is not the primary overhang-calling method.

## Why analyse the full overhang spectrum

The analysis must not define Dicer processing solely as:

```text
2-nt 3′ overhang = Dicer
anything else = non-Dicer
```

Instead, official stepRNA is used to examine the complete:

* 5′ overhang/underhang distribution
* 3′ overhang/underhang distribution
* passenger-length distribution

The previously observed ~2-nt Dicer-like geometry is a **pre-specified biologically important signature**, but it is evaluated within the complete distribution.

## File A — reference population

Run separately for:

```text
23-nt sense
23-nt antisense
24-nt sense
24-nt antisense
```

within each eligible sample-virus unit.

## File B — potential passengers

Canonical Varroa setting:

```text
15–30 nt
```

from:

* the same sample
* the same virus
* the opposite mapped strand

The opposite-strand restriction is a Varroa project-specific biological filter.

The 15–30-nt range is a project parameter supported by the previous analysis and is also within ranges used in published stepRNA demonstrations.

## Passenger-range sensitivity analysis

Repeat the key population-level conclusions with:

```text
18–28 nt
```

as a sensitivity analysis.

This tests whether the inferred Dicer signal depends strongly on the outer passenger-length boundaries.

The canonical analysis remains 15–30 nt unless results reveal a clear methodological problem.

## Alignment behaviour

Use official stepRNA's conservative exact reverse-complementary alignment behaviour.

Do not independently introduce mismatches into the canonical run.

Any mismatch-tolerant analysis is exploratory and must be separately labelled.

## Collapsed and non-collapsed analysis

Run both:

### Collapsed / unique-reference analysis

Each distinct reference sequence is represented once.

Question:

> Is Dicer-like geometry distributed across diverse sequence species?

### Non-collapsed / expression-weighted analysis

Repeated reference sequences retain abundance.

Question:

> How much of the accumulated sequenced population carries the geometry?

## stepRNA outputs retained

Retain official outputs including:

* complete 5′ overhang/underhang counts
* complete 3′ overhang/underhang counts
* unique-overhang counts
* passenger number
* passenger length
* overhang type
* relevant classified BAMs
* official enrichment/log-odds statistics
* official Wald Z-scores

stepRNA calculates enrichment relative to the mean end-distance count and uses Wald Z-scores for its overhang-enrichment inference.

## Canonical pre-specified Varroa signature

The existing Varroa analysis found a joint geometry represented relative to the reference strand as approximately:

```text
5′: +2 nt underhang
3′: -2 nt overhang
```

corresponding to the project label:

```text
5p_underhang_2__3p_overhang_2
```

This remains a pre-specified secondary statistic for comparison with previous results.

It does not replace the full official stepRNA distribution.

## Outputs

```text
results/03_steprna/
    inputs/
    raw_outputs/
    summaries/
    sensitivity/
```

Summary tables should include:

```text
steprna_overhang_by_pair.tsv
steprna_unique_overhang_by_pair.tsv
steprna_passenger_length_by_pair.tsv
steprna_canonical_signature_by_pair.tsv
steprna_23_vs_24_summary.tsv
```

## Interpretation limit

A reproducible enriched duplex-end geometry supports Dicer/Dicer-like processing.

The analysis does not:

* directly observe intact biological duplexes
* identify the responsible Dicer paralogue
* prove every RNA of that length was produced by Dicer

Failure to recover a passenger is weaker evidence than successful recovery because passenger strands can disappear during small-RNA maturation.

---

# 04 — Dicer evidence aggregation and Dicer-conditioned sequence features

This stage has two separate purposes.

---

## 04A — Population-level Dicer evidence

### Question

How reproducible is the stepRNA Dicer-like signal across the Varroa dataset, and how does the signal compare between 23 and 24 nt?

### Primary input

Official stepRNA results from Stage 03.

### Primary summaries

For each sample-virus unit and for each:

```text
23 sense
23 antisense
24 sense
24 antisense
```

retain:

* dominant 5′ distance
* dominant 3′ distance
* stepRNA Z-score
* canonical 2-nt support
* fraction with predicted passenger
* passenger-length profile

### Sample-aware aggregation

Calculate sample-level summaries first where multiple viruses occur within one sample.

Use sample-clustered resampling for primary confidence intervals.

### 23-versus-24 comparison

Calculate paired differences wherever both populations are available within the same biological unit.

Then aggregate these differences with sample clustering preserved.

### Existing custom statistic

The previous analysis used:

```text
Δ_Dicer
=
support at the pre-specified Dicer distance
−
mean support at alternative tested distances
```

This may be retained as a **secondary project-specific validation statistic**.

It must always be labelled:

```text
custom Varroa summary statistic
```

rather than an official stepRNA statistic.

### Sensitivity outputs

Also report:

* pair-balanced summaries
* virus-balanced summaries
* leave-one-virus-out where informative
* 15–30 versus 18–28 passenger-range comparison

---

## 04B — Dicer-conditioned sequence-feature analysis

### Question

Do small RNAs with strong Dicer-like support have sequence features that are different from the general 23/24 small-RNA population?

This is the analysis that could eventually justify an additional candidate-ranking metric.

### Pre-specified Dicer-supported subset

For the first canonical analysis, use a pre-specified Dicer-compatible geometry based on the prior Varroa result rather than selecting whichever distance happens to maximize enrichment in the same dataset.

Primary strong subset:

```text
canonical joint 2-nt geometry
```

Additional subsets based on other reproducibly enriched stepRNA peaks may be explored but must be labelled exploratory until independently validated.

### Features examined initially

Pre-specify a limited feature set:

```text
5p1
5p2
3p2
3p1
guide length
reference strand
passenger length
```

Additional sequence features should be added only with clear biological justification.

### Absolute Dicer-conditioned enrichment

For a sequence feature:

```text
E_Dicer_absolute
=
frequency among Dicer-supported RNAs
/
viral-sequence expected frequency
```

This asks:

> What features characterize the Dicer-supported subset relative to sequence availability?

### Dicer-specific contrast

Also calculate:

```text
Dicer_specific_contrast
=
E_Dicer_absolute
/
E_all_observed
```

where:

```text
E_all_observed
```

is the general enrichment from Stage 02 for the same:

* length
* strand
* position
* nucleotide
* biological unit

This asks the more useful question:

> Is this feature specifically associated with the Dicer-supported subset beyond the nucleotide preference already present in all observed siRNAs?

Interpretation:

```text
≈1
Dicer-supported RNAs do not add much information beyond general enrichment

>1
feature is more strongly represented among Dicer-supported RNAs

<1
feature is less represented among Dicer-supported RNAs
```

### Stratification

Do not pool 23 and 24 nt.

Do not pool sense and antisense before the strand-specific analysis.

Passenger recovery differs between strand-biased populations, particularly the strongly antisense-biased 24-nt class.

### Redundancy test

Quantify correlation between:

```text
general enrichment metric
```

and:

```text
Dicer-specific metric
```

If the new Dicer metric is highly redundant with the existing enrichment metric, it must not later be treated as an independent scoring feature without justification.

### Outputs

```text
results/04_dicer_features/
```

including:

```text
dicer_population_summary_by_pair.tsv
dicer_population_summary_by_sample.tsv
dicer_23_vs_24.tsv
dicer_passenger_length_summary.tsv
dicer_conditioned_terminal_features.tsv
dicer_specific_contrasts.tsv
dicer_vs_general_enrichment.tsv
dicer_sensitivity_summary.tsv
```

---

# 05 — Viral spatial/transitivity-consistency analysis

## Biological question

Is the secondary-associated 24-nt antisense population spatially related to primary-like 23-nt processing in a pattern consistent with secondary/transitive RNAi?

Because these are natural viral infections, this stage is deliberately called:

> **viral spatial/transitivity-consistency analysis**

rather than proof of transitivity.

Viral replication itself can produce complementary RNA substrates and therefore complicates mechanistic interpretation.

---

# 05A — Coordinate reconstruction and QC

## Inputs

Use existing:

```text
exact competitive SAMs
read-level feature tables
eligibility table
viral polarity metadata
final viral references
```

No remapping.

## Required reconstruction

For every selected exact 23/24-nt read:

recover:

* sample
* assigned virus
* strand
* genomic coordinate
* sequence
* length
* mapping location information

## QC

Report:

* total selected reads
* successfully coordinate-recovered reads
* distinct recovered sequences
* strand mismatches
* sequence mismatches
* metadata conflicts
* missing alignments

The previous strengthened run recovered all selected reads with no strand or metadata conflicts; the canonical pipeline should independently verify this rather than hard-code the result.

---

# 05B — Eligible viral units

The directional analysis is restricted to suitable positive-sense viral genomes for which upstream/downstream interpretation is biologically defined in a consistent reference orientation.

The canonical pipeline must write the exact eligible set to:

```text
eligible_positive_sense_units.tsv
```

The previous strengthened analysis contained:

```text
14 samples
19 positive-sense sample-virus units
```

but the new implementation must derive this from the validated metadata.

---

# 05C — Spatial tracks

Construct position-wise tracks for:

```text
23 sense
23 antisense
24 sense
24 antisense
```

for each eligible sample-virus unit.

Maintain separate:

* abundance tracks
* true unique-sequence tracks

## True unique-sequence requirement

A unique-sequence analysis must genuinely deduplicate sequence identities at the defined analysis level.

It must not reproduce the earlier error where a nominal unique-sequence mode remained effectively read-weighted.

---

# 05D — Primary-like 23-nt anchor tracks

Maintain the two established anchor concepts.

## `balanced23`

```text
balanced23
=
sqrt(23_sense × 23_antisense)
```

This preferentially identifies positions/regions with evidence from both strands.

## `combined23`

```text
combined23
=
23_sense + 23_antisense
```

This is a broader measure of local 23-nt activity.

They remain separate analyses.

---

# 05E — Implementation lock before Stage 05 coding

The archived project documentation confirms the final v1.4.1 endpoints and statistical safeguards, but the currently recovered documentation does **not specify with sufficient precision**:

1. the exact hotspot/anchor-selection algorithm;
2. the exact threshold used to call a hotspot;
3. whether neighbouring anchor positions were merged and how;
4. the exact treatment of within-virus positional multimappers;
5. the exact boundary behaviour of the positional-shift null;
6. the exact family of tests used for the previous BH correction.

These details must **not be guessed**.

Before implementing Stage 05, recover the final v1.4.1 script/result package and lock these definitions into this specification.

Until this is done:

```text
Stage 05 scientific logic = specified
Stage 05 exact implementation = NOT YET LOCKED
```

This is an intentional reproducibility safeguard.

## If the final v1.4.1 implementation cannot be recovered

A new prospective method must be defined explicitly before looking at new results.

A conservative option could then use:

* unique-position mappings as the primary spatial analysis
* multimappers as a sensitivity analysis
* a pre-specified hotspot threshold
* a pre-specified merging rule
* a pre-specified permutation boundary rule

However, this would constitute a **new method**, not a reconstruction of v1.4.1, and must be labelled accordingly.

---

# 05F — Canonical distances

Retain the previously tested distances:

```text
100 nt
250 nt
500 nt
```

These are the confirmatory canonical windows.

Any additional distance is exploratory unless separately pre-specified before examining results.

---

# 05G — Primary endpoint: antisense-specific directionality

Define downstream-minus-upstream signal for each 24-nt strand.

Then use:

```text
D_24AS − D_24S
```

where:

```text
D_24AS =
24AS_downstream − 24AS_upstream

D_24S =
24S_downstream − 24S_upstream
```

Purpose:

> test whether downstream behaviour is specifically stronger for the antisense 24-nt population rather than reflecting a generic directional property of the viral genome.

---

# 05H — Composition endpoint

Define:

```text
F24_AS
=
24AS
/
(23AS + 24AS)
```

Then:

```text
ΔF24_AS
=
F24_AS_downstream
−
F24_AS_upstream
```

This asks:

> Does the downstream antisense population become relatively more 24-nt dominated?

## Zero-denominator rule

If:

```text
23AS + 24AS = 0
```

then:

```text
F24_AS = NA
```

There is no biological composition to estimate.

Do not generate a 0.5 value using a pseudocount.

---

# 05I — Aggregation

Primary inference must preserve sample clustering.

Report:

### Primary

* sample-aware/clustered summary

### Sensitivity

* sample-virus pair-balanced
* biological-virus-balanced
* abundance mode
* unique-sequence mode
* balanced23 anchors
* combined23 anchors

---

# 05J — Leave-one-virus-out analysis

Repeat the relevant higher-level summary while omitting each biological virus in turn.

Purpose:

> determine whether one virus is responsible for the overall result.

This is a robustness analysis, not a replacement for the primary estimate.

---

# 05K — Permutation null

Preserve the established principle that:

```text
24 sense
and
24 antisense
```

must move together under the spatial null.

This preserves their internal strand relationship while disrupting their positional relationship with the 23-nt anchor structure.

The exact positional-shift implementation and boundary behaviour must be recovered from v1.4.1 before coding.

Record:

* random seed
* number of permutations
* valid permutations
* observed statistic
* empirical P-value

---

# 05L — Multiple testing

Use BH-FDR correction over the pre-specified confirmatory test family.

The exact family used in the reconstructed v1.4.1 analysis must be recovered before implementation.

Do not choose the FDR family after inspecting which tests are significant.

Retain:

```text
raw P
BH-adjusted P
effect size
uncertainty interval
```

---

# 05M — Required outputs

```text
results/05_viral_transitivity/
```

including:

```text
coordinate_qc.tsv
eligible_positive_sense_units.tsv

spatial_tracks/
anchor_summary.tsv

transitivity_by_pair.tsv
transitivity_by_sample.tsv

sample_clustered_results.tsv
pair_balanced_results.tsv
virus_balanced_results.tsv

leave_one_virus_out.tsv
permutation_results.tsv
multiple_testing_summary.tsv

final_transitivity_summary.tsv
```

plus figures.

---

# 05N — Interpretation limits

The viral analysis may support statements such as:

* 23/24 spatial co-localization
* downstream/upstream asymmetry
* an antisense-specific directional effect
* a shift in antisense composition toward 24 nt
* consistency with secondary/amplification-associated biology

It must not by itself establish:

* that every 23-mer is primary
* that every 24-mer is secondary
* that RdRP directly synthesizes 24-mers
* that 24-mers are Dicer-independent
* a universal secondary propagation distance
* a specific Dicer/Ago/RdRP paralogue
* host-mRNA transitivity

---

# 6. Reproducibility requirements

## External data

The validated legacy core remains read-only.

## New results

All new files must be written under:

```text
results/
```

inside the canonical Git project.

## Configuration

Analysis parameters should live in configuration rather than being hidden inside scripts.

At minimum record:

```text
target_lengths = [23, 24]

steprna_passenger_range = [15, 30]
steprna_sensitivity_range = [18, 28]

transitivity_distances = [100, 250, 500]

bootstrap_iterations
permutation_iterations
random_seed
```

Once recovered, also record:

```text
hotspot threshold
anchor merge rule
multimapper policy
permutation boundary rule
FDR test family
```

## Software versions

Record exact versions of:

* Python
* Snakemake
* stepRNA
* Bowtie2 used internally by stepRNA
* pysam
* pandas
* numpy
* scipy
* statistical packages
* plotting packages

Only software actually used by this downstream build needs to be installed.

## Provenance

Each run should record:

* pipeline Git commit
* configuration
* software versions
* legacy-core path
* relevant input-file identity
* run date
* random seed

---

# 7. Required tests

Deterministic tests should cover at least:

## Terminal enrichment

* 5p1 extraction
* 5p2 extraction
* 3p2 extraction
* 3p1 extraction
* reverse-complement antisense expectation
* observed-strand-weighted combined expectation

## Unique-sequence logic

* repeated sequence counted once
* different sequences counted separately
* abundance mode unchanged

## Dicer analysis

* correct File-A strand
* opposite-strand File-B selection
* correct 23/24 separation
* passenger-length filtering
* canonical sign interpretation
* correct parsing of official stepRNA outputs

## Transitivity

* positive-sense coordinate orientation
* upstream/downstream assignment
* true sequence deduplication
* zero-denominator → NA
* paired movement of 24S/24AS in null
* deterministic seed behaviour

## General

* missing required input → explicit failure
* no silent fallback to legacy downstream results

---

# 8. Explicitly superseded analyses

The following old outputs are historical references only:

```text
legacy fixed_length_23_24 summaries
legacy custom Dicer-overhang summaries
older viral transitivity implementations
```

They may be used for numerical comparison during validation.

They must not silently become inputs to the new canonical results.

---

# 9. Explicitly outside the current scope

This build does not contain:

* CHH analysis
* Pero analysis
* Ago2 host analysis
* host transitivity
* inferred host-trigger analysis
* ViennaRNA accessibility
* vdCHIBIN thermodynamic asymmetry
* vdCHIBIN scoring
* construct architecture comparison
* Nectar Designer integration

Those are separate future modules.

---

# 10. Definition of success

The first canonical viral pipeline is complete when:

1. a fresh Git clone can point to the validated external legacy core;
2. Stage 00 passes without changing the legacy directory;
3. Stages 01–04 regenerate all 23/24, enrichment and Dicer outputs from frozen inputs;
4. official stepRNA is the primary Dicer-overhang caller;
5. the full stepRNA overhang distributions are retained rather than reducing the analysis to a single 2-nt number;
6. sample clustering is respected in primary cross-dataset uncertainty estimates;
7. the corrected strand-weighted background is used for combined enrichment;
8. Dicer-conditioned sequence features are explicitly tested for information beyond general nucleotide enrichment;
9. the exact final v1.4.1 anchor/multimapper/permutation definitions are recovered and locked before Stage 05 is implemented;
10. Stage 05 reproduces the corrected analysis logic rather than an older transitivity version;
11. all major calculations have deterministic tests;
12. all metrics are defined in `METRIC_DICTIONARY.md`;
13. one Snakemake entry point regenerates the complete downstream viral analysis;
14. no manual copying of intermediate files is required;
15. any difference from historical results is reported rather than hidden.

---

# 11. Academic method basis

Primary methodological references include:

* Murcott B. et al. (2022). *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. **Frontiers in Bioinformatics 2:994871.** DOI: 10.3389/fbinf.2022.994871.
* Saravanan V., Berman G.J., Sober S.J. (2020). *Application of the hierarchical bootstrap to multi-level data in neuroscience*. Used here as general methodological support for respecting clustered/nested biological observations during resampling.
* Benjamini Y. & Hochberg Y. (1995). *Controlling the false discovery rate: a practical and powerful approach to multiple testing.*
* The corrected Varroa v1.2.1 core and v1.3 matched 23/24 analysis remain the project-specific source for eligibility, strict mappings, depth-masked backgrounds and strand-weighted nucleotide expectations.

---

# 12. Canonical workflow summary

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
```

This is the complete scope of the first canonical downstream Varroa viral small-RNA pipeline.
