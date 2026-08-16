# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.21  
**Status:** Stages 00–09 implemented and validated; Stage 10A specified as the canonical individual-window integration and ranking layer; Stage 11 reserved for future region selection and visualization  
**Scope:** Canonical viral small-RNA analysis, generic transcript candidate enumeration, empirical Varroa guide-sequence association, candidate biophysics, Stage 09 three-layer evidence synthesis, and Stage 10A individual-window integration  
**Host transitivity:** Excluded from current canonical build  
**Individual candidate-window ranking:** Stage 10A; **region/construct selection:** deferred to Stage 11 or later

---

# 1. Purpose

This repository contains the canonical, reproducible downstream analysis of *Varroa destructor* viral small-RNA sequencing data and the design-facing calculations used to evaluate antisense candidate windows for target transcripts such as Vd-CHIBIN.

The expensive upstream work—read preprocessing, virus discovery, sample-specific consensus generation, strict mapping, and read-level feature extraction—has already been completed and audited. The legacy project is therefore treated as a **read-only validated data core**.

Canonical repository:

```text
/Users/patrickmod/Documents/GitHub/varroa-sirna-model
```

Frozen legacy core:

```text
/Users/patrickmod/Desktop/varroa_all_samples_pipeline_v1.0.0
```

Core principle:

> Rebuild biological conclusions from validated frozen upstream inputs in a clean sequence. Historical analyses are regression/reference targets, not analytical inputs.

Conceptual workflow:

```text
00  validate frozen legacy core

01  viral 15–35-nt length landscape
    + focused 23/24-nt population analysis

02  terminal nucleotide enrichment
    against matched viral sequence opportunity

03  official stepRNA duplex geometry

04  sample-aware duplex-geometry aggregation
    + geometry-conditioned sequence analysis

05  viral spatial/transitivity-consistency analysis

06  generic mature-transcript preparation
    + exhaustive candidate enumeration

07  Varroa empirical guide-sequence association landscape
    + single-position features
    + fixed-width regional GC
    + Wang/Bartel feature synthesis

08  generic candidate biophysics
    + whole-site target accessibility
    + guide g2–g8 target accessibility
    + duplex-end thermodynamic asymmetry
    + guide self-folding

09A Layer 1 — Varroa small-RNA accumulation propensity
09B Layer 2 — guide competence / strand-selection biophysics
09C Layer 3 — target engagement / predicted accessibility

10A individual-window integration
    + equal-layer ranking
    + Pareto-front diagnostics

11  future user-configurable region selection
    + visualization
    + long-dsRNA / construct analysis only when separately specified
```

No historical arbitrary `0.6 / 0.3 / 0.1` weighting is inherited.

---

# 2. Workflow and repository safety

The frozen legacy core is immutable unless explicit regeneration is requested.

Mandatory workflow for new implementation:

1. edit narrowly;
2. run targeted deterministic tests only;
3. dry-run the **exact intended target**;
4. the dry run must show only the intended stage;
5. if unexpected upstream Stage 00–08 jobs appear, **STOP**;
6. do not repair/reconstruct/delete upstream outputs automatically;
7. execute only the exact target after dry-run PASS;
8. do not casually run broad `snakemake` or `snakemake all`;
9. generated outputs remain under canonical `results/`;
10. machine-specific paths are not committed.

Before expensive execution, explicitly distinguish:

```text
READS
WRITES
MUST NOT RUN
```

---

# 3. General analysis principles

## 3.1 Evidence-building principle

For each stage:

1. define the analysis from this specification;
2. calculate from validated upstream inputs;
3. interpret the newly generated result;
4. only then compare with historical results.

Historical result tables must never be used as analytical inputs to recreate the result they are supposed to validate.

## 3.2 Analysis levels

```text
sample
sample-virus unit
sample-virus-contig unit
candidate
target × candidate_length stratum
```

A sample is the top-level biological clustering unit for cross-dataset inference because multiple viruses can occur in one library.

## 3.3 Weighting modes

Two empirical viral weighting modes remain distinct:

```text
abundance
unique_sequence
```

Abundance mode uses observed read abundance.

Unique-sequence mode gives each distinct RNA sequence total weight 1 within the explicitly defined analysis unit.

## 3.4 Uncertainty

Canonical cross-dataset uncertainty must respect biological samples as top-level clusters.

Where appropriate:

```text
sample-virus metric
→ median across viruses within each sample
→ median across samples
→ sample-clustered bootstrap
```

## 3.5 Randomness

All bootstrap/permutation analyses must record:

- seed;
- requested replicate count;
- valid replicate count;
- aggregation rule;
- multiple-testing family.

## 3.6 Mechanistic restraint

The pipeline distinguishes:

```text
observed small-RNA association
predicted biophysical property
mechanistic consistency
measured efficacy
```

These are not interchangeable.

---

# 4. Frozen legacy data core

The external legacy path is supplied locally through:

```text
config/paths.local.yaml
```

Validated reusable layers include:

- corrected processed small-RNA FASTQs;
- preprocessing audit records;
- virus-discovery mappings;
- approved virus metadata and selected sample-virus manifest;
- sample-specific final viral consensuses;
- depth-masked background consensuses;
- competitive exact and one-mismatch mapping files;
- read-level feature tables;
- eligibility/mapping-summary tables.

The canonical pipeline must not silently reuse old downstream 23/24, Dicer, sequence-enrichment, or transitivity result summaries.

---

# 00 — Validate frozen legacy core

## Purpose

Confirm that all required frozen inputs and schemas exist without remapping or biological inference.

Minimum required dependencies include:

```text
results/descriptive/eligibility.tsv
config/virus_catalog.tsv
tables/<sample>/<sample>.read_level_features.tsv.gz
alignments/<sample>.all_viruses.exact.sam
references/consensus/<sample>.<analysis_unit>.final.fa
references/consensus/<sample>.<analysis_unit>.final.background_masked.fa
```

Required checks:

- expected libraries present;
- corrected provenance present;
- required read-level columns present;
- exact mapping files structurally readable;
- required consensus/background records present;
- identifiers agree across inputs;
- output paths do not resolve inside legacy core;
- frozen input identity recorded where practical.

Outputs:

```text
results/00_validation/
    legacy_core_validation.tsv
    legacy_core_validation.md
```

Any missing/inconsistent required dependency is a hard failure.

---

# 01 — Viral length landscape and focused 23/24 analysis

## 01.0 Purpose

Reconstruct the complete retained viral small-RNA length spectrum before focusing on 23 and 24 nt.

Canonical retained range:

```text
15–35 nt
```

Primary population:

```text
20 primary samples
54 sample-virus units
```

No pathway identity is assigned solely from length.

---

## 01A — 15–35-nt length landscape

Canonical inclusion:

```text
primary_eligible sample × analysis_unit
mapping_mode = exact
virus_assignment = assigned
strand ∈ {sense, antisense}
length ∈ [15,35]
```

### Abundance mode

```text
length_count_abundance(L)
=
sum(read-level count for eligible rows of length L)
```

### Unique-sequence mode

Distinct identity:

```text
sample × analysis_unit × length × strand × sequence
```

Each distinct sequence contributes 1.

### Primary comparable quantity

```text
length_fraction(L)
=
length_count(L)
/
Σ length_count(k), k=15..35
```

Rank lengths within each sample-virus unit using competition ranking with minimum rank for ties.

Across dataset, report separately by weighting mode:

- sample-balanced median length fraction;
- sample-clustered 95% CI;
- median rank;
- contributing samples/units;
- descriptive top-1/top-3 frequencies.

Current canonical abundance landscape includes:

```text
24 nt ≈ 0.535672
23 nt ≈ 0.212265
22 nt ≈ 0.080003
25 nt ≈ 0.067678
21 nt ≈ 0.031949
```

These are descriptive population properties.

---

## 01B — focused 23/24 population analysis

Retain:

```text
23 sense
23 antisense
24 sense
24 antisense
```

under both abundance and unique-sequence weighting.

Calculate:

```text
count_23
count_24
antisense_fraction_23
antisense_fraction_24
delta_antisense_fraction_24_minus_23
length23_fraction_among_23_24
length24_fraction_among_23_24
```

Interpretation limits:

Stage 01 must not conclude:

```text
23 nt = primary products
24 nt = secondary products
23 nt = Dicer
24 nt = RdRP
```

from length/strand bias alone.

Outputs:

```text
results/01_viral_23_24/
    qc/
    length_spectrum/
    fixed_23_24/
    figures/
```

---

# 02 — Length-matched terminal nucleotide enrichment

## Purpose

Measure whether physical terminal bases occur more/less often than expected from the depth-supported viral sequence opportunity available in each infection.

Primary lengths:

```text
23 nt
24 nt
```

Primary positions:

```text
5p1
5p2
3p2
3p1
```

Observed antisense reads are already sequenced 5′→3′ and must **not** be reverse-complemented again.

Expected antisense background is generated by reverse-complementing matched reference windows, then reading physical antisense termini.

Canonical observed inclusion:

```text
primary_eligible
mapping_mode = exact
virus_assignment = assigned
strand ∈ {sense, antisense}
length ∈ {23,24}
```

Matched background windows must:

- have the same length;
- lie fully within one background FASTA record;
- contain only A/C/G/T;
- never cross record boundaries.

For base `b`, position `p`:

```text
enrichment_ratio(b,p)
=
observed_fraction(b,p)
/
expected_fraction(b,p)
```

No pseudocount.

If expected frequency is zero:

```text
NA
```

Combined-strand expectation uses the observed sense/antisense mixture rather than imposing 50:50.

Primary cross-dataset result:

```text
sample-balanced median enrichment
+
sample-clustered bootstrap CI
```

Historical pair-balanced median enrichment remains available as regression/context only.

Current 23/24 terminal landscapes are highly concordant (Spearman ~0.97–0.98) and show recurring associations including:

```text
3p1 T/U enrichment
3p2 C/G enrichment
5p1 G depletion
3p1 A depletion
3p2 A depletion
```

These are empirical associations, not proof of one molecular mechanism.

---

# 03 — Official stepRNA duplex geometry

## 03.0 Purpose

Use the official stepRNA method to reconstruct complementary focal/passenger relationships and signed duplex-end geometry.

Canonical software:

```text
stepRNA 1.0.6
Bowtie2 2.5.5
```

Primary focal classes:

```text
23S
23AS
24S
24AS
```

Primary passenger-length search:

```text
15–30 nt
```

Sensitivity:

```text
18–28 nt
```

Official signed distance convention:

```text
negative = reference overhang
positive = reference underhang
0        = blunt
```

Stage 03 must retain both:

```text
marginal 5′ distance spectrum
marginal 3′ distance spectrum
full same-duplex joint (d5,d3) spectrum
```

A strong marginal distance-0 peak must not be described as evidence that most duplexes are `(0,0)` without inspecting the joint spectrum.

Pre-specified Varroa geometry of interest:

```text
(+2,-2)
```

This remains a pre-specified secondary feature inside the full landscape and must not be redefined after observing the spectrum.

Abundance weighting must use focal small-RNA abundance from canonical read-level counts, not number of passenger alignments.

Outputs:

```text
results/03_steprna/
    qc/
    provenance/
    inputs/
    raw/
    parsed/
    sensitivity/
```

Stage 03 does not directly identify a nuclease, prove cleavage, prove RdRP amplification, or create a candidate score.

---

# 04 — Sample-aware duplex-geometry aggregation and geometry-conditioned sequence features

## Purpose

Stage 04 does not rerun stepRNA.

It asks:

1. which geometry features are reproducible across samples;
2. whether focal RNAs supporting the pre-specified `(+2,-2)` geometry have terminal sequence properties beyond the general Stage 02 population.

Primary geometry aggregation is sample-balanced.

Relevant named quantities include:

```text
sample_balanced_steprna_log_ratio
sample_balanced_joint_duplex_fraction(d5,d3)
sample_balanced_joint_00_duplex_fraction
sample_balanced_varroa_2nt_joint_fraction
```

Current validated interpretation:

- marginal distance 0 can be prominent;
- fully blunt `(0,0)` same-duplex geometry is nevertheless a minority;
- pre-specified `(+2,-2)` geometry is also a minority;
- the marginal 3′ `-2` component is more reproducibly represented than the matching 5′ `+2`;
- geometry alone does not justify classifying 23 nt as primary/Dicer and 24 nt as secondary/RdRP.

Geometry-conditioned subsets:

```text
all focal
passenger-recovered
(+2,-2)-supporting focal
```

Geometry-conditioned terminal analysis reuses Stage 02 physical terminal definitions and matched expected backgrounds.

No geometry/Dicer score is carried forward to candidate ranking because geometry-conditioned effects are not sufficiently independent/consistent to justify an additional candidate-level feature.

Historical custom `Δ_Dicer` is not a canonical metric and must not be approximately reconstructed if its exact historical definition is unavailable.

---

# 05 — Viral spatial/transitivity-consistency analysis

## Purpose

Test whether viral 23/24 small-RNA spatial patterns are consistent with downstream changes expected under an amplification/transitivity-associated hypothesis.

Stage 05 is **analysis-only** and does not generate a target-candidate feature.

Canonical positive-sense analysis subset:

```text
14 samples
19 eligible positive-sense sample-virus units
BMLV
VDV-5
VDV-9
```

Key configuration:

```text
bin_size_nt = 10
windows_nt = [100,250,500]
anchor_percentile = 90
anchor_min_separation_nt = 50
minimum_anchors = 3
permutations = 5000
random_seed = 20260810
```

Anchor percentile is calculated over non-zero bins only.

Primary endpoints include:

```text
D_24AS
D_24S
delta_D_24AS_minus_24S

F24_AS_upstream
F24_AS_downstream
delta_F24_AS
```

`delta_F24_AS` is a **composition shift** within the antisense 23+24 population, not an absolute 24-nt abundance increase.

Canonical inference uses sample-balanced aggregation and sample-aware uncertainty.

Current interpretation:

- a modest downstream shift toward 24 nt appears over approximately 250–500 nt;
- no positive antisense-specific absolute 24-nt downstream directionality effect survives the pre-specified multiple-testing framework;
- therefore Stage 05 supports a **relative compositional shift**, not a demonstrated absolute secondary-siRNA wave.

Stage 05 must not retroactively assign the 23/24 populations to fixed biochemical pathways.

---

# 06 — Generic transcript target preparation and exhaustive candidate enumeration

## 06.1 Purpose

Stage 06 is target-agnostic.

It accepts a mature/spliced transcript and requested candidate lengths and enumerates every complete window.

The computational target is the transcript sequence, not a genomic interval.

Each transcript isoform is a separate target.

Canonical target registry:

```text
resources/targets/target_manifest.tsv
```

Required fields include:

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

`candidate_lengths_nt` is a parameter and must not be hard-coded.

For transcript length `L` and candidate length `w`:

```text
start_1based = 1 ... L-w+1
end_1based   = start_1based+w-1
n_candidates = L-w+1
```

Sequence definitions:

```text
target_sequence_dna
=
exact mature-transcript slice, 5′→3′

target_sequence_rna
=
target_sequence_dna with T→U

antisense_guide_sequence_rna
=
reverse complement of target_sequence_rna, 5′→3′
```

Candidate identifier:

```text
TARGET_ID__LENGTHnt__START_END
```

No accessibility, thermodynamic, empirical, geometry, transitivity, efficacy, or ranking filter is applied.

Current Vd-CHIBIN fixture:

```text
target_id = Vd_CHIBIN
transcript_id = XM_022792159.1
length = 710 nt
candidate lengths = 23,24

23 nt = 688
24 nt = 687
total = 1,375
```

Vd-CHIBIN is a regression fixture, not a hard-coded special case.

Outputs:

```text
results/06_targets/
    target_reference_summary.tsv
    target_candidates.tsv
    qc/stage06_accounting.tsv
    provenance/stage06_manifest.tsv
```

---

# 07 — Varroa empirical guide-sequence association landscape

## 07.0 Purpose

Stage 07 asks which antisense 23/24-nt sequence characteristics are associated with:

```text
representation
accumulation
```

in total Varroa viral small-RNA sequencing.

It is **not** an AGO-loading or efficacy assay.

Primary biological scope:

```text
antisense
```

Sense is retained as a comparator.

Matched expected sequence opportunity is generated from validated depth-supported viral backgrounds.

Stage 07 maintains separate endpoints for:

```text
unique representation
abundance representation
accumulation
```

and respects sample-level clustering.

## 07A — single-position landscape

Analyse every nucleotide at every guide position.

Terminal positions must exactly reproduce Stage 02:

```text
position 1   ↔ 5p1
position 2   ↔ 5p2
position L-1 ↔ 3p2
position L   ↔ 3p1
```

Failure of this regression blocks interpretation.

Internal positional discovery uses a conservative dependent-test FDR strategy.

## 07B — literature-specified validation

Pre-specified literature-guided variables include:

```text
A10
continuous GC9–14
```

These are evaluated separately from exploratory positional discovery.

## 07C — fixed-width regional GC landscape

Scan all six-nucleotide GC windows.

`GC9_14` remains literature-specified.

Other six-base windows are exploratory and use dependent-test correction.

Overlapping supported windows must be interpreted as a broad regional signal rather than independent discoveries.

## 07D — Wang/Bartel feature synthesis

Canonical named guide descriptors:

```text
W7  = A/U at g7
R10 = A/G at g10
W17 = A/U at g17
```

Stage 07 carry-forward evidence ultimately supports:

```text
A3p3
low GC / high AU at guide 3p5–10
W17
R10
```

as candidate features for later multivariable evaluation.

`W7` is not retained as a default Stage 09 predictor.

Broad early/central G-rich effects are descriptive context rather than a default score term.

## 07E — redundancy principle

Single-position and regional-GC views derive from the same sequences and may overlap strongly.

Stage 07 does not assign weights or double-count them.

Stage 07 also does not force agreement with later thermodynamic asymmetry.

Notably:

```text
A3p3 is empirically associated with accumulation
but tends to oppose classically favourable terminal asymmetry

GC_3p5_10, W17 and R10 are mostly independent of asymmetry
```

This disagreement is preserved for later layer analysis.

---

# 08 — Generic candidate biophysics

## 08.0 Purpose

Stage 08 calculates **raw predicted physical descriptors** for every Stage 06 candidate.

It is target-agnostic and independent of Stage 07 empirical feature scoring.

Stage 08 has four components:

1. whole-site predicted target accessibility;
2. target-side g2–g8 predicted accessibility;
3. guide/passenger duplex-end thermodynamic asymmetry;
4. antisense-guide self-folding.

Stage 08 does **not** rank, filter, gate, reward, penalize, or remove candidates.

## 08.1 Structural software

Canonical ViennaRNA:

```text
ViennaRNA Package 2.7.2
temperature = 37 °C
```

Primary RNAplfold parameters:

```text
W = 150
L = 100
u >= max(candidate length, 7)
```

Sensitivity:

```text
W100/L80
W200/L150
```

## 08.2 Whole-site accessibility

```text
target_whole_p_unpaired
```

Probability that the complete target interval is simultaneously unpaired in the local structural ensemble.

Higher = greater predicted intrinsic target accessibility.

Length-dependent; 23 and 24 nt are not directly interchangeable raw scales.

## 08.3 Seed-side accessibility

```text
target_seed_g2_8_p_unpaired
```

The target interval complementary to guide g2–g8.

For target end coordinate `e`:

```text
[e-7, e-1]
```

Higher = greater predicted local seed-side accessibility.

Secondary descriptor; does not supersede whole-site accessibility.

## 08.4 Terminal thermodynamic asymmetry

Canonical thermodynamic source:

```text
Zuber et al. 2022
DOI 10.1093/nar/gkac261
```

Assume a perfect complementary duplex with **no imposed Dicer overhang**.

Primary local end size:

```text
4 paired nucleotides
```

Sensitivity:

```text
5 paired nucleotides
```

For 4 bp:

```text
guide_5p_terminal_dg_4bp
=
3 NN stack increments
+
applicable correction at the real outer physical helix end

passenger_5p_terminal_dg_4bp
=
same definition
```

Do not add:

```text
isolated-duplex initiation term
symmetry term
artificial inner-end correction
Dicer overhang
```

Canonical asymmetry:

```text
asymmetry_ddg_4bp
=
guide_5p_terminal_dg_4bp
-
passenger_5p_terminal_dg_4bp
```

Interpretation:

```text
positive = guide 5′ relatively less stable
           = classical guide-favouring direction

zero     = balanced

negative = guide 5′ relatively more stable
```

`asymmetry_ddg_5bp` is sensitivity only.

## 08.5 Guide self-folding

Fold the isolated antisense guide with RNAfold-equivalent ViennaRNA MFE at 37 °C.

Record:

```text
guide_self_fold_mfe_kcal_mol
guide_self_fold_structure
```

More negative MFE = stronger predicted self-folding.

No hairpin threshold is imposed.

## 08.6 Stage 08 output

```text
results/08_candidate_biophysics/
    candidate_biophysics.tsv
    stage08_parameters.tsv
    stage08_qc.tsv
```

Current Vd-CHIBIN Stage 08 regression:

```text
23 nt = 688
24 nt = 687
total = 1,375
candidate loss = 0
```

---

# 09 — Candidate evidence layers

## 09.0 Purpose

Stage 09 converts Stage 07 empirical sequence evidence and Stage 08 raw biophysical descriptors into three explicitly separated candidate-level layers:

```text
09A Layer 1 — Varroa small-RNA accumulation propensity
09B Layer 2 — guide competence / strand-selection biophysics
09C Layer 3 — target engagement / predicted target accessibility
```

Stage 09 is an evidence-layer stage, **not** the final overall ranking stage.

It must not:

- claim measured RNAi efficacy;
- choose the final candidate;
- select the final long-dsRNA region;
- assign a final biological weight between Layers 1/2/3;
- introduce historical arbitrary cross-layer weights.

Final robust inter-layer ranking belongs to Stage 10.

---

## 09.0.1 Mandatory 23/24 separation

**23-nt and 24-nt candidates remain separate analytical strata throughout Stage 09.**

This remains true if Stage 09A selects a statistical model with shared coefficients.

Shared coefficients mean only that parameter estimation borrows information across lengths.

They do not mean:

```text
pooled validation
pooled normalization
pooled score distributions
pooled QC
pooled biological interpretation
```

The following are always separate for 23 and 24 nt:

- accounting;
- descriptive summaries;
- model performance;
- score distributions;
- percentile normalization;
- correlation/redundancy analyses;
- sensitivity analyses;
- QC.

Therefore:

```text
shared coefficients != pooled biological analysis
23 nt and 24 nt are never normalized against each other
```

No automatic 23- or 24-nt bonus is applied.

---

# 09A — Layer 1: Varroa small-RNA accumulation propensity

## 09A.1 Biological question

Stage 09A asks:

> Which antisense sequence characteristics are associated with greater accumulation in the validated Varroa viral small-RNA population, and do those associations transfer in the same direction to an unseen virus?

The primary interpretation is:

```text
Varroa small-RNA processing / recovery / accumulation propensity
```

Stage 09A must not be described as:

```text
AGO2 loading probability
RISC incorporation probability
target cleavage probability
RNAi efficacy
```

because the training libraries contain total viral small RNAs rather than purified AGO-bound guides or candidate-level knockdown measurements.

---

## 09A.2 Simplified canonical architecture

Stage 09A deliberately uses a **simple, transparent multivariable accumulation model**.

Development analyses showed that:

```text
representation
accumulation
and a representation × accumulation hurdle score
```

all contain sequence-associated signal, but accumulation propensity was the most useful primary quantity for recovering high-abundance products.

Therefore v0.20 defines:

```text
PRIMARY Layer 1
= abundance / accumulation propensity

representation
= supporting Stage 07 evidence and Stage 09A accounting/QC only

hurdle combination
= non-canonical development analysis only
```

Stage 09A does **not** fit a second representation model and does **not** fit a hurdle model.

This avoids unnecessary model duplication and avoids double-counting overlapping sequence information.

---

## 09A.3 Mandatory separate 23-nt and 24-nt models

Canonical Stage 09A fits:

```text
one 23-nt accumulation model
one 24-nt accumulation model
```

The two lengths are analysed independently.

They have:

- separate training matrices;
- separate fitted coefficients;
- separate leave-one-virus-out validation;
- separate validation summaries;
- separate candidate score distributions;
- separate candidate percentile normalization;
- separate QC.

No shared-coefficient model-selection exercise is performed in v0.20.

This choice is intentionally simple and is consistent with the project requirement that 23-nt and 24-nt products remain separate biological analyses.

No automatic bonus or penalty is applied to either length.

---

## 09A.4 Frozen training population

Reconstruct Stage 09A only from validated frozen viral inputs.

Required inclusion:

```text
primary_eligible sample × virus unit
mapping_mode = exact
virus_assignment = assigned
strand = antisense
length = 23 nt or 24 nt
sequence fully A/C/G/T
sequence corresponds to a valid depth-supported background opportunity
```

Observed antisense sequences are already in physical sequenced 5′→3′ orientation and must **not** be reverse-complemented again.

The matched theoretical antisense opportunity is the reverse complement of the corresponding valid reference window.

Current frozen-data regression fixtures:

```text
primary biological samples              20
eligible sample-virus units             54

23-nt supported opportunities       411,079
24-nt supported opportunities       408,148
total supported opportunities       819,227

23-nt represented opportunities     121,592
24-nt represented opportunities     175,564
total represented opportunities     297,156

supported observed abundance      3,445,943
```

Current represented fractions:

```text
23 nt = 29.5787%
24 nt = 43.0148%
```

Observed exact antisense sequence/group entries outside the validated background opportunity are excluded from model fitting but retained in QC.

Current audit expectation:

```text
outside-background distinct sequence/group entries = 3,616
outside-background abundance                        = 3,973
outside-background abundance fraction               ≈ 0.1152%
```

These values are regression fixtures for the frozen dataset and must be reconstructed rather than copied from historical Stage 09 outputs.

---

## 09A.5 Training unit and positive accumulation outcome

The theoretical opportunity unit is:

```text
sample × virus × length × antisense sequence
```

For the accumulation model, retain only represented supported opportunities:

```text
abundance > 0
```

For represented sequence `i`:

```text
count_i
=
summed canonical abundance
for that exact sample × virus × length × sequence
```

Primary response:

```text
y_i = ln(count_i)
```

No pseudocount is required.

The log transform reduces domination by the strongly right-skewed positive abundance distribution.

The modelling objective is **relative accumulation ranking**, not literal read-count prediction.

---

## 09A.6 Sequence predictors

Each length-specific model uses the same eight pre-specified descriptors:

```text
5p1
5p2
3p2
3p1
A3p3
GC_3p5_10
W17
R10
```

Definitions:

### `5p1`

Physical first nucleotide of the antisense guide.

### `5p2`

Physical second nucleotide of the antisense guide.

### `3p2`

Physical penultimate nucleotide of the antisense guide.

### `3p1`

Physical final nucleotide of the antisense guide.

Terminal nucleotide variables are categorical:

```text
A / C / G / U
```

with fixed regression reference category:

```text
A
```

### `A3p3`

```text
1 = third nucleotide from physical guide 3′ end is A
0 = otherwise
```

### `GC_3p5_10`

Continuous GC fraction across exactly six guide positions:

```text
3p5, 3p6, 3p7, 3p8, 3p9, 3p10
```

Equivalent sequence slice:

```text
guide[-10:-4]
```

### `W17`

```text
1 = guide position 17 is A or U
0 = guide position 17 is G or C
```

### `R10`

```text
1 = guide position 10 is A or G
0 = guide position 10 is C or U
```

Stage 02 terminal enrichment ratios are not inserted as separate predictors because the terminal categories already carry that sequence information.

Stage 03/04 geometry, Stage 05 transitivity, W7, Stage 08 accessibility, Stage 08 asymmetry and Stage 08 self-folding are excluded from Layer 1.

---

## 09A.7 Canonical model

For each length independently, fit a weighted fixed-effect linear model:

```text
ln(count_i)
=
a_g
+
βᵀX_i
+
ε_i
```

where:

```text
g
=
sample × virus group
for the current candidate length

a_g
=
unpenalized nuisance group intercept

X_i
=
the eight pre-specified sequence descriptors
after fixed categorical encoding

β
=
sequence-effect coefficients
```

No elastic-net penalty is used.

No hyperparameter search is used.

No model-structure search is used.

The canonical implementation may use weighted least squares with explicit group fixed effects or an exactly equivalent weighted within-group transformation.

The fitted nuisance group effects account for baseline differences among viral/library groups.

They are not candidate biology and are never transferred to target candidates.

The candidate-facing sequence score is:

```text
β_hatᵀX_candidate
```

without any viral nuisance intercept.

---

## 09A.8 Sample-aware fitting weights

For a given length, let:

```text
G_s
=
number of eligible sample-virus groups
contributed by biological sample s

N_g
=
number of represented training sequence species
in sample-virus group g
```

For represented training observation `i` in sample `s` and group `g`:

```text
w_i_raw
=
1 / (G_s × N_g)
```

Rescale inside the current training set:

```text
w_i
=
w_i_raw
×
n_training
/
Σ_i w_i_raw
```

so mean training weight is 1.

Therefore:

```text
each biological sample contributes equal total weight

within a sample:
each eligible virus group contributes equal total weight

within a group:
each represented sequence species contributes equal fitting weight
```

Observed read abundance is the response and is **not** also used as a fitting weight.

This weighting is recalculated within every leave-one-virus-out training set using training data only.

---

## 09A.9 Primary validation: leave one virus/family out

Stage 09A uses one simple external-style internal validation:

```text
leave-one-biological-virus/family-out
```

For each of the five current biological virus/family groups:

1. hold out all samples belonging to that virus/family;
2. fit the 23-nt model on the remaining viruses/families;
3. fit the 24-nt model on the remaining viruses/families;
4. calculate sequence-only predictions for the held-out opportunities;
5. evaluate held-out performance separately for 23 nt and 24 nt.

There is:

```text
no nested CV
no leave-one-sample-out model search
no hyperparameter tuning
no repeated model-family selection
```

Current canonical validation therefore requires approximately:

```text
5 held-out-virus fits for 23 nt
5 held-out-virus fits for 24 nt

+ 1 final full-data fit for 23 nt
+ 1 final full-data fit for 24 nt
```

for a total of approximately:

```text
12 primary accumulation fits
```

This simplicity is intentional.

---

## 09A.10 Held-out evaluation

For one held-out virus/family and one length, evaluate within each available:

```text
sample × held-out-virus × length
```

group across **all valid theoretical opportunities**.

Unrepresented opportunities receive:

```text
observed abundance = 0
```

Required per-group metrics:

```text
Spearman rho(
    predicted sequence score,
    observed abundance
)

top10_abundance_share

top10_abundance_lift
```

where:

```text
top10_abundance_share
=
abundance contained in the highest-scoring
ceil(0.10 × n) opportunities
/
total observed abundance
```

and:

```text
top10_abundance_lift
=
top10_abundance_share
/
(ceil(0.10 × n) / n)
```

Interpretation:

```text
rho > 0
=
higher predicted sequence propensity tends to
correspond to greater observed abundance

top10 lift > 1
=
predicted top-scoring opportunities capture more
abundance than random selection of the same size
```

If a held-out group has no observed abundance, its abundance-share/lift are `NA`.

---

## 09A.11 Virus-level and cross-virus summaries

Within each held-out virus/family and length:

```text
virus_holdout_median_rho
=
median of valid sample-group Spearman rho values

virus_holdout_median_top10_lift
=
median of valid sample-group top10 lifts
```

Across the five held-out virus/family analyses, report separately for 23 nt and 24 nt:

```text
median virus-holdout rho

range / individual virus-holdout rho values

number of virus holdouts with rho > 0

median virus-holdout top10 lift

range / individual virus-holdout top10 lift values

number of virus holdouts with top10 lift > 1
```

These are descriptive generalization diagnostics.

Stage 09A does **not** impose an arbitrary hard PASS threshold such as requiring a fixed number of positive folds.

The results must instead be reported transparently and interpreted according to their magnitude and consistency.

---

## 09A.12 Coefficient stability diagnostic

The five leave-one-virus-out fits already produced for validation are reused to assess coefficient stability.

For each length and fitted sequence coefficient report:

```text
final_full_data_coefficient

median_holdout_fit_coefficient

minimum_holdout_fit_coefficient

maximum_holdout_fit_coefficient

n_holdout_fits_same_sign_as_final
```

This introduces no additional model fits.

It is a diagnostic of whether an estimated sequence association is strongly dependent on one virus/family.

It is not an additional scoring layer.

---

## 09A.13 Final full-data models

After leave-one-virus-out validation, fit:

```text
one final 23-nt model
one final 24-nt model
```

using all frozen eligible viral training data for the corresponding length.

No validation result is used to tune model complexity because there is no hyperparameter/model-structure search.

The final candidate-facing raw scores are:

```text
layer1_accumulation_linear_predictor_23nt

layer1_accumulation_linear_predictor_24nt
```

implemented in the joined candidate table as the generic field:

```text
layer1_accumulation_linear_predictor
```

with interpretation determined by `candidate_length_nt`.

Higher values mean greater predicted Varroa viral-small-RNA accumulation propensity **within that fitted length-specific model**.

The raw 23-nt and 24-nt scales must not be treated as directly interchangeable.

---

## 09A.14 Candidate-facing normalization

Within each:

```text
target × candidate_length
```

convert the raw Layer-1 sequence score to the canonical favourable percentile.

Thus:

```text
23-nt candidates are ranked only against 23-nt candidates
24-nt candidates are ranked only against 24-nt candidates
```

No cross-length normalization is allowed.

No automatic 23- or 24-nt bonus is introduced.

---

## 09A.15 Representation evidence

Stage 09A does not fit a new representation model.

Representation evidence remains available from:

- Stage 07 representation analyses;
- Stage 09A training/accounting tables;
- represented fractions of the frozen theoretical opportunity universe.

This evidence supports the biological interpretation that sequence composition affects whether small-RNA products are recovered, but it is not independently added to the primary Layer-1 accumulation score.

---

## 09A.16 External validation

Muita observations must not be used for:

- fitting Stage 09A;
- predictor selection;
- modifying coefficients after inspection;
- choosing candidate scores.

The existing Muita CHH analysis remains exploratory because the true author-defined administered trigger sequence is not publicly available.

When exact Muita trigger sequences become available:

```text
freeze the Stage 09A models first
apply the appropriate length-specific model unchanged
evaluate external treatment libraries
do not refit against Muita
```

The same principle applies to future Damayo synthetic-dsRNA datasets.

# 09B — Layer 2: guide competence / strand-selection biophysics

## 09B.1 Inputs

Primary Stage 08 inputs:

```text
asymmetry_ddg_4bp
guide_self_fold_mfe_kcal_mol
```

Sensitivity-only Stage 08 input:

```text
asymmetry_ddg_5bp
```

The 4-bp asymmetry is canonical. The 5-bp asymmetry is retained only to assess sensitivity to the terminal-window definition.

## 09B.2 Favourable directions

```text
higher / more positive asymmetry_ddg_4bp
= guide 5′ end relatively less stable
= more classically guide-favouring strand-selection asymmetry

higher / less-negative guide_self_fold_mfe_kcal_mol
= weaker predicted guide self-folding
= more favourable
```

The same favourable direction applies to `asymmetry_ddg_5bp`, but it is not used in the canonical Layer-2 reference score.

## 09B.3 Separate length analyses

For each target, calculate separately:

```text
23-nt distributions
24-nt distributions
23-nt percentiles
24-nt percentiles
23-nt sensitivity correlations
24-nt sensitivity correlations
```

23 nt and 24 nt are never pooled for normalization.

## 09B.4 Normalized canonical components

Within `target × candidate_length` calculate:

```text
layer2_asymmetry_percentile
layer2_self_fold_percentile
```

using the common favourable-percentile transform.

Also calculate the sensitivity-only:

```text
layer2_asymmetry_5bp_percentile
```

from `asymmetry_ddg_5bp`.

## 09B.5 Neutral Layer-2 reference

Because no Varroa efficacy dataset identifies the correct within-layer weight:

```text
layer2_reference_score
=
0.5 × layer2_asymmetry_percentile
+
0.5 × layer2_self_fold_percentile
```

This is a **neutral reference summary**, not a learned efficacy weighting.

The 5-bp asymmetry contributes zero weight to the canonical score.

## 09B.6 4-bp / 5-bp sensitivity

Retain:

```text
layer2_asymmetry_5bp_percentile

layer2_asymmetry_4bp_5bp_percentile_difference
=
layer2_asymmetry_percentile
-
layer2_asymmetry_5bp_percentile
```

and, separately for each `target × candidate_length`:

```text
Spearman(asymmetry_ddg_4bp, asymmetry_ddg_5bp)
```

The current Vd-CHIBIN regression fixtures are approximately:

```text
23 nt rho = 0.916831
24 nt rho = 0.915505
```

No 4/5-bp disagreement gate is applied.

## 09B.7 Interpretation

Layer 2 summarizes two related but non-identical physical considerations:

```text
duplex-end strand-selection asymmetry
guide self-folding
```

It is a predicted biophysical evidence layer.

It is not an AGO-loading probability and is not calibrated to Varroa knockdown efficacy.

---
# 09C — Layer 3: target engagement / predicted target accessibility

## 09C.1 Inputs

Canonical Stage 08 RNAplfold parameterization:

```text
W150 / L100
```

Canonical inputs:

```text
target_whole_p_unpaired
target_seed_g2_8_p_unpaired
```

Sensitivity-only parameterizations:

```text
W100 / L80
W200 / L150
```

No RNAplfold calculation is rerun in Stage 09.

## 09C.2 Favourable direction

For both whole-site and seed-side accessibility:

```text
higher P(unpaired)
= more predicted target accessibility
= more favourable
```

## 09C.3 Separate length analyses

23 nt and 24 nt remain distinct.

This is particularly important for whole-site simultaneous-unpaired probabilities because candidate interval length affects the raw probability scale.

## 09C.4 Normalized canonical components

Within `target × candidate_length` calculate:

```text
layer3_whole_accessibility_percentile
layer3_seed_accessibility_percentile
```

using the common favourable-percentile transform.

## 09C.5 Neutral Layer-3 reference

```text
layer3_reference_score
=
0.5 × layer3_whole_accessibility_percentile
+
0.5 × layer3_seed_accessibility_percentile
```

This is a neutral reference summary, not a learned efficacy weighting.

## 09C.6 RNAplfold parameter sensitivity

The alternative W/L settings are robustness diagnostics only.

Within each `target × candidate_length`, report canonical-versus-sensitivity Spearman correlations for:

```text
whole-site accessibility:
    W150/L100 vs W100/L80
    W150/L100 vs W200/L150

seed-side accessibility:
    W150/L100 vs W100/L80
    W150/L100 vs W200/L150
```

Current Vd-CHIBIN regression fixtures:

```text
23 nt whole:
    W100/L80  rho = 0.895335
    W200/L150 rho = 0.926373

23 nt seed:
    W100/L80  rho = 0.913162
    W200/L150 rho = 0.947187

24 nt whole:
    W100/L80  rho = 0.900316
    W200/L150 rho = 0.924570

24 nt seed:
    W100/L80  rho = 0.913367
    W200/L150 rho = 0.947307
```

The alternative parameterizations contribute zero weight to the canonical Layer-3 reference score.

## 09C.7 Interpretation

Layer 3 represents predicted target accessibility.

It is not a direct measurement of Argonaute binding, cleavage, or RNAi efficacy.

No accessibility-based candidate filtering occurs in Stage 09.

---
# 09.1 Common favourable-percentile transform

Within:

```text
target × candidate_length
```

orient a metric so larger = more favourable.

Let:

```text
r_i = average ascending rank
n = number of candidates
```

Then:

```text
Q_i = (r_i - 0.5) / n
```

Properties:

```text
higher = more favourable
ties = average rank
n=1 -> 0.5
```

This is target-relative and not a biological probability.

23 and 24 nt are never normalized together.

---

# 09.2 Cross-layer diagnostics

Within each target × candidate_length calculate correlations among:

```text
layer1_accumulation_linear_predictor
asymmetry_ddg_4bp
guide_self_fold_mfe_kcal_mol
target_whole_p_unpaired
target_seed_g2_8_p_unpaired
```

and normalized forms.

Purpose:

- redundancy;
- independence;
- antagonism;
- candidate-level disagreement.

Do not force distinct evidence layers to agree.

---

# 09.3 Stage 09 stops before cross-layer ranking

Stage 09 must not define:

```text
overall_stage09_score
efficacy_probability
final_candidate_rank
best_candidate
```

Stage 09 produces the three evidence layers.

Stage 10A is the first stage permitted to integrate those layers into an individual-window ranking.

Stage 10A must still avoid claims of calibrated RNAi efficacy because no measured Varroa efficacy outcome currently identifies empirical cross-layer weights.

---

# 09.4 Canonical outputs

```text
results/09_feature_layers/
│
├── 09A_layer1_accumulation/
│   ├── layer1_training_accounting.tsv
│   ├── layer1_coefficients_23nt.tsv
│   ├── layer1_coefficients_24nt.tsv
│   ├── layer1_leave_one_virus_out.tsv
│   ├── layer1_cv_summary_23nt.tsv
│   ├── layer1_cv_summary_24nt.tsv
│   ├── layer1_coefficient_stability.tsv
│   ├── layer1_model_provenance.tsv
│   └── candidate_layer1.tsv
│
├── 09B_layer2_guide_competence/
│   ├── candidate_layer2.tsv
│   └── layer2_correlations.tsv
│
├── 09C_layer3_target_engagement/
│   ├── candidate_layer3.tsv
│   └── layer3_correlations.tsv
│
├── candidate_stage09_layers.tsv
├── stage09_feature_correlations.tsv
├── stage09_parameters.tsv
└── stage09_qc.tsv
```

The joined candidate table may contain multiple candidate lengths, but every specified normalization/statistical analysis remains length-stratified.

---

# 10A — Individual-window evidence integration and ranking

## 10A.1 Purpose

Stage 10A asks:

> Which individual candidate windows are consistently favourable across empirical Varroa accumulation, guide-competence biophysics, and predicted target accessibility?

Stage 10A integrates the three Stage 09 evidence layers for each individual candidate window.

It does **not**:

- estimate a probability of RNAi efficacy;
- select a long dsRNA region;
- simulate construct architecture;
- apply junction penalties;
- fit a new biological model.

The current project lacks sufficiently large, independent, measured Varroa candidate-level efficacy data from which defensible cross-layer coefficients could be learned.

Therefore Stage 10A uses:

1. a transparent **equal-layer reference score** for a complete practical ranking;
2. a **Pareto-front analysis** as a complementary weight-free multi-objective diagnostic;
3. a **minimum-layer score** as a weakest-layer diagnostic only.

This is intentionally different from transferring coefficients from another species. For example, dsRIP could fit feature contributions against experimentally measured insecticidal outcomes in *Tribolium*; those coefficients are not transferred into the Varroa model.

---

## 10A.2 Inputs

Required Stage 09 inputs for each candidate:

```text
layer1_accumulation_percentile
layer2_reference_score
layer3_reference_score
```

Current Vd-CHIBIN regression fixture:

```text
23 nt candidates = 688
24 nt candidates = 687
total            = 1,375
```

Candidate identity, target coordinates, guide sequence and target sequence must be preserved exactly.

No raw Stage 08 metric enters Stage 10A directly.

Stage 08 contributes only through the already defined Stage 09B and Stage 09C evidence layers.

---

## 10A.3 Equalize the three layer-level integration scales

`layer1_accumulation_percentile` is already a favourable within-target, within-length percentile.

However:

```text
layer2_reference_score
layer3_reference_score
```

are averages of favourable component percentiles and therefore need not have the same empirical distribution as Layer 1.

Before equal-layer integration, define:

```text
stage10_layer1_percentile
=
layer1_accumulation_percentile
```

and calculate, within each:

```text
target × candidate_length
```

```text
stage10_layer2_percentile
=
favourable percentile of layer2_reference_score

stage10_layer3_percentile
=
favourable percentile of layer3_reference_score
```

using the canonical percentile transform:

```text
Q_i = (r_i - 0.5) / n
```

where `r_i` is the average ascending rank and higher values are more favourable.

This step places all three **layer-level summaries on the same rank-percentile scale** before assigning equal layer weights.

It is a rank normalization, not a probability calibration.

23 nt and 24 nt are never pooled.

---

## 10A.4 Canonical equal-layer reference score

For candidate `i`, define:

```text
L1_i = stage10_layer1_percentile
L2_i = stage10_layer2_percentile
L3_i = stage10_layer3_percentile
```

Canonical Stage 10A score:

```text
stage10_equal_layer_score
=
(L1_i + L2_i + L3_i) / 3
```

Thus:

```text
Layer 1 weight = 1/3
Layer 2 weight = 1/3
Layer 3 weight = 1/3
```

The equal thirds are not claimed to be empirically optimal.

They are a neutral integration rule chosen because no Varroa efficacy outcome currently supports learned cross-layer weights.

Because Stage 09B and Stage 09C themselves use equal two-component reference summaries, the canonical arithmetic can be described as:

```text
Layer 1 accumulation                         = 1/3

within Layer 2:
    4-bp thermodynamic asymmetry             = 1/6
    guide self-folding                       = 1/6

within Layer 3:
    whole-site accessibility                 = 1/6
    seed-side accessibility                  = 1/6
```

This decomposition is descriptive arithmetic only; Stage 10A still treats L1, L2 and L3 as the three conceptual evidence layers.

The 5-bp asymmetry and alternative RNAplfold W/L parameterizations remain sensitivity diagnostics and contribute zero weight to the canonical Stage 10A score.

---

## 10A.5 Canonical complete ranking

Within each:

```text
target × candidate_length
```

calculate:

```text
stage10_equal_layer_rank
```

as the average descending rank of `stage10_equal_layer_score`:

```text
rank 1 = most favourable
```

Ties receive average rank.

Do not break equal-score ties by transcript coordinate, candidate ID, or arbitrary ordering.

Also calculate:

```text
stage10_equal_layer_percentile
```

by applying the common favourable-percentile transform directly to `stage10_equal_layer_score`.

Interpretation:

```text
higher stage10_equal_layer_score       = more favourable neutral integrated evidence
lower  stage10_equal_layer_rank        = better rank
higher stage10_equal_layer_percentile  = better relative position
```

The equal-layer rank is the **primary practical complete ranking** in Stage 10A.

It is not a calibrated efficacy ranking.

---

## 10A.6 Pareto dominance

For candidates within the same:

```text
target × candidate_length
```

compare the three layer percentiles:

```text
(L1, L2, L3)
```

Candidate `A` dominates candidate `B` if and only if:

```text
L1_A >= L1_B
L2_A >= L2_B
L3_A >= L3_B
```

and at least one inequality is strict.

Equivalent mathematical definition:

```text
A dominates B
iff
A is no worse than B in every layer
and strictly better in at least one layer
```

Identical three-layer vectors do not dominate one another.

Use the full stored numerical values; do not round before dominance comparisons.

Pareto dominance is evaluated separately for each candidate length.

---

## 10A.7 Pareto fronts

Assign iterative non-dominated fronts:

```text
front 1
=
all candidates not dominated by any other candidate

front 2
=
non-dominated candidates after removing front 1

front 3
=
non-dominated candidates after removing fronts 1 and 2

...
```

Output:

```text
stage10_pareto_front
```

with:

```text
1 = first / non-dominated front
2 = second front
3 = third front
...
```

Lower front number indicates a stronger multi-objective position.

Every candidate must receive exactly one front.

Pareto front is **not** a complete ranking because many candidates may occupy the same front and represent different trade-offs.

Therefore:

```text
equal-layer rank
=
primary practical complete ranking

Pareto front
=
complementary multi-objective diagnostic
```

Do not combine Pareto front number with the equal-layer score to create another score.

---

## 10A.8 Minimum-layer diagnostic

For each candidate:

```text
stage10_minimum_layer_score
=
min(L1, L2, L3)
```

Higher values indicate a stronger weakest evidence layer.

Purpose:

- expose candidates whose high equal-layer average masks one weak layer;
- provide an easily interpreted balance diagnostic.

It is not:

- a gate;
- a filter;
- a tie-break;
- the canonical ranking score.

---

## 10A.9 Layer correlation diagnostics

Within each:

```text
target × candidate_length
```

report Spearman correlations for:

```text
L1 vs L2
L1 vs L3
L2 vs L3
```

Purpose:

- identify redundancy;
- identify complementary information;
- identify antagonistic layer behaviour.

These correlations do not modify ranking.

---

## 10A.10 Pareto summary diagnostics

For each:

```text
target × candidate_length × pareto_front
```

report:

```text
n_candidates
fraction_candidates
```

This describes how densely candidates occupy successive fronts.

It does not alter candidate ranking.

---

## 10A.11 Separate 23-nt and 24-nt ranking populations

For the current Varroa model:

```text
23 nt
24 nt
```

remain completely separate Stage 10A ranking populations.

Never calculate:

```text
one combined 23/24 percentile
one combined 23/24 rank
one Pareto front containing both lengths
```

No automatic 23-nt or 24-nt bonus is introduced.

Comparison/integration across lengths, if required for longer-region design, belongs to a later explicitly specified stage.

---

## 10A.12 Canonical outputs

```text
results/10_candidate_integration/
├── candidate_stage10a.tsv
├── stage10a_layer_correlations.tsv
├── stage10a_pareto_summary.tsv
├── stage10_parameters.tsv
└── stage10_qc.tsv
```

Minimum candidate-facing columns:

```text
candidate_id
target_id
candidate_length_nt
target_start_1based
target_end_1based
target_sequence
guide_sequence

layer1_accumulation_percentile
layer2_reference_score
layer3_reference_score

stage10_layer1_percentile
stage10_layer2_percentile
stage10_layer3_percentile

stage10_equal_layer_score
stage10_equal_layer_rank
stage10_equal_layer_percentile

stage10_pareto_front
stage10_minimum_layer_score
```

Existing canonical identifier/coordinate field names should be reused exactly where they differ from the display names above.

---

## 10A.13 Required QC

Required checks:

```text
current Vd-CHIBIN:
    23 nt = 688
    24 nt = 687
    total = 1,375

candidate IDs exactly preserved
candidate coordinates/sequences exactly preserved

no NA/Inf in required Stage 10A numeric fields

stage10_layer1_percentile
exactly equals the Stage 09A accumulation percentile

Stage 10 L2/L3 percentiles
use the exact favourable-percentile formula

23/24 normalization completely separate

equal-layer score arithmetic exact to machine precision

equal-layer rank tie handling exact

Pareto fronts exhaustive and mutually exclusive

identical vectors do not dominate one another

minimum-layer arithmetic exact

no candidate filtering or gating

no Stage 08 raw metric direct leakage

no region aggregation

no long-dsRNA / construct scoring

no efficacy-probability output
```

---

## 10A.14 Efficiency and workflow safety

Stage 10A is a small deterministic table transformation.

It must not:

- fit a model;
- perform cross-validation;
- optimize weights;
- rerun Stage 09A;
- rerun RNAfold/RNAplfold/ViennaRNA;
- rerun any Stage 00–09 upstream computation.

Implementation should:

1. read the required Stage 09 candidate outputs once;
2. join candidates by exact candidate identifier with one-to-one assertions;
3. calculate percentile transforms vectorially;
4. calculate equal-layer scores/ranks vectorially;
5. calculate Pareto fronts deterministically within each target-length stratum;
6. write the canonical outputs.

With approximately 1,375 current candidates, expected runtime should be seconds.

A dry run of the exact Stage 10A target must show Stage 10A only.

If an upstream stage is unexpectedly scheduled, stop before execution.

---

## 10A.15 Scientific interpretation and limitation

The canonical Stage 10A equal-layer score is:

```text
Project-specific
Canonical
Neutral integration
```

It is not:

```text
predicted knockdown percentage
RNAi efficacy probability
mortality probability
AGO-loading probability
```

A high score means:

> the candidate ranks favourably, on average, across the three predefined evidence layers under an explicitly equal-layer integration rule.

The Pareto front answers a different question:

> is another candidate at least as favourable in every layer and strictly more favourable in at least one?

These two views are intentionally retained side by side.

If sufficiently large, independent Varroa efficacy datasets become available in the future, empirical cross-layer coefficients may be learned against measured outcomes using a separately specified and properly validated model.

Until then:

- do not fit weights to the current target candidate rankings;
- do not borrow numerical coefficients from dsRIP or another species;
- retain the equal-layer score as a transparent reference.

---

# 11. Reserved future stage — region selection and visualization

Stage 11 is deliberately not implemented or fully frozen in v0.21.

Its intended question is:

> Given the individual-window rankings from Stage 10A, which longer target regions contain the most favourable population of candidate windows?

Expected user-facing parameters include:

```text
desired_region_length_nt
selection_direction = highest / lowest
number_of_regions_to_return
```

Possible future transcript-region restrictions may also be supported.

The current design intention is to summarize contained Stage 10 window scores using the **mean** as the primary region-level statistic, consistent with the project's earlier sliding-region analyses and conceptually similar to long-region aggregation used by dsRIP.

However, the exact treatment of:

- 23-nt versus 24-nt information within a region;
- overlapping selected regions;
- lowest-score negative-control regions;
- region boundary rules;
- graphical summaries;
- construct architecture;

must be specified before Stage 11 implementation.

Stage 11 must consume Stage 10A outputs rather than recompute Stage 08–10 biology.

---

# 12. Future website / app extensibility principle

The current empirical Varroa model is scientifically characterized primarily for:

```text
23 nt
24 nt
```

The eventual software/website must nevertheless be **parameter-driven**, not hard-coded to these lengths.

A future user may request:

```text
candidate_length = L nt
```

For Varroa, the interface should recommend 23/24 nt because these are currently best supported by empirical data, but it should not prevent another choice such as 20 or 21 nt.

The software should then run every scientifically applicable analysis for the requested length.

Outputs should distinguish:

```text
validated empirical model
generic mechanistic / biophysical calculation
extrapolated model
unsupported / unavailable model
```

For example:

- candidate enumeration can accept arbitrary positive requested lengths that fit within the transcript;
- accessibility and folding calculations can be parameterized for other lengths subject to method constraints;
- the current Stage 09A empirical accumulation model must not be silently described as validated at an untrained length;
- Stage 10A can integrate a requested length only when the required Stage 09 layer inputs exist and their evidence status is reported.

If multiple user-requested lengths are analysed, each length remains its own normalization/analysis stratum unless a future validated specification explicitly says otherwise.

Recommended settings should guide users without unnecessarily restricting them.

---

# 13. Reproducibility parameters

Configuration/provenance must record relevant values including:

```text
target_lengths / candidate_lengths
steprna_passenger_range
steprna_sensitivity_range
transitivity_bin_size_nt
transitivity_windows_nt
transitivity_anchor_percentile
transitivity_anchor_min_separation_nt
transitivity_min_anchors
transitivity_permutations
bootstrap_replicates
random_seed
ViennaRNA version
RNAplfold parameter sets
thermodynamic parameter resource/version
Stage 09A length-specific model family
Stage 09A predictor encoding
Stage 09A sample-weighting rule
Stage 09A nuisance fixed-effect definition
Stage 09A leave-one-virus/family-out validation grouping
Stage 10A layer-level percentile transform
Stage 10A equal-layer weights = 1/3,1/3,1/3
Stage 10A Pareto dominance definition
Stage 10A ranking stratum = target × candidate_length
```

Each run records:

- Git commit;
- configuration;
- software versions;
- frozen-core path;
- input identity/checksums where practical;
- run date;
- random seed.

---

# 14. Required deterministic tests

## Stage 01

- exact/assigned inclusion;
- abundance uses `count`;
- true sequence deduplication;
- length fractions sum correctly;
- competition rank ties;
- 23/24 strand fractions;
- zero denominators -> `NA`.

## Stage 02

- physical 5p1/5p2/3p2/3p1 extraction;
- observed antisense not reverse-complemented;
- expected antisense reverse-complemented correctly;
- strand-weighted combined expectation;
- zero expected frequency -> `NA`.

## Stage 03/04

- correct focal class;
- opposite-strand passenger selection;
- passenger-length filters;
- official stepRNA signed-distance convention;
- parser correctness;
- passenger recovery denominator;
- joint geometry fractions sum to 1;
- `(+2,-2)` regression;
- sample-balanced aggregation.

## Stage 05

- alignment midpoint;
- 10-nt bin assignment;
- multimapping weight conservation;
- anchor percentile/separation;
- window boundaries;
- exact D and F24 formulas;
- permutation shift logic;
- nonzero Monte Carlo P-value formula;
- sample-level aggregation;
- multiple-testing family;
- fixed-seed reproducibility.

## Stage 06

- transcript hash/length;
- exhaustive `L-w+1` enumeration;
- first/last interval;
- all legal starts exactly once;
- target slice exact;
- RNA conversion exact;
- guide reverse complement exact;
- candidate IDs unique;
- no candidate filtering.

## Stage 07

- matched background construction;
- terminal regression to Stage 02;
- positional feature extraction;
- fixed 6-nt GC windows;
- A10/GC9–14 validation family;
- BY/BH multiple-testing definitions;
- sample-aware aggregation;
- no efficacy/ranking output.

## Stage 08

- Stage 06 row preservation;
- RNAplfold whole-site interval;
- target g2–g8 interval;
- primary/sensitivity parameters;
- exact 4-bp thermodynamic sign convention;
- 5-bp sensitivity;
- no Dicer overhang;
- RNAfold guide MFE;
- no scoring/filtering/ranking.

## Stage 09

- frozen Layer-1 accounting fixtures;
- exact guide orientation;
- exact eight predictors;
- no Stage 08 leakage into Layer 1;
- exact sample-aware weight calculation;
- weighted fixed-effect regression equivalence on synthetic data;
- five leave-one-virus/family-out folds only;
- no nested CV or hyperparameter search;
- held-out virus/family leakage prevention;
- 23/24 models fitted separately;
- 23/24 validation summaries separate;
- 23/24 normalization separate;
- coefficient-stability summary reuses existing holdout fits;
- representation not refitted or double-counted;
- Layer-2 favourable directions;
- 4-bp asymmetry canonical and 5-bp asymmetry sensitivity-only;
- Layer-2 neutral 0.5 reference exact;
- Layer-3 favourable directions;
- W150/L100 canonical and alternative W/L settings sensitivity-only;
- Layer-3 neutral 0.5 reference exact;
- Stage 06/08/09 row counts identical;
- no Stage 09 final overall rank.

## Stage 10A

- Stage 09 candidate row/ID preservation;
- Stage 10 Layer-1 percentile exact-copy regression;
- Stage 10 Layer-2/Layer-3 favourable-percentile transforms exact;
- average-rank tie handling;
- `n=1 -> 0.5`;
- 23/24 ranking populations separate;
- equal-layer arithmetic exactly `1/3 + 1/3 + 1/3`;
- equal-score ties use average rank;
- equal-layer percentile formula exact;
- Pareto dominance synthetic examples;
- identical vectors do not dominate one another;
- every candidate assigned exactly one Pareto front;
- minimum-layer score exact;
- no filtering/gating;
- no direct Stage 08 metric in Stage 10 integration;
- no region aggregation or Stage 11 logic;
- no efficacy-probability output.

---

# 15. Explicitly superseded / excluded analyses

Historical downstream results remain references only.

Not canonical inputs:

- legacy fixed-length summaries;
- legacy custom Dicer-overhang summaries;
- historical approximate Dicer scores;
- v1.4.0 transitivity implementation;
- old arbitrary 0.6/0.3/0.1 candidate weighting;
- historical Stage 09/selection outputs if any.

Currently outside scope:

- host transitivity;
- formal CHH/Pero Muita validation before true trigger release;
- empirically fitted Stage 10 cross-layer efficacy weights;
- Stage 11 Vd-CHIBIN region selection and visualization;
- construct concatenation/junction scoring;
- full user-facing Nectar Designer implementation.

---

# 16. Definition of current pipeline success

The canonical build is successful when:

1. a fresh clone can point to the frozen validated core;
2. Stage 00 validates without modifying it;
3. Stages 01–05 regenerate their canonical viral analyses from frozen inputs;
4. Stage 06 generically enumerates requested transcript candidate lengths;
5. Stage 07 reconstructs matched-background empirical guide-sequence associations;
6. Stage 08 preserves every candidate while computing raw biophysics;
7. Stage 09A learns two transparent length-specific empirical accumulation models using the eight specified predictors;
8. 23 nt and 24 nt remain separate analyses for fitting, validation and normalization;
9. Stage 09B and 09C generate length-stratified normalized evidence layers;
10. Stage 10A preserves all candidates and produces the canonical equal-layer individual-window ranking;
11. Stage 10A also reports Pareto fronts and minimum-layer diagnostics without changing the primary ranking;
12. no Muita/Damayo external result influences training or cross-layer weights;
13. no borrowed or efficacy-fitted cross-layer coefficient is introduced without suitable Varroa outcome data;
14. all important metrics are documented in `docs/METRIC_DICTIONARY.md`;
15. all major calculations have deterministic tests;
16. parameters/provenance are machine-readable;
17. no upstream frozen input is modified;
18. Stage 11 region selection remains separate until explicitly specified.

---
---

# 17. Methodological references

- Murcott B, Pawluk RJ, Protasio AV, Akinmusola RY, Lastik D, Hunt VL. 2022. *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. Frontiers in Bioinformatics 2:994871. DOI: `10.3389/fbinf.2022.994871`.
- Benjamini Y, Hochberg Y. 1995. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. JRSS B 57:289–300.
- Benjamini Y, Yekutieli D. 2001. *The control of the false discovery rate in multiple testing under dependency*. Annals of Statistics 29:1165–1188.
- Phipson B, Smyth GK. 2010. *Permutation P-values Should Never Be Zero*. Statistical Applications in Genetics and Molecular Biology 9:Article 39.
- Schwarz DS et al. 2003. *Asymmetry in the assembly of the RNAi enzyme complex.* Cell 115:199–208. DOI: `10.1016/S0092-8674(03)00759-1`.
- Tomari Y et al. 2004. *A protein sensor for siRNA asymmetry.* Science 306:1377–1380. DOI: `10.1126/science.1102755`.
- Zuber J et al. 2022. *Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects.* Nucleic Acids Research 50:5251–5262. DOI: `10.1093/nar/gkac261`.
- Wang PY, Bartel DP. 2024. *The guide-RNA sequence dictates the slicing kinetics and conformational dynamics of the Argonaute silencing complex.* Molecular Cell 84:2918–2934.e11. DOI: `10.1016/j.molcel.2024.06.026`.
- Ruijtenberg S et al. 2020. *mRNA structural dynamics shape Argonaute-target interactions.* Nature Structural & Molecular Biology 27:790–801. DOI: `10.1038/s41594-020-0461-1`.
- Cedden D, Güney G, Rostás M, Bucher G. 2025. *Optimizing dsRNA sequences for RNAi in pest control and research with the dsRIP web platform.* BMC Biology 23:114. DOI: `10.1186/s12915-025-02219-6`.

Evidence from non-Varroa systems is used as mechanistic/comparative support only. Their numerical coefficients are not transferred into the canonical Varroa empirical model.
