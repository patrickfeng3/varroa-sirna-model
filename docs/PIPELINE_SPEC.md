# Canonical Varroa vsiRNA Pipeline Specification

**Specification version:** 0.19.1  
**Status:** Stages 00–08 implemented/validated as previously frozen; Stage 09A/09B/09C implementation-ready specification  
**Scope:** Canonical viral small-RNA analysis, generic transcript candidate enumeration, empirical Varroa guide-sequence association, candidate biophysics, and Stage 09 three-layer candidate evidence synthesis  
**Host transitivity:** Excluded from current canonical build  
**Final candidate / construct ranking:** Deferred to Stage 10

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

10  future robust inter-layer ranking
    + long-dsRNA region / construct selection
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

> Which antisense sequence characteristics predict greater accumulation in the validated Varroa viral small-RNA population?

Interpretation:

```text
Varroa small-RNA processing / recovery / accumulation propensity
```

Do not call it AGO loading or efficacy.

## 09A.2 Architectural decision

Development-stage comparisons evaluated:

```text
representation only
accumulation only
hurdle = representation × conditional accumulation
```

Accumulation propensity gave the strongest overall held-out abundance recovery.

Therefore:

```text
PRIMARY Layer 1      = accumulation propensity
SECONDARY DIAGNOSTIC = representation probability
BENCHMARK ONLY        = hurdle combination
```

Representation must not be added as an independent equal-weight second Layer-1 contribution.

## 09A.3 Frozen training population

Inclusion:

```text
primary_eligible sample-virus
mapping_mode = exact
virus_assignment = assigned
strand = antisense
length ∈ {23,24}
fully A/C/G/T
valid depth-supported background opportunity
```

Current frozen fixtures:

```text
20 primary samples
54 sample-virus units
108 sample-virus-length groups

23-nt opportunities = 411,079
24-nt opportunities = 408,148
total opportunities = 819,227

23-nt represented = 121,592
24-nt represented = 175,564
total represented = 297,156

supported abundance = 3,445,943
```

Outside-background exact antisense observations are excluded from model fitting but retained in QC.

## 09A.4 Training unit

```text
sample × virus × length × antisense sequence
```

## 09A.5 Primary outcome

Among represented supported sequences:

```text
count_i = summed canonical abundance
y_i = ln(count_i)
```

No pseudocount.

The objective is ranking relative accumulation propensity.

## 09A.6 Representation diagnostic outcome

Across all valid opportunities:

```text
Y_rep = 1 if exact sequence observed
        0 otherwise
```

## 09A.7 Predictors

Use exactly:

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

```text
A3p3 = A at third nucleotide from physical guide 3′ end

GC_3p5_10
= GC fraction across exactly six guide bases 3p5–3p10

W17 = A/U at g17

R10 = A/G at g10
```

Do not add Stage 02 enrichment ratios as separate predictors.

Do not include Stage 03/04 geometry, Stage 05 transitivity, W7, or Stage 08 variables in Layer 1.

## 09A.8 Canonical accumulation model

The canonical Stage 09A accumulation model is a **fixed-effect elastic-net regression of positive log abundance**.

For represented sequence `i` in training group `g`:

```text
y_i = ln(count_i)

y_i
=
a_g
+
βᵀX_i
+
ε_i
```

where:

```text
a_g
= unpenalized nuisance intercept for
  sample × virus × candidate_length

X_i
= pre-specified sequence-feature vector

β
= penalized sequence-effect coefficients
```

The nuisance intercepts account for baseline differences between viral/library/length groups.

They are not candidate biology and are never transferred to new target candidates.

### Exact nuisance-intercept treatment

For the linear accumulation model, nuisance intercepts are removed by the fixed-effect within-group transformation inside each training fold.

For every training group `g`:

```text
y*_i = y_i - weighted_mean_g(y)

X*_ij = X_ij - weighted_mean_g(X_j)
```

The elastic-net model is then fitted to:

```text
y* ~ X*
```

with:

```text
fit_intercept = false
```

This is the canonical implementation of unpenalized group fixed effects.

Held-out groups are never used to calculate these training-group means.

For candidate application, the frozen sequence model is evaluated from the original candidate sequence features using the stored training-fold/final-model encoding and scaling parameters; no viral group intercept is added.

## 09A.9 Predictor encoding and scaling

Terminal nucleotide predictors are one-hot encoded with fixed reference category:

```text
A
```

Binary features:

```text
A3p3
W17
R10
```

use `0/1` coding.

`GC_3p5_10` remains continuous.

For penalized fitting, every encoded sequence-feature column is standardized using **training data only**:

```text
z_j = (x_j - mean_training_j) / sd_training_j
```

Rules:

- scaling parameters are learned inside each inner/outer training set only;
- held-out data never contribute to means or SDs;
- columns with zero training SD are removed for that fold and recorded in QC;
- final full-data scaling parameters are frozen and stored;
- candidate application uses the final frozen scaling parameters.

After standardization, the fixed-effect within-group transformation is applied.

## 09A.10 Exact elastic-net objective

For training observations with sample-aware weights `w_i`, fit sequence coefficients by minimizing:

```text
(1 / (2 Σ_i w_i))
Σ_i w_i (y*_i - X*_i β)^2

+
alpha [
    l1_ratio × ||β||_1
    +
    0.5 × (1-l1_ratio) × ||β||_2^2
]
```

Canonical hyperparameter grid:

```text
alpha ∈ {
    1e-5,
    3e-5,
    1e-4,
    3e-4,
    1e-3,
    3e-3,
    1e-2,
    3e-2,
    1e-1,
    3e-1,
    1
}

l1_ratio ∈ {
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    1.00
}
```

If the selected `alpha` is on the minimum or maximum grid boundary in the final tuning procedure, record:

```text
hyperparameter_boundary_warning = true
```

and inspect whether a wider grid is required before freezing coefficients.

The canonical model is the log-abundance elastic-net model.

A count-native overdispersed positive-count model may be explored later as a **non-blocking sensitivity analysis**, but it does not replace the canonical model in v0.19.1 unless the specification is deliberately revised before inspecting candidate rankings.

## 09A.11 Exact sample-aware fitting weights

Let:

```text
s(i) = biological sample containing observation i
g(i) = sample × virus × length group containing observation i

G_s = number of eligible sample-virus-length groups
      contributed by biological sample s
      in the current training set

N_g = number of represented training sequence species
      in group g
```

Define the unscaled observation weight:

```text
w_i_raw
=
1 / (G_s(i) × N_g(i))
```

Then rescale within the current training set:

```text
w_i
=
w_i_raw
×
n_training
/
Σ_i w_i_raw
```

so that mean training weight is 1.

Consequences:

```text
each biological sample contributes equal total weight

within a sample:
each eligible sample-virus-length group contributes
equal total weight

within a group:
each represented sequence species contributes
equal fitting weight
```

Observed read abundance is the response and is **not** also used as a model-fitting weight.

Sensitivity analysis:

```text
equal_group_weighting
```

where every sample-virus-length group contributes equal total weight regardless of how many eligible groups occur in a sample.

## 09A.12 Candidate 23/24 model structures

Three pre-specified structures are compared.

### Structure A — shared effects

```text
one β coefficient set shared by 23 nt and 24 nt
```

Length-specific baseline differences remain absorbed by the nuisance sample-virus-length fixed effects.

### Structure B — shared effects plus length interactions

```text
shared main effects
+
feature × candidate_length interactions
```

Candidate length is encoded only as the 23/24 interaction indicator required for this comparison.

### Structure C — separate models

```text
independent 23-nt model
independent 24-nt model
```

No structure is declared correct in advance.

Regardless of the selected structure:

```text
23-nt validation remains separate
24-nt validation remains separate

23-nt candidate normalization remains separate
24-nt candidate normalization remains separate
```

Shared coefficients never imply pooled biological analysis.

## 09A.13 Nested cross-validation

### Outer primary validation

Use:

```text
leave-one-biological-virus/family-out
```

All units belonging to the held-out virus/family are excluded from:

- sequence-feature scaling;
- fixed-effect centering;
- hyperparameter tuning;
- model-structure choice;
- model fitting.

Where multiple analysis units are recognized as closely related strains of one biological virus, they must be held out together.

### Inner tuning

Within each outer training set, use leave-one-virus/family-out CV among the remaining training viruses/families.

For every candidate model structure and hyperparameter configuration:

1. fit using only the inner-training viruses/families;
2. predict the inner-held-out virus/family;
3. calculate performance separately for 23 nt and 24 nt.

### Secondary validation

After primary virus/family CV, also report:

```text
leave-one-biological-sample-out
```

using the already specified model-selection procedure.

This is a robustness analysis and does not replace the primary virus/family transfer test.

## 09A.14 Exact model-selection objective

For every inner-held-out fold, calculate performance separately for 23 nt and 24 nt.

For each candidate configuration calculate:

```text
M23 = median within-group Spearman rho for 23 nt
M24 = median within-group Spearman rho for 24 nt
```

Primary selection score:

```text
selection_score_rho
=
(M23 + M24) / 2
```

Thus 23 nt and 24 nt contribute equal weight to model selection even if their opportunity counts differ.

Secondary tie-break:

```text
L23 = median top10 abundance lift for 23 nt
L24 = median top10 abundance lift for 24 nt

selection_score_top10
=
(L23 + L24) / 2
```

Selection order:

1. maximize `selection_score_rho`;
2. if tied to numerical tolerance `1e-6`, maximize `selection_score_top10`;
3. if still tied, choose the larger `alpha` (stronger regularization);
4. if still tied, choose the larger `l1_ratio`;
5. if still tied, prefer the simpler structure in order:

```text
A shared
B shared + interactions
C separate
```

A candidate structure/configuration must not be selected solely because one length performs strongly while the other is substantially degraded without this being visible in the separate 23/24 outputs.

All separate 23/24 metrics are retained regardless of the combined selection statistic.

## 09A.15 Held-out prediction metrics

Evaluate predictions within each held-out:

```text
sample × virus × candidate_length
```

group across **all valid theoretical opportunities**, assigning:

```text
observed abundance = 0
```

to unrepresented opportunities.

Required:

```text
Spearman rho(predicted score, observed abundance)

top10_abundance_share

top10_abundance_lift

conditional_positive_spearman_rho
```

where:

```text
top10_abundance_share
=
abundance in highest-scoring ceil(0.10 × n) opportunities
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

If total observed abundance is zero, abundance-share/lift are `NA`.

All validation summaries are reported independently for 23 nt and 24 nt.

## 09A.16 Representation diagnostic model

The representation model remains secondary.

Response:

```text
Y_rep = 1 if exact supported opportunity is represented
        0 otherwise
```

Use the same eight sequence predictors.

The model is a regularized logistic regression with:

```text
sample × virus × length nuisance intercepts unpenalized
sequence-effect coefficients penalized
```

The implementation must use a solver that supports coefficient-specific penalty factors, with:

```text
penalty_factor(group nuisance intercepts) = 0
penalty_factor(sequence effects)          = 1
```

or an exactly equivalent implementation.

Representation tuning must occur inside the same leakage-free grouped CV framework.

Required diagnostics, separately for 23 nt and 24 nt:

```text
ROC-AUC
average precision
AP lift relative to prevalence
top-decile representation enrichment
```

The representation prediction is not independently added to the primary Layer-1 accumulation score.

## 09A.17 Hurdle benchmark

For diagnostic comparison only:

```text
hurdle_score
=
P(represented)
×
conditional accumulation prediction
```

Evaluate under the same held-out scheme.

The hurdle result is retained as a benchmark and does not become canonical without a future pre-specified revision.

## 09A.18 Final frozen Layer-1 fit

After completing nested outer evaluation:

1. choose the model structure using the pre-specified CV rule;
2. using all frozen viral training data, tune `alpha` and `l1_ratio` with leave-one-virus/family-out CV;
3. fit the final model on all frozen viral training data;
4. freeze sequence encoding, scaling parameters, structure, hyperparameters and coefficients;
5. store exact software/package versions.

Candidate-facing raw score:

```text
layer1_accumulation_linear_predictor
=
β_hatᵀ Z_candidate
```

where `Z_candidate` is the candidate sequence-feature vector after applying the frozen final encoding/scaling.

No viral nuisance intercept is added.

Higher = greater predicted Varroa viral-small-RNA accumulation propensity.

This is not a predicted read count.

## 09A.19 Candidate-facing normalization

Normalize the final raw Layer-1 prediction separately within:

```text
target × candidate_length
```

to obtain:

```text
layer1_accumulation_percentile
```

23 nt and 24 nt are never normalized together.

No automatic length bonus is introduced.

## 09A.20 External validation

Muita and Damayo synthetic treatment datasets are external validation only.

Current Muita CHH work is exploratory because the true administered trigger sequence is not publicly available.

No external validation result may alter fitted Stage 09A coefficients post hoc.

---

# 09B — Layer 2: guide competence / strand-selection biophysics

## 09B.1 Inputs

Primary:

```text
asymmetry_ddg_4bp
guide_self_fold_mfe_kcal_mol
```

Sensitivity:

```text
asymmetry_ddg_5bp
```

## 09B.2 Favourable directions

```text
higher / more positive asymmetry_ddg_4bp
= more classically guide-favouring

higher / less-negative guide_self_fold_mfe
= weaker predicted self-folding
= favourable direction
```

## 09B.3 Separate length analyses

For each target, calculate separately:

```text
23-nt distributions
24-nt distributions
23-nt percentiles
24-nt percentiles
23-nt correlations
24-nt correlations
```

## 09B.4 Normalized components

Within `target × candidate_length`:

```text
layer2_asymmetry_percentile
layer2_self_fold_percentile
```

## 09B.5 Neutral reference

Because no Varroa efficacy data identify the correct within-layer weight:

```text
layer2_reference_score
=
0.5 × layer2_asymmetry_percentile
+
0.5 × layer2_self_fold_percentile
```

This is a **neutral reference**, not a learned biological weighting.

## 09B.6 Sensitivity family

For:

```text
alpha = 0.0, 0.1, ..., 1.0
```

calculate:

```text
L2(alpha)
=
alpha × asymmetry_percentile
+
(1-alpha) × self_fold_percentile
```

No alpha is biologically canonical at Stage 09.

## 09B.7 Redundancy

Within each target × length:

```text
Spearman(asymmetry_4bp, self_fold_MFE)
Spearman(asymmetry_4bp, asymmetry_5bp)
```

Also retain:

```text
layer2_component_difference
=
asymmetry_percentile
-
self_fold_percentile
```

No disagreement gate.

---

# 09C — Layer 3: target engagement / predicted target accessibility

## 09C.1 Inputs

Primary Stage 08:

```text
target_whole_p_unpaired
target_seed_g2_8_p_unpaired
```

Sensitivity:

```text
W100/L80
W200/L150
```

## 09C.2 Favourable direction

Higher predicted accessibility = more favourable.

## 09C.3 Separate length analyses

23 nt and 24 nt remain distinct.

This is especially important for whole-site simultaneous-unpaired probabilities because the interval length itself affects the probability scale.

## 09C.4 Normalized components

Within `target × candidate_length`:

```text
layer3_whole_accessibility_percentile
layer3_seed_accessibility_percentile
```

## 09C.5 Neutral reference

```text
layer3_reference_score
=
0.5 × whole_accessibility_percentile
+
0.5 × seed_accessibility_percentile
```

Neutral reference only.

## 09C.6 Sensitivity family

For:

```text
gamma = 0.0, 0.1, ..., 1.0
```

calculate:

```text
L3(gamma)
=
gamma × whole_accessibility_percentile
+
(1-gamma) × seed_accessibility_percentile
```

No gamma is biologically canonical.

## 09C.7 Robustness

Within each target × length report:

```text
whole vs seed rho
whole main vs W100/L80 rho
whole main vs W200/L150 rho
seed main vs W100/L80 rho
seed main vs W200/L150 rho
```

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

# 09.3 No final inter-layer weighting

Stage 09 must not define:

```text
overall_stage09_score
overall_candidate_score
efficacy_probability
final_candidate_rank
best_candidate
```

Stage 10 will address:

- inter-layer weight uncertainty;
- robust overall ranking;
- integration of 23/24 information at region level;
- long-dsRNA scoring;
- construct architecture;
- junction penalties;
- final design selection.

---

# 09.4 Canonical outputs

```text
results/09_feature_layers/
│
├── 09A_layer1_accumulation/
│   ├── layer1_model_coefficients.tsv
│   ├── layer1_model_preprocessing.tsv
│   ├── layer1_model_selection.tsv
│   ├── layer1_cv_by_group.tsv
│   ├── layer1_cv_summary_23nt.tsv
│   ├── layer1_cv_summary_24nt.tsv
│   ├── layer1_representation_diagnostic.tsv
│   ├── layer1_architecture_benchmarks.tsv
│   └── candidate_layer1.tsv
│
├── 09B_layer2_guide_competence/
│   ├── candidate_layer2.tsv
│   ├── layer2_weight_sensitivity_23nt.tsv
│   ├── layer2_weight_sensitivity_24nt.tsv
│   └── layer2_correlations.tsv
│
├── 09C_layer3_target_engagement/
│   ├── candidate_layer3.tsv
│   ├── layer3_weight_sensitivity_23nt.tsv
│   ├── layer3_weight_sensitivity_24nt.tsv
│   └── layer3_correlations.tsv
│
├── candidate_stage09_layers.tsv
├── stage09_feature_correlations.tsv
├── stage09_parameters.tsv
└── stage09_qc.tsv
```

The joined candidate table may contain both lengths, but every specified normalization/statistical analysis remains length-stratified.

---

# 10. Future website / app extensibility principle

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
- the current Stage 09A empirical accumulation model must not be silently described as validated at an untrained length.

If multiple user-requested lengths are analysed, each length remains its own normalization/analysis stratum unless a future validated specification explicitly says otherwise.

Recommended settings should guide users without unnecessarily restricting them.

---

# 11. Reproducibility parameters

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
Stage 09 model family
Stage 09 predictor encoding/scaling
Stage 09 sample-weighting rule
Stage 09 nuisance-intercept treatment
Stage 09 CV scheme
Stage 09 model-structure choice
Stage 09 hyperparameters
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

# 12. Required deterministic tests

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
- held-out-group leakage prevention;
- exact sample-aware weight regression;
- exact training-only feature scaling;
- exact fixed-effect within-group centering for accumulation model;
- unpenalized nuisance intercept handling;
- hyperparameter-grid and tie-break regression;
- 23/24 validation separate;
- 23/24 normalization separate;
- shared coefficients do not imply pooled analysis;
- representation remains diagnostic;
- Layer-2 favourable directions;
- Layer-3 favourable directions;
- neutral 0.5 references exact;
- alpha/gamma sensitivity arithmetic;
- Stage 06/08/09 row counts identical;
- no final overall rank.

---

# 13. Explicitly superseded / excluded analyses

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
- final Stage 10 inter-layer weighting;
- final Vd-CHIBIN region/construct selection;
- construct concatenation scoring;
- full user-facing Nectar Designer implementation.

---

# 14. Definition of current pipeline success

The canonical build is successful when:

1. a fresh clone can point to the frozen validated core;
2. Stage 00 validates without modifying it;
3. Stages 01–05 regenerate their canonical viral analyses from frozen inputs;
4. Stage 06 generically enumerates requested transcript candidate lengths;
5. Stage 07 reconstructs matched-background empirical guide-sequence associations;
6. Stage 08 preserves every candidate while computing raw biophysics;
7. Stage 09A learns a leakage-free empirical accumulation layer using the eight specified predictors;
8. 23 nt and 24 nt remain separate analyses even if Stage 09A shares coefficients;
9. Stage 09B and 09C generate length-stratified normalized evidence layers;
10. no Muita/Damayo external result influences training;
11. no arbitrary cross-layer efficacy weighting is introduced in Stage 09;
12. all important metrics are documented in `docs/METRIC_DICTIONARY.md`;
13. all major calculations have deterministic tests;
14. parameters/provenance are machine-readable;
15. no upstream frozen input is modified.

---

# 15. Methodological references

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
- Zou H, Hastie T. 2005. *Regularization and variable selection via the elastic net.* JRSS B 67:301–320. DOI: `10.1111/j.1467-9868.2005.00503.x`.

Evidence from non-Varroa systems is used as mechanistic/comparative support only. Their numerical coefficients are not transferred into the canonical Varroa empirical model.
