# Varroa vsiRNA Metric Dictionary

**Version:** 0.10  
**Scope:** Canonical viral pipeline through viral spatial/transitivity-consistency analysis

---

## 1. Purpose and metric classes

This document gives every important metric one fixed meaning.

Each metric is labelled as one of:

- **Standard/descriptive** — conventional mathematical/statistical quantity.
- **Published-method-derived** — output or definition taken directly from a published method such as stepRNA.
- **Project-specific** — designed for this Varroa analysis to answer a defined biological question.
- **Historical** — exact definition needed to reproduce v1.4.1.
- **Canonical** — preferred definition/aggregation for the new pipeline.

A project-specific metric can be scientifically useful; the label simply prevents us from presenting it as a universally established RNAi statistic.

---

## 2. Analysis units

### `sample`

One sequencing library/sample identifier.

For canonical cross-dataset inference, this is the top-level clustering unit because multiple viruses can be observed within the same library.

### `sample_virus_unit`

One virus analysed within one sample.

### `sample_virus_contig_unit`

One viral reference contig analysed within one sample-virus unit.

The historical v1.4.1 output happened to contain one analysed contig per included sample-virus unit, but code must not assume that this is universal.

---

## 3. Weighting modes

### `abundance_weighted`

**Class:** Standard/descriptive

A sequence contributes according to observed read abundance.

Question answered:

> What dominates the accumulated sequenced small-RNA population?

### `unique_sequence`

**Class:** Standard concept with project-specific analysis-unit definition

Each distinct RNA sequence contributes total weight 1 within the explicitly stated unit.

Question answered:

> Is the pattern represented across many distinct RNA sequence species rather than being driven by a few very abundant reads?

For Stage 01, unique-sequence identity is:

```text
sample × analysis_unit × length × strand × sequence
```

Equivalently, each distinct `sequence` contributes once within each:

```text
sample × analysis_unit × length × strand
```

unit. The same sequence can therefore contribute again in another sample, analysis unit, length, or mapped strand.

For Stage 05, where positional multimapping matters, identity is:

```text
virus × strand × length × sequence
```

and total sequence weight 1 is divided across all exact compatible loci for that sequence within the relevant sample-virus analysis.

---

# 4. Stage 01 length-landscape and 23/24 population metrics

## `length_count(L)`

**Class:** Standard/descriptive

Count/weight of eligible viral small RNAs of length `L` within one sample-virus unit and one weighting mode.

For abundance weighting:

```text
length_count(L)
    = sum(read-level count for eligible rows of length L)
```

For unique-sequence weighting:

```text
length_count(L)
    = number of distinct Stage 01 sequences of length L
```

after strand-specific deduplication.

Canonical Stage 01 evaluates `L = 15, ..., 35`.

Raw abundance counts are strongly influenced by viral load and sequencing depth and should not be interpreted as directly normalized cross-sample expression measurements.

## `length_fraction(L)`

**Class:** Standard proportion used as the primary comparable length-spectrum quantity

Within one sample-virus unit and weighting mode:

```text
length_fraction(L)
    = length_count(L)
      / Σ length_count(k),  k = 15,...,35
```

Range:

```text
0 to 1
```

If the 15–35-nt denominator is zero, report `NA`.

Interpretation:

> What fraction of the canonical retained viral small-RNA population belongs to length L?

This normalizes the length spectrum within a sample-virus unit; it does not normalize absolute viral small-RNA load between samples.

## `length_rank(L)`

**Class:** Standard descriptive ranking

Lengths are ordered by descending `length_fraction` within each sample-virus unit and weighting mode.

Use standard competition ranking with minimum rank for ties:

```text
1, 2, 2, 4, ...
```

Thus tied lengths receive the same best applicable rank and no arbitrary tie-breaking is introduced.

Interpretation:

```text
1 = most represented length in that unit
larger value = lower relative representation
```

## `top1_indicator(L)` and `top3_indicator(L)`

**Class:** Standard/descriptive

```text
top1_indicator(L) = 1 if length_rank(L) <= 1 else 0
top3_indicator(L) = 1 if length_rank(L) <= 3 else 0
```

Because tied lengths share rank, more than one length may be classified as top 1 or top 3 in a tied unit.

Across sample-virus units, the corresponding fractions are descriptive robustness summaries only; they are not formal hypothesis tests and do not replace sample-balanced summaries.

## `count_23`

**Class:** Standard/descriptive

Total eligible 23-nt count/weight within the stated sample-virus unit and weighting mode:

```text
count_23 = count_23_sense + count_23_antisense
```

## `count_24`

Equivalent quantity for 24-nt viral small RNAs:

```text
count_24 = count_24_sense + count_24_antisense
```

## `antisense_fraction_23`

**Class:** Standard/descriptive

```text
antisense_fraction_23
    = count_23_antisense
      / (count_23_sense + count_23_antisense)
```

If the denominator is zero, report `NA`.

## `antisense_fraction_24`

**Class:** Standard/descriptive

```text
antisense_fraction_24
    = count_24_antisense
      / (count_24_sense + count_24_antisense)
```

If the denominator is zero, report `NA`.

Interpretation for either metric:

```text
0.5  = equal sense and antisense contribution
>0.5 = antisense-biased
<0.5 = sense-biased
```

**Naming rule:** `antisense_fraction_24` must not be abbreviated to `F24_AS`. Stage 05 reserves `F24_AS`-type notation for the different concept of 24-nt composition within the antisense 23+24 population.

## `sense_fraction_23` and `sense_fraction_24`

```text
sense_fraction_23
    = count_23_sense / (count_23_sense + count_23_antisense)

sense_fraction_24
    = count_24_sense / (count_24_sense + count_24_antisense)
```

When defined:

```text
sense_fraction_23 + antisense_fraction_23 = 1
sense_fraction_24 + antisense_fraction_24 = 1
```

## `delta_antisense_fraction_24_minus_23`

**Class:** Standard difference in proportions used as a project-specific descriptive effect size

```text
delta_antisense_fraction_24_minus_23
    = antisense_fraction_24
      - antisense_fraction_23
```

Interpretation:

```text
>0 = 24 nt is more antisense-biased than 23 nt
 0 = equal antisense fraction
<0 = 24 nt is less antisense-biased than 23 nt
```

Range when both fractions are defined:

```text
-1 to +1
```

A value of `+0.30` means the 24-nt population is 30 percentage points more antisense-biased than the 23-nt population in that analysis unit.

This is a population contrast, not evidence that an individual 24-mer is secondary or an individual 23-mer is primary.

## `length23_fraction_among_23_24`

**Class:** Standard/descriptive composition

```text
length23_fraction_among_23_24
    = count_23 / (count_23 + count_24)
```

## `length24_fraction_among_23_24`

```text
length24_fraction_among_23_24
    = count_24 / (count_23 + count_24)
```

When defined:

```text
length23_fraction_among_23_24
+ length24_fraction_among_23_24
= 1
```

If `count_23 + count_24 = 0`, both are `NA`.

These are not the Stage 05 `F24_AS` quantities.

## `sample_balanced_median(metric)`

**Class:** Project-specific canonical aggregation rule

For a continuous Stage 01 pair-level metric:

1. calculate the metric separately for each eligible sample-virus unit;
2. within each sample, take the median across eligible virus units;
3. across samples, take the median of those sample-level values.

Primary uncertainty is obtained by resampling biological samples with replacement and recomputing the summary while retaining each selected sample's virus observations together.

This prevents samples containing more eligible viruses from automatically receiving more independent weight.

For length-spectrum metrics, this rule is applied separately for every length and weighting mode.

## `sample_clustered_CI95`

**Class:** Standard cluster/bootstrap uncertainty approach applied canonically

A 95% bootstrap confidence interval obtained by resampling top-level biological samples with replacement and recomputing the pre-specified sample-balanced statistic.

The implementation must record:

- number of requested bootstrap replicates;
- number of valid replicates;
- random seed;
- interval construction method.

Stage 01 emphasizes effect sizes and confidence intervals; it does not require P-values simply to describe prominence of lengths or strand bias.

---

# 5. Terminal nucleotide coordinates

Terminal positions are defined relative to the **physical RNA sequence in its own 5′→3′ orientation**:

```text
5′ N1 N2 ................ N(n-1) Nn 3′
   ↑  ↑                    ↑     ↑
 5p1 5p2                  3p2   3p1
```

Canonical Stage 02 TSVs use the sequencing alphabet:

```text
A, C, G, T
```

where `T` corresponds biologically to uridine (`U`) in RNA.

Observed antisense read sequences are already evaluated in their sequenced 5′→3′ orientation and are **not** reverse-complemented again. Expected antisense sequences are created by reverse-complementing reference-orientation viral windows and then evaluating the resulting antisense sequence 5′→3′.

---

# 6. Terminal nucleotide enrichment

## `observed_terminal_weight(b,p)`

**Class:** Standard/descriptive

Within one:

```text
sample × analysis_unit × length × strand × weighting_mode
```

this is the total observed weight of eligible RNAs carrying nucleotide `b ∈ {A,C,G,T}` at terminal position `p ∈ {5p1,5p2,3p2,3p1}`.

Abundance mode uses the numeric read-level `count` field.

Unique-sequence mode gives each distinct:

```text
sample × analysis_unit × length × strand × sequence
```

total observed weight 1.

## `observed_total_weight`

Total eligible observed weight in the same population.

## `observed_fraction(b,p)`

**Class:** Standard/descriptive

```text
observed_fraction(b,p)
    = observed_terminal_weight(b,p)
      / observed_total_weight
```

Range:

```text
0 to 1
```

If `observed_total_weight = 0`, report `NA`.

For every valid population and terminal position:

```text
Σ_b observed_fraction(b,p) = 1
```

within numerical tolerance.

Question answered:

> Among the observed viral small-RNA population under the stated weighting mode, how common is nucleotide b at terminal position p?

---

## `valid_background_window_count(L)`

**Class:** Project-specific matched-background accounting quantity

For target length `L`, count all windows of length `L` that lie entirely within one FASTA record of the sample-specific depth-masked viral consensus and contain only:

```text
A, C, G, T
```

Windows containing `N` or another non-ACGT base are excluded. Windows never cross FASTA-record boundaries.

Every valid genomic start position contributes one background opportunity; repeated identical sequence strings at different positions remain separate opportunities.

---

## `expected_fraction_sense(b,p)`

**Class:** Project-specific matched-background quantity

```text
expected_fraction_sense(b,p)
    = number of fully supported reference-orientation windows
      carrying nucleotide b at position p
      / total fully supported windows of the same length
```

If no fully supported window exists, report `NA`.

## `expected_fraction_antisense(b,p)`

For each fully supported reference window, first calculate its reverse complement and then read terminal position `p` in that antisense sequence's own 5′→3′ orientation:

```text
expected_fraction_antisense(b,p)
    = number of fully supported reverse-complement windows
      carrying nucleotide b at position p
      / total fully supported windows of the same length
```

This is not obtained by simply relabelling reference ends; the reverse-complement transformation must be correct.

For every valid expected population and position:

```text
Σ_b expected_fraction(b,p) = 1
```

within numerical tolerance.

---

## `observed_strand_weight_sense` and `observed_strand_weight_antisense`

**Class:** Standard proportions used for matched combined expectation

For one sample-virus × length × weighting mode:

```text
wS
    = observed sense weight
      / (observed sense weight + observed antisense weight)

wAS
    = observed antisense weight
      / (observed sense weight + observed antisense weight)
```

When defined:

```text
wS + wAS = 1
```

The weights are calculated separately for abundance and unique-sequence modes because their observed strand mixtures may differ.

---

## `expected_fraction_combined(b,p)`

**Class:** Project-specific strand-matched background quantity

```text
expected_fraction_combined(b,p)
    = wS  × expected_fraction_sense(b,p)
      + wAS × expected_fraction_antisense(b,p)
```

The combined expected background is therefore matched to the observed strand mixture and is **not** forced to 50:50.

If the combined observed denominator is zero, report `NA`.

---

## `expected_fraction(b,p)`

Generic field name for the matched expected frequency under the stated strand scope:

```text
sense     → expected_fraction_sense
antisense → expected_fraction_antisense
combined  → expected_fraction_combined
```

The background is a positional viral-sequence opportunity background. The same valid genomic windows are used for abundance and unique-sequence observed modes; the weighting modes differ on the observed side and, for combined scope, in the observed strand-mixture weights.

---

## `enrichment_ratio(b,p)`

**Class:** Project-specific empirical effect size

Within one sample-virus unit:

```text
enrichment_ratio(b,p)
    = observed_fraction(b,p)
      / expected_fraction(b,p)
```

Interpretation:

```text
1   = observed as often as matched sequence availability predicts
>1  = enriched
<1  = depleted
0   = not observed although expected frequency is positive
```

If `expected_fraction = 0`, report `NA`.

If the observed population denominator is zero, report `NA`.

No arbitrary pseudocount is added.

The metric measures an empirical terminal sequence association; it does not identify the molecular mechanism that produced the association.

---

## `pair_median_enrichment_ratio`

**Class:** Project-specific pair-balanced / historical-regression summary

For a fixed:

```text
length × strand_scope × weighting_mode × terminal_position × nucleotide
```

```text
pair_median_enrichment_ratio
    = median of finite enrichment_ratio values
      across eligible sample-virus units
```

This corresponds most closely to the project's historical design-facing `median_enrichment_ratio` logic.

It is not the canonical primary cross-dataset inference because several pair observations can come from the same biological sample.

### Historical name `median_enrichment_ratio`

If a historical-compatible output field is named exactly:

```text
median_enrichment_ratio
```

it is treated as an alias of the pair-balanced historical/regression quantity, **not** as an alias of the canonical sample-balanced result.

This prevents old design references from being silently redefined.

---

## `sample_enrichment_median`

**Class:** Canonical intermediate aggregation

Within one sample, for a fixed terminal feature:

```text
sample_enrichment_median(sample)
    = median of finite enrichment_ratio values
      across that sample's eligible viruses
```

The median is taken over already abundance-weighted or unique-sequence-weighted **pair-level enrichment ratios**. Abundance information has therefore already entered the pair-level estimate before the sample balancing occurs.

---

## `sample_balanced_median_enrichment_ratio`

**Class:** Project-specific canonical primary summary

```text
sample_balanced_median_enrichment_ratio
    = median of sample_enrichment_median(sample)
      across contributing biological samples
```

Equivalent procedure:

1. calculate matched enrichment per sample-virus;
2. median across viruses within each sample;
3. median across samples.

This is explicitly **not**:

```text
median(observed_fraction) / median(expected_fraction)
```

and is not a pooled-read estimate.

Question answered:

> What terminal enrichment is typical across biological samples, while preventing deeper libraries or samples containing more viruses from automatically dominating the final estimate?

Primary uncertainty uses `sample_clustered_CI95`.

---

## `pooled_abundance_observed_fraction(b,p)`

**Class:** Project-specific secondary descriptive molecular-pool quantity

Defined only for abundance weighting.

For eligible sample-virus units `u` with defined matched backgrounds, let:

```text
N_u = total observed abundance in the stated length/strand scope
O_u = observed abundance carrying nucleotide b at position p
```

Then:

```text
pooled_abundance_observed_fraction(b,p)
    = Σ_u O_u / Σ_u N_u
```

This allows high-abundance/high-depth infections to contribute proportionally more molecular weight.

---

## `pooled_abundance_expected_fraction(b,p)`

**Class:** Project-specific abundance-matched expected quantity

Let `E_u = expected_fraction_u(b,p)` for the same unit and feature.

```text
pooled_abundance_expected_fraction(b,p)
    = Σ_u (N_u × E_u) / Σ_u N_u
```

For combined scope, `E_u` is the abundance-mode strand-weighted combined expected fraction.

This preserves each sample-virus unit's matched viral sequence composition while matching the overall expected background to the molecular abundance contributing to the pooled observed population.

---

## `pooled_abundance_enrichment_ratio(b,p)`

**Class:** Project-specific secondary descriptive effect size

```text
pooled_abundance_enrichment_ratio(b,p)
    = pooled_abundance_observed_fraction(b,p)
      / pooled_abundance_expected_fraction(b,p)
```

If the pooled observed denominator is zero or the pooled expected fraction is zero, report `NA`.

Question answered:

> Across the total accumulated eligible molecular pool, what terminal enrichment is observed when molecular abundance is allowed to weight infections unequally?

This is intentionally different from `sample_balanced_median_enrichment_ratio`.

Interpretation rule:

- useful secondary description of the accumulated molecular pool;
- not an estimate of the typical biological sample;
- no read-level P-value;
- no automatic inferential confidence interval;
- never treated as millions of independent biological replicates;
- not automatically selected as the later vdCHIBIN design reference.

The output must retain total contributing observed abundance, number of contributing sample-virus units, and number of contributing biological samples.

---

## `spearman_rho_23_24`

**Class:** Standard statistical metric applied to canonical Stage 02 features

Spearman rank correlation between matched 23-nt and 24-nt terminal enrichment landscapes.

Canonical input vector for each length consists of the same:

```text
4 terminal positions × 4 nucleotides = 16 features
```

using `sample_balanced_median_enrichment_ratio`.

At minimum calculate separately for:

```text
combined strand scope
antisense strand scope
```

and separately for abundance and unique-sequence weighting.

Interpretation:

```text
+1 = identical rank order of terminal enrichment features
 0 = no monotonic association
-1 = opposite rank order
```

This measures similarity of empirical enrichment patterns, not shared enzymatic origin.

---

## Relationship between Stage 02 abundance and median aggregation

The following quantities answer different levels of the biological question and must not be conflated:

```text
abundance weighting within sample-virus
    → which terminal features dominate accumulated molecules in that infection

sample-balanced median across biological samples
    → whether that abundance-weighted preference is reproducible across infections

pooled-abundance secondary summary
    → which feature dominates the total sequenced molecular pool when deep/high-load infections receive more molecular weight
```

Thus using a median at the cross-sample level does **not** remove abundance information from the abundance-mode analysis; abundance has already determined each sample-virus enrichment ratio.

---


# 7. Stage 03 official stepRNA geometry metrics

These definitions preserve the published stepRNA conventions while distinguishing official stepRNA counts from project-specific Varroa abundance weighting.

## `focal_reference_sequence`

**Class:** Project-specific input unit built for a published method

One distinct File-A sequence within:

```text
sample × analysis_unit × focal_length × focal_strand
```

Stage 03 uses one FASTA record per distinct focal sequence.

The focal sequence's observed molecular abundance is retained separately as `focal_abundance`.

## `focal_abundance`

**Class:** Standard/descriptive

For one Stage 03 focal reference sequence:

```text
focal_abundance
    = sum(read-level count for all canonical rows
          carrying that focal sequence in the same
          sample × analysis_unit × focal_length × focal_strand)
```

This field is used for project abundance-weighted reference-support summaries after stepRNA geometry reconstruction.

It is not inferred from the number of stepRNA alignments.

## `passenger_candidate_sequence`

**Class:** Project-specific input unit built for a published method

One distinct File-B sequence within the same sample-virus unit as File A, from the opposite mapped viral strand, with canonical passenger length:

```text
15–30 nt inclusive
```

Passenger FASTA sequences remain in observed physical 5′→3′ orientation.

---

## `steprna_5p_distance`

**Class:** Published-method-derived

Signed 5′ distance relative to the File-A reference RNA:

```text
negative = File-A reference overhang
0        = blunt end
positive = File-A reference underhang
```

No sign transformation is applied in canonical storage.

## `steprna_3p_distance`

**Class:** Published-method-derived

Signed 3′ distance relative to the File-A reference RNA using the same convention:

```text
negative = File-A reference overhang
0        = blunt end
positive = File-A reference underhang
```

The 5′ and 3′ values must refer to the same reconstructed focal/passenger duplex when used in a joint-geometry metric.

---

## `steprna_joint_geometry`

**Class:** Canonical pairing of published-method-derived quantities

For one recovered focal/passenger duplex:

```text
steprna_joint_geometry
    = (steprna_5p_distance, steprna_3p_distance)
```

Both distances must come from the **same official classified alignment**.

Examples:

```text
(0,0)     fully blunt at both analysed ends
(+2,-2)   pre-specified Varroa 2-nt joint geometry
(0,-1)    blunt/flush at the 5′ marginal end, -1 at the 3′ end
(-1,0)    -1 at the 5′ end, blunt/flush at the 3′ marginal end
```

This metric exists specifically to prevent marginal end-distance peaks from being misinterpreted as complete duplex geometry.

## `joint_geometry_duplex_count(d5,d3)`

**Class:** Canonical descriptive count

Number of recovered focal/passenger duplex alignments in one Stage 03 run with exactly:

```text
steprna_5p_distance = d5
steprna_3p_distance = d3
```

The counts across all observed `(d5,d3)` combinations must sum exactly to `n_recovered_duplexes`.

## `joint_duplex_fraction(d5,d3)`

**Class:** Canonical descriptive metric

```text
joint_duplex_fraction(d5,d3)
    = joint_geometry_duplex_count(d5,d3)
      / n_recovered_duplexes
```

For every non-empty run:

```text
Σ joint_duplex_fraction(d5,d3) = 1
```

within numerical precision.

If no duplex is recovered, report `NA`.

## `joint_00_duplex_fraction`

**Class:** Canonical named special case

```text
joint_00_duplex_fraction
    = joint_duplex_fraction(0,0)
```

Interpretation:

> Fraction of reconstructed duplex relationships that are simultaneously distance `0` at both analysed ends.

This must **not** be inferred from the marginal 5′ and 3′ distance-0 frequencies.

---

## `steprna_official_duplex_count`

**Class:** Published-method-derived/native software output

Distance-specific count from the official stepRNA overhang output in which a File-A reference may be represented multiple times through recovered duplex relationships.

This quantity must be preserved using the semantics of the installed official stepRNA version.

It must **not** automatically be renamed `abundance_weighted`, because the project-wide abundance definition is based on upstream observed focal read counts, not on the number of recovered passenger alignments.

## `steprna_official_unique_reference_count`

**Class:** Published-method-derived/native software output

Distance-specific count from the official stepRNA unique-reference overhang output.

Use the official software result directly; do not silently reimplement a different uniqueness rule under the same name.

---

## `passenger_length`

**Class:** Published-method-derived

Length of a complementary File-B passenger sequence recovered by official stepRNA.

Canonical primary File B permits passenger lengths:

```text
15–30 nt
```

## `passenger_count_per_reference`

**Class:** Published-method-derived/native software output

Number of recovered passenger relationships associated with a File-A focal reference as represented by official stepRNA passenger-number output.

This is not equivalent to passenger molecular abundance.

---

## `passenger_recovery_fraction_unique`

**Class:** Canonical descriptive metric based on stepRNA output

Within one sample-virus × focal-length × focal-strand run:

```text
passenger_recovery_fraction_unique
    = number of distinct focal_reference_sequence values
      with at least one recovered passenger
      / number of distinct eligible focal_reference_sequence values
```

Range:

```text
0 to 1
```

Interpretation:

> What fraction of focal sequence diversity has at least one recoverable complementary passenger?

A low value does not by itself demonstrate absence of Dicer processing.

## `passenger_recovery_fraction_abundance`

**Class:** Canonical descriptive metric

```text
passenger_recovery_fraction_abundance
    = Σ focal_abundance for focal references with ≥1 recovered passenger
      / Σ focal_abundance for all eligible focal references
```

Range:

```text
0 to 1
```

Interpretation:

> What fraction of accumulated focal small-RNA abundance belongs to focal sequences for which at least one passenger is recoverable?

This uses upstream focal abundance and does not expand FASTA records according to read count.

---

## `steprna_log_ratio`

**Class:** Published-method-derived

The official stepRNA distance-enrichment/log-ratio quantity produced by the pinned software version.

The published method describes this as the logarithm of a ratio comparing a distance-specific overhang count with the mean count across end distances.

The canonical pipeline archives and parses the official value; it does not silently substitute another odds-ratio or enrichment formula under the same name.

## `steprna_wald_z`

**Class:** Published-method-derived

Official stepRNA Wald Z-score associated with distance enrichment.

Use the value produced by the pinned official stepRNA implementation.

Stage 03 records this statistic; Stage 04 determines how population-level evidence should be summarized across biological samples.

---

# 8. Pre-specified Varroa joint Dicer-like geometry

The term “canonical Dicer geometry” is avoided as a universal biological label because Dicer/Dicer-like end geometry can vary across organisms, pathways, and substrates.

## `varroa_2nt_joint_geometry`

**Class:** Project-specific, pre-specified pathway feature

Defined under the official stepRNA sign convention as:

```text
steprna_5p_distance = +2
steprna_3p_distance = -2
```

Equivalent interpretation:

```text
File-A 5′ underhang = 2 nt
File-A 3′ overhang  = 2 nt
```

Equivalent historical label:

```text
5p_underhang_2__3p_overhang_2
```

The two distances must come from the **same reconstructed focal/passenger duplex**.

This feature is pre-specified before inspecting the canonical Stage 03 results and is not a universal definition of Dicer cleavage.

## `n_joint_geometry_duplexes`

**Class:** Canonical descriptive count

Number of recovered focal/passenger duplex alignments in one Stage 03 run for which:

```text
5p = +2
and
3p = -2
```

## `varroa_2nt_joint_duplex_fraction`

**Class:** Canonical descriptive metric

```text
varroa_2nt_joint_duplex_fraction
    = n_joint_geometry_duplexes
      / n_recovered_duplexes
```

where `n_recovered_duplexes` is the number of reconstructed focal/passenger duplex relationships in the same run.

This is a **duplex-level** quantity. A focal reference with many recovered passengers can contribute multiple duplexes.

If no duplex is recovered, report `NA`.

## `n_focal_references_supporting_joint_geometry`

**Class:** Canonical descriptive count

Number of distinct File-A focal references with at least one recovered passenger duplex satisfying:

```text
5p = +2
3p = -2
```

A focal reference may also support other geometries through other recovered passengers.

## `varroa_2nt_reference_fraction_all`

**Class:** Canonical descriptive metric

```text
varroa_2nt_reference_fraction_all
    = n_focal_references_supporting_joint_geometry
      / n_focal_references
```

This combines passenger recoverability and joint-geometry support.

It therefore must not be interpreted alone as the geometry distribution among recovered duplexes.

## `varroa_2nt_reference_fraction_recovered`

**Class:** Canonical descriptive metric

```text
varroa_2nt_reference_fraction_recovered
    = n_focal_references_supporting_joint_geometry
      / n_recovered_focal_references
```

Interpretation:

> Among distinct focal sequences for which at least one passenger is recoverable, what fraction support the pre-specified Varroa joint geometry at least once?

If no focal reference has a recovered passenger, report `NA`.

## `varroa_2nt_reference_fraction_abundance_all`

**Class:** Canonical project-specific abundance-support metric

```text
varroa_2nt_reference_fraction_abundance_all
    = Σ focal_abundance for focal references supporting joint geometry
      / Σ focal_abundance for all eligible focal references
```

Interpretation:

> What fraction of accumulated focal abundance is represented by focal sequences that support the pre-specified joint geometry?

This is a reference-support metric, not a mutually exclusive geometry distribution.

## `varroa_2nt_reference_fraction_abundance_recovered`

**Class:** Canonical project-specific abundance-support metric

```text
varroa_2nt_reference_fraction_abundance_recovered
    = Σ focal_abundance for focal references supporting joint geometry
      / Σ focal_abundance for focal references with ≥1 recovered passenger
```

If the denominator is zero, report `NA`.

Interpretation:

> Conditional on focal sequences having a recoverable passenger, what fraction of accumulated focal abundance belongs to sequences supporting the pre-specified joint geometry?

---

## Full-spectrum versus joint-geometry rule

The **full signed marginal 5′ and 3′ distance spectra and the full same-duplex `(d5,d3)` spectrum are complementary canonical Stage 03 outputs**.

A strong marginal peak at distance `0` does not imply that `(0,0)` is common in the same duplex. The joint spectrum must be consulted before using words such as “fully blunt duplex.”

`varroa_2nt_joint_geometry` is a pre-specified secondary feature of interest inside the full same-duplex landscape.

A strong +2/−2 result must not cause other reproducibly enriched distances to be discarded, and a different dominant distance must not be relabelled +2/−2 after the fact.

---

## Stage 03 abundance-weighting rule

Three concepts must remain separate:

```text
official stepRNA duplex/alignment count
official stepRNA unique-reference count
canonical focal-abundance-weighted reference support
```

The canonical project does **not** treat the number of passenger alignments as a substitute for upstream small-RNA abundance.

Abundance weighting uses `focal_abundance` derived from the canonical read-level `count` field.

---


# 9. Stage 04 sample-aware duplex-geometry metrics

## `sample_steprna_log_ratio_median`

**Class:** Canonical aggregation of a published-method-derived metric

For one fixed:

```text
sample
× focal_length
× focal_strand
× end
× signed_distance
× official_view
```

```text
sample_steprna_log_ratio_median
    = median of finite official stepRNA log-ratios
      across eligible sample-virus units within that sample
```

`official_view` is one of:

```text
duplex
unique_reference
```

Purpose: prevent a sequencing library containing several eligible viruses from automatically contributing several independent biological votes.

## `sample_balanced_steprna_log_ratio`

**Class:** Canonical

```text
sample_balanced_steprna_log_ratio
    = median across contributing samples of
      sample_steprna_log_ratio_median
```

Primary uncertainty is the canonical sample-clustered percentile bootstrap 95% CI.

This is the preferred cross-dataset effect summary for the full official stepRNA distance spectrum.

It is not an official stepRNA statistic itself; it is a canonical aggregation of official stepRNA per-run log-ratios.

## `sample_balanced_steprna_wald_z_descriptive`

**Class:** Canonical descriptive aggregation of a published-method-derived statistic

If reported, this is the median-across-viruses then median-across-samples summary of official per-run Wald Z values.

It is **descriptive only**.

Do not interpret it as a newly derived population-level Wald test, and do not derive a population P-value from it.

---

## `sample_balanced_joint_duplex_fraction(d5,d3)`

**Class:** Canonical same-duplex geometry aggregation

For each fixed focal length, focal strand, and same-duplex geometry `(d5,d3)`:

```text
pair-level joint_duplex_fraction(d5,d3)
→ median across viruses within sample
→ median across samples
```

with sample-clustered percentile bootstrap 95% CI.

This is the primary cross-dataset estimator for the prevalence of a specific **complete same-duplex geometry**.

## `sample_balanced_joint_00_duplex_fraction`

**Class:** Canonical named special case

```text
sample_balanced_joint_00_duplex_fraction
    = sample_balanced_joint_duplex_fraction(0,0)
```

This measures fully blunt reconstructed duplex relationships and is distinct from the sample-balanced marginal distance-0 log-ratio.

## `joint_geometry_mode_by_pair`

**Class:** Canonical descriptive

The `(d5,d3)` combination with the largest `joint_geometry_duplex_count` within one sample-virus/focal-class run.

A mode is the single most common category; it does not imply a majority.

Ties must be represented deterministically and transparently rather than silently broken as biological evidence.

## `joint_00_is_mode`

**Class:** Canonical descriptive indicator

```text
1 if (0,0) is a joint-geometry mode in the run
0 otherwise
```

Used to count how often fully blunt geometry is the most common complete geometry across sample-virus runs.

## `joint_00_majority`

**Class:** Canonical descriptive indicator

```text
1 if joint_00_duplex_fraction > 0.5
0 otherwise
```

This directly answers whether most reconstructed duplexes in a run are fully blunt.

In the validated current dataset this is `0/54` for all four focal classes.

## `pooled_joint_duplex_fraction(d5,d3)`

**Class:** Secondary descriptive

```text
Σ joint_geometry_duplex_count(d5,d3)
/
Σ n_recovered_duplexes
```

across the selected sample-virus units.

This is useful for intuitive reporting but is not the primary biological estimator because high-depth runs contribute more duplex relationships.

---

## `paired_delta_24_minus_23(M)`

**Class:** Canonical paired effect size

For the same:

```text
sample
analysis_unit
focal_strand
```

and the same metric definition `M`:

```text
paired_delta_24_minus_23(M)
    = M_24 - M_23
```

Examples include:

```text
paired_delta_24_minus_23(varroa_2nt_reference_fraction_recovered)
paired_delta_24_minus_23(sample-virus official 3p -2 log-ratio)
```

Canonical dataset aggregation:

```text
pair-level 24-minus-23 difference
→ median across viruses within sample
→ median across samples
```

with sample-clustered bootstrap 95% CI.

Positive values mean the specified metric is larger for 24 nt; negative values mean it is larger for 23 nt.

Sense and antisense focal-strand comparisons are kept separate.

---

## `sample_balanced_passenger_recovery_fraction_unique`

**Class:** Canonical descriptive

Sample-balanced summary of Stage 03 `passenger_recovery_fraction_unique`.

Passenger recovery measures observability of complementary partners, not Dicer activity.

## `sample_balanced_passenger_recovery_fraction_abundance`

**Class:** Canonical descriptive

Sample-balanced summary of Stage 03 `passenger_recovery_fraction_abundance`.

---

## `sample_balanced_varroa_2nt_joint_duplex_fraction`

**Class:** Canonical descriptive

Sample-balanced summary of the Stage 03 duplex-level `varroa_2nt_joint_duplex_fraction`.

This remains distinct from focal-reference support.

## `sample_balanced_varroa_2nt_reference_fraction_recovered`

**Class:** Canonical project-specific geometry-support metric

Sample-balanced summary of:

```text
n_focal_references_supporting_joint_geometry
/ n_recovered_focal_references
```

This is the preferred unique-reference joint-support view because it conditions on a passenger having been recoverable.

## `sample_balanced_varroa_2nt_reference_fraction_abundance_recovered`

**Class:** Canonical project-specific geometry-support metric

Sample-balanced summary of the abundance-weighted recovered-reference joint-support fraction.

This uses upstream `focal_abundance`, not passenger-alignment multiplicity.

---

## `pair_balanced_geometry_summary`

**Class:** Secondary descriptive/sensitivity aggregation

Median across sample-virus units for a fixed geometry metric.

Useful for regression with historical pair-balanced analyses but not the primary cross-dataset estimator.

## `virus_balanced_geometry_summary`

**Class:** Secondary sensitivity aggregation

For a fixed geometry metric:

```text
median across samples within biological_virus
→ median across biological viruses
```

This asks whether a result is broadly shared across virus identities rather than dominated by a repeatedly observed virus.

---

## 9A. Historical custom `Δ_Dicer`

## `historical_delta_dicer`

Written historically as `Δ_Dicer`.

**Class:** Historical regression-only metric

The earlier Varroa pipeline defined a project-specific statistic contrasting support at a pre-specified Dicer-like overhang feature against support at other tested distances and assessed it using a custom permutation framework.

The canonical v0.8 project does **not** assign a new formula to this historical name from memory or from the new stepRNA output.

It may be reproduced only from the exact archived v1.4.0 implementation/configuration, preserving:

- target geometry/distance definition;
- comparison-distance set;
- support definition;
- weighting mode;
- aggregation;
- permutation/null construction;
- random seed and replicate count when recorded.

If the exact definition is unavailable:

```text
historical_delta_dicer = NA
historical_delta_dicer_status = exact_definition_unavailable
```

Do not invent `D0` or substitute the official stepRNA log-ratio.

This metric is not a candidate-window score and is not part of the primary canonical Stage 04 inference.

---

# 10. Stage 04 geometry-conditioned terminal sequence metrics

The canonical names use **geometry** rather than **Dicer** because the subset is defined by reconstructed stepRNA geometry, not direct biochemical assignment to a nuclease.

## `joint_observed_fraction(f)`

**Class:** Canonical descriptive

For terminal feature `f` within one sample-virus × focal-length × focal-strand × weighting-mode unit:

```text
joint_observed_fraction(f)
    = weighted frequency of f among focal references
      supporting the pre-specified (+2,-2) joint geometry
```

Unique-sequence mode gives each distinct focal reference sequence weight 1.

Abundance mode weights each focal reference by upstream `focal_abundance`.

Do not weight by passenger multiplicity.

If no joint-supporting focal reference exists, report `NA`.

## `recovered_observed_fraction(f)`

**Class:** Canonical descriptive/control quantity

```text
recovered_observed_fraction(f)
    = weighted frequency of f among all focal references
      with at least one recovered passenger
```

This provides the passenger-recovery-conditioned comparison population.

## `E_joint_absolute(f)`

**Class:** Canonical exploratory/candidate-development metric

```text
E_joint_absolute(f)
    = joint_observed_fraction(f)
      / Stage02_expected_fraction(f)
```

The expected fraction must be the exact matching Stage 02 pair-level viral background for the same:

```text
sample × analysis_unit × focal_length × focal_strand
× terminal_position × nucleotide
```

Interpretation:

> Absolute enrichment of feature `f` among focal RNAs supporting the pre-specified joint geometry relative to matched viral sequence opportunity.

If the subset is empty or the expected fraction is zero/undefined, report `NA`.

No pseudocount.

## `E_recovered_absolute(f)`

**Class:** Canonical control metric

```text
E_recovered_absolute(f)
    = recovered_observed_fraction(f)
      / Stage02_expected_fraction(f)
```

This measures terminal enrichment among all passenger-recovered focal RNAs.

## `E_all(f)`

**Class:** Canonical upstream metric reference

The exact matching Stage 02 pair-level `enrichment_ratio` for all observed focal RNAs of the same sample-virus × length × strand × weighting mode × terminal feature.

Do not substitute the Stage 02 across-dataset median when calculating a pair-level contrast.

---

## `joint_vs_all_log2_contrast(f)`

**Class:** Canonical exploratory/candidate-development metric

```text
joint_vs_all_log2_contrast(f)
    = log2(E_joint_absolute(f) / E_all(f))
```

Interpretation:

```text
0  = joint-supporting subset resembles the overall observed enrichment
>0 = relatively more enriched in the joint-supporting subset
<0 = relatively less enriched in the joint-supporting subset
```

If either required enrichment is non-positive or undefined, report `NA`.

No pseudocount.

## `joint_vs_recovered_log2_contrast(f)`

**Class:** Canonical preferred geometry-specific sequence contrast

```text
joint_vs_recovered_log2_contrast(f)
    = log2(E_joint_absolute(f) / E_recovered_absolute(f))
```

When the same background is defined in numerator and denominator:

```text
= log2(joint_observed_fraction(f)
       / recovered_observed_fraction(f))
```

Interpretation:

> Conditional on a complementary passenger being recoverable, is feature `f` relatively more or less represented among focal sequences supporting the pre-specified (+2,-2) geometry?

This reduces confounding by passenger recoverability and is preferred over `joint_vs_all_log2_contrast` when asking whether the terminal feature is specifically associated with the reconstructed joint geometry.

It does not eliminate every possible sequencing/selection bias.

If a required quantity is zero or undefined, report `NA`; no pseudocount.

---

## `sample_balanced_E_joint_absolute(f)`

**Class:** Canonical aggregation

For each fixed terminal feature and analysis stratum:

```text
pair-level E_joint_absolute
→ median across viruses within sample
→ median across samples
```

with sample-clustered bootstrap 95% CI.

## `sample_balanced_joint_vs_all_log2_contrast(f)`

**Class:** Canonical aggregation

Sample-balanced median of pair-level `joint_vs_all_log2_contrast(f)` with sample-clustered bootstrap 95% CI.

## `sample_balanced_joint_vs_recovered_log2_contrast(f)`

**Class:** Canonical aggregation

Sample-balanced median of pair-level `joint_vs_recovered_log2_contrast(f)` with sample-clustered bootstrap 95% CI.

This is the primary terminal-sequence effect estimate for assessing geometry-specific information beyond passenger recoverability.

---

## `rho_joint_vs_general`

**Class:** Standard Spearman correlation applied descriptively to canonical metrics

Across the 16 matched terminal features for one focal length × focal strand × weighting mode:

```text
rho_joint_vs_general
    = Spearman rho(
        sample_balanced_E_joint_absolute,
        matching Stage 02 sample-balanced enrichment
      )
```

Purpose: quantify redundancy between the absolute joint-geometry enrichment landscape and the general Stage 02 terminal-enrichment landscape.

A high positive rho suggests the two landscapes have similar rank structure; it does not prove identical biological mechanism.

No formal correlation P-value is required for the canonical use.

## `rho_joint_contrast_abundance_vs_unique`

**Class:** Standard descriptive correlation

Across the same 16 terminal features:

```text
rho_joint_contrast_abundance_vs_unique
    = Spearman rho(
        sample_balanced_joint_vs_recovered_log2_contrast_abundance,
        sample_balanced_joint_vs_recovered_log2_contrast_unique
      )
```

Purpose: assess whether the geometry-specific terminal pattern is concordant between accumulated molecules and distinct sequence diversity.

---

## Deprecated pre-v0.7 names

The following names are deprecated before Stage 04 implementation because they imply stronger mechanistic assignment than the data justify:

```text
E_Dicer_absolute
dicer_specific_log2_contrast
dicer_general_correlation
```

They must not be used in new canonical output schemas.

Historical documents containing those names remain valid records of earlier planning, but the v0.7 canonical replacements are:

```text
E_joint_absolute
joint_vs_all_log2_contrast / joint_vs_recovered_log2_contrast
rho_joint_vs_general
```

---


# 11. Stage-05 coordinate and track metrics

## `alignment_midpoint_nt`

**Class:** Historical/project-specific coordinate convention

```text
alignment_midpoint_nt
    = alignment_start_0based + (read_length - 1) / 2
```

This is the coordinate used by v1.4.1 to place a 23/24-nt read on the viral genome.

## `bin_index`

With historical/default `bin_size_nt = 10`:

```text
bin_index = floor(alignment_midpoint_nt / 10)
```

## `abundance_locus_weight`

**Class:** Project-specific multimapper handling

For one QNAME with read abundance `a` and `k` unique exact compatible loci:

```text
abundance_locus_weight = a / k
```

Total weight over its loci is `a`.

## `unique_sequence_locus_weight`

**Class:** Project-specific corrected v1.4.1 definition

For one distinct sequence identity:

```text
virus × strand × length × sequence
```

with `k` unique exact compatible loci:

```text
unique_sequence_locus_weight = 1 / k
```

Total weight over all compatible loci is exactly 1.

---

# 12. 23-nt anchor metrics

## `balanced23_anchor_score`

Previously `balanced23`.

**Class:** Project-specific

Per 10-nt bin:

```text
balanced23_anchor_score = sqrt(23S × 23AS)
```

This is the geometric mean of the local sense and antisense 23-nt signals.

High values require both strands to be represented.

It is a **23-nt spatial anchor score**.

High values identify bins with concurrent 23S and 23AS signal for spatial analysis. They do not establish that the locus is primary, Dicer-derived, or causally upstream of 24-nt production.

## `combined23_anchor_score`

Previously `combined23`.

**Class:** Project-specific

```text
combined23_anchor_score = 23S + 23AS
```

Measures total local 23-nt signal without requiring strand balance.

This is also a spatial anchor score only; it does not assign biochemical pathway origin.

---

# 13. Historical anchor-selection parameters

These are parameters, not biological metrics.

```text
anchor_percentile        = 90th percentile
percentile_population    = non-zero anchor-score bins only
anchor_min_separation_nt = 50 nt
minimum_anchors          = 3
bin_size_nt              = 10 nt
```

## `anchor_threshold`

```text
anchor_threshold
    = percentile_90(anchor scores among bins with score > 0)
```

A candidate bin must satisfy:

```text
score >= anchor_threshold
```

Candidates are processed from strongest to weakest and greedily retained if at least 50 nt from already selected anchors.

Canonical implementation uses genomic coordinate as the deterministic tie-breaker for equal scores.

---

# 14. Anchor-window mean density

This is essential for interpreting all Stage-05 spatial metrics.

For each anchor and window `W`, the anchor bin itself is excluded. Available bins are collected separately upstream and downstream; reference boundaries truncate the window.

## `M_X_down(W)`

For track `X`:

```text
M_X_down(W)
    = sum of X over all valid anchor-specific downstream bin observations
      / number of valid anchor-specific downstream bin observations
```

## `M_X_up(W)`

```text
M_X_up(W)
    = sum of X over all valid anchor-specific upstream bin observations
      / number of valid anchor-specific upstream bin observations
```

If two anchor windows overlap, the same genomic bin can appear once for each anchor neighbourhood containing it.

Therefore these are **anchor-window pooled mean densities**, not means over unique genomic territory.

Canonical windows are:

```text
100 nt
250 nt
500 nt
```

With 10-nt bins these correspond to 10, 25, and 50 bins on each side.

---

# 15. Normalized 24-nt directionality

## `D_24AS`

**Class:** Project-specific historical/canonical endpoint component

```text
D_24AS
    = (M_24AS_down - M_24AS_up)
      / (M_24AS_down + M_24AS_up)
```

If the denominator is zero or either mean is not finite, report `NA`.

Range when defined with non-negative tracks:

```text
-1 to +1
```

Interpretation:

```text
+1  = signal only downstream
 0  = equal upstream/downstream mean density
-1  = signal only upstream
```

## `D_24S`

```text
D_24S
    = (M_24S_down - M_24S_up)
      / (M_24S_down + M_24S_up)
```

Same interpretation for the sense 24-nt control track.

## `antisense_specific_directionality`

Historical field:

```text
D24_antisense_minus_sense
```

**Class:** Project-specific spatial endpoint

```text
antisense_specific_directionality
    = D_24AS - D_24S
```

Possible range is approximately `-2` to `+2`.

Interpretation:

```text
>0 = 24-AS is more downstream-biased than 24-S
≈0 = little antisense-specific directional difference
<0 = 24-AS is less downstream-biased than 24-S
```

The subtraction is intended to control for generic positional asymmetry that affects both 24-nt strands.

This endpoint measures strand-controlled spatial directionality. It is not a direct classifier of primary versus secondary biogenesis or Dicer versus RdRP origin.

---

# 16. Antisense 23→24 composition

## `F24_AS_down`

**Class:** Standard proportion used in a project-specific spatial comparison

```text
F24_AS_down
    = M_24AS_down
      / (M_23AS_down + M_24AS_down)
```

## `F24_AS_up`

```text
F24_AS_up
    = M_24AS_up
      / (M_23AS_up + M_24AS_up)
```

If either denominator is zero, the corresponding fraction is `NA`.

Range:

```text
0 = antisense 23+24 signal is entirely 23 nt
1 = antisense 23+24 signal is entirely 24 nt
```

## `delta_F24_AS`

Historical field:

```text
downstream_minus_upstream_24_fraction_AS
```

**Class:** Project-specific spatial composition endpoint

```text
delta_F24_AS = F24_AS_down - F24_AS_up
```

Interpretation:

```text
>0 = downstream antisense signal is relatively more 24-nt dominated
≈0 = little upstream/downstream compositional change
<0 = downstream antisense signal is relatively more 23-nt dominated
```

A positive value does **not** necessarily mean there are more total small RNAs downstream. It means the antisense 23:24 composition shifted toward 24 nt.

For intuition, `+0.003` corresponds to a +0.3 percentage-point change in the 24-nt fraction.

---

# 17. Cross-correlation metrics

These are descriptive support, not the primary transitivity endpoints.

## `crosscorr_23_to_24(lag)`

**Class:** Standard correlation applied to project tracks

Historical v1.4.1 calculates Pearson correlation between:

```text
log1p(23-nt anchor-score track)
```

and:

```text
log1p(24-nt target-strand track)
```

at lags from `-500` to `+500 nt` in 10-nt steps.

Historical convention:

```text
positive lag = target 24-nt signal displaced downstream of the 23-nt anchor signal
```

A lag is `NA` when fewer than 10 overlapping bins remain or either vector has zero variance.

## `lag_asymmetry_strand`

```text
lag_asymmetry_strand
    = mean(correlation at positive lags)
      - mean(correlation at negative lags)
```

## `lag_asymmetry_AS_minus_S`

```text
lag_asymmetry_AS_minus_S
    = lag_asymmetry_antisense - lag_asymmetry_sense
```

These statistics describe directional correlation structure; they are not treated as proof of transitivity.

---

# 18. Historical circular-shift null

## `allowed_circular_shift`

**Class:** Historical project-specific null definition

For a track with `n` bins and exclusion distance `e` bins:

```text
shift ∈ {1, 2, ..., n-1}
```

is preferred when:

```text
min(shift, n - shift) > e
```

Under default parameters:

```text
e = max_window_nt / bin_size_nt = 500 / 10 = 50 bins
```

If no preferred shifts exist, historical v1.4.1 falls back to all non-zero circular shifts.

Within one contig and permutation replicate, **the same shift is applied to 24-AS and 24-S**.

The 23-nt tracks and anchors remain fixed.

The circular operation is a statistical randomization device used to preserve the internal spatial structure/autocorrelation of the shifted tracks while breaking their registration relative to the 23-nt anchors. It does not imply that the viral genome is biologically circular.

---

# 19. Empirical permutation P-value

## `p_shift`

**Class:** Standard Monte-Carlo permutation-P construction applied to the project null

For the pre-specified upper-tail alternative:

```text
p_shift
    = (b + 1) / (m + 1)
```

where:

```text
b = number of valid null statistics >= observed statistic
m = number of valid null statistics
```

The +1 correction prevents a finite random-permutation analysis from reporting `p = 0`.

The test direction must be fixed before viewing results.

---

# 20. Historical aggregation metrics

## `pair_balanced_median`

**Class:** Historical/project-specific summary

For a fixed weighting × anchor × window × endpoint:

```text
pair_balanced_median
    = median of finite sample-virus-contig endpoint values
```

This gives each row equal influence regardless of read depth.

It does **not** account for several virus rows originating from one sample.

## `virus_balanced_median`

**Class:** Historical/project-specific sensitivity summary

```text
virus_balanced_median
    = median across viruses of
      [median endpoint within each virus]
```

Purpose: reduce domination by viruses represented in many sample-virus units.

---

# 21. Canonical sample-balanced aggregation

## `sample_balanced_median`

**Class:** Canonical project summary

For a fixed weighting × anchor × window × endpoint:

```text
sample_stat(sample)
    = median of finite endpoint values across that sample's eligible virus-contigs

sample_balanced_median
    = median of sample_stat(sample) across samples
```

Purpose: prevent a sample containing several eligible viruses from behaving like several independent top-level samples.

## `sample_clustered_CI95`

**Class:** Canonical cluster-bootstrap uncertainty interval

Bootstrap sample IDs with replacement. All eligible observations belonging to a sampled sample are kept together, the within-sample median is recomputed, and then the across-sample median is recomputed.

Report:

```text
point estimate
2.5th percentile
97.5th percentile
number of bootstrap replicates
seed
```

---

# 22. Permutation aggregation

## `pair_global_shift_null`

Historical v1.4.1 takes, at each permutation index, the median of finite per-contig null statistics.

## `virus_balanced_global_shift_null`

At each permutation index:

1. median null statistic within each virus;
2. median across virus medians.

## `sample_balanced_global_shift_null`

**Class:** Canonical null aggregation

At each permutation index:

1. calculate each contig's shifted endpoint;
2. median across eligible virus-contigs within each sample;
3. median across sample medians.

Canonical `p_shift` compares the observed `sample_balanced_median` with this sample-balanced null distribution.

---

# 23. Benjamini-Hochberg adjusted values

## `q_BH_historical`

**Class:** Historical multiple-testing summary

Historical v1.4.1 applies BH across the three windows:

```text
100, 250, 500 nt
```

separately for each:

```text
weighting × anchor definition × endpoint × aggregation
```

Each historical family therefore contains three P-values.

## `q_BH_canonical`

**Class:** Canonical multiple-testing summary

For the primary sample-balanced observational analysis, one pre-specified inferential family is defined for each endpoint across:

```text
2 weighting modes × 2 anchor definitions × 3 windows = 12 tests
```

Separate 12-test families are used for:

```text
delta_F24_AS
antisense_specific_directionality
```

Robustness/descriptive analyses are not silently used as additional routes to a primary inferential claim.

---

# 24. Leave-one-virus-out analysis

## `leave_one_virus_out_effect`

**Class:** Sensitivity analysis

Repeatedly calculate the higher-level effect after excluding one virus.

Purpose:

> Determine whether the conclusion is heavily dependent on one virus.

A stable sign/magnitude across exclusions supports broader robustness; strong changes indicate virus dependence.

---

# 25. Interpretation hierarchy for Stage 05

Stage 05 is an **observational spatial analysis**. It does not assume that 23 nt is a primary/Dicer class or that 24 nt is a secondary/RdRP class.

The two main endpoints answer different questions.

## `antisense_specific_directionality`

Asks:

> Is 24-AS more downstream-biased than the 24-S control around predefined 23-nt spatial hotspots?

This is about **absolute spatial directionality after strand control**.

A positive value is consistent with antisense-specific downstream displacement of 24-nt signal. It does not establish the biochemical source of that signal.

## `delta_F24_AS`

Asks:

> Does the downstream antisense 23+24 population contain a larger fraction of 24-mers than the upstream population?

This is about **length composition**, not total abundance.

A positive value means a relative compositional shift toward 24 nt. It does not necessarily mean that total 24-nt abundance increased downstream.

A dataset can therefore show:

```text
positive delta_F24_AS
+
near-zero antisense_specific_directionality
```

without showing an absolute downstream wave of 24-AS molecules.

Such a pattern may be **consistent with** amplification/transitivity-associated biology, but it is not sufficient to prove RdRP-dependent secondary siRNA production or a 23→24 precursor-product relationship.

## Stage 05 ranking rule

Stage 05 metrics are analysis outputs only.

Do not automatically convert any of the following into an intrinsic vdCHIBIN window score:

```text
balanced23_anchor_score
combined23_anchor_score
antisense_specific_directionality
delta_F24_AS
crosscorr_23_to_24(lag)
lag_asymmetry_AS_minus_S
```

If the spatial analysis later informs design, it should first be considered at a regional/construct level and only after a separately documented decision.

---

# 26. Historical v1.4.1 regression values

These values are not metric definitions and must never be used as analytical inputs. They are stored only to verify exact historical replication.

For `unique_sequence × balanced23`, archived pair-balanced `delta_F24_AS` values are approximately:

```text
100 nt  -0.000113
250 nt  +0.002682
500 nt  +0.003662
```

Equivalent percentage-point shifts are approximately:

```text
100 nt  -0.0113 percentage points
250 nt  +0.2682 percentage points
500 nt  +0.3662 percentage points
```

Archived pair-level BH values for this endpoint are approximately:

```text
100 nt  0.471706
250 nt  0.018596
500 nt  0.000600
```

The archived `antisense_specific_directionality` endpoint did not show convincing evidence of an absolute downstream 24-AS wave.

---

# 27. Metric classes by biological role

## Stage 01 descriptive population metrics

```text
length_count(L)
length_fraction(L)
length_rank(L)
top1_indicator(L)
top3_indicator(L)
count_23
count_24
sense_fraction_23
antisense_fraction_23
sense_fraction_24
antisense_fraction_24
delta_antisense_fraction_24_minus_23
length23_fraction_among_23_24
length24_fraction_among_23_24
sample_balanced_median(metric)
sample_clustered_CI95
```

## Empirical sequence metrics

```text
observed_terminal_weight
observed_total_weight
observed_fraction
valid_background_window_count
expected_fraction_sense
expected_fraction_antisense
expected_fraction_combined
expected_fraction
enrichment_ratio
pair_median_enrichment_ratio
sample_enrichment_median
sample_balanced_median_enrichment_ratio
pooled_abundance_observed_fraction
pooled_abundance_expected_fraction
pooled_abundance_enrichment_ratio
spearman_rho_23_24
```

## Published duplex-geometry analysis metrics

```text
steprna_5p_distance
steprna_3p_distance
passenger_length
passenger_recovery_fraction_unique
passenger_recovery_fraction_abundance
steprna_log_ratio
steprna_wald_z
```

## Canonical Stage 04 population summaries

```text
sample_balanced_steprna_log_ratio
sample_balanced_joint_duplex_fraction(d5,d3)
sample_balanced_joint_00_duplex_fraction
joint_geometry_mode_by_pair
paired_delta_24_minus_23(M)
sample_balanced_varroa_2nt_joint_duplex_fraction
sample_balanced_varroa_2nt_reference_fraction_recovered
sample_balanced_varroa_2nt_reference_fraction_abundance_recovered
pair_balanced_geometry_summary
virus_balanced_geometry_summary
```

## Historical regression-only geometry statistic

```text
historical_delta_dicer
```

## Geometry-conditioned sequence metrics

```text
joint_observed_fraction
recovered_observed_fraction
E_joint_absolute
E_recovered_absolute
joint_vs_all_log2_contrast
joint_vs_recovered_log2_contrast
rho_joint_vs_general
rho_joint_contrast_abundance_vs_unique
```

These sequence metrics remain exploratory/candidate-development quantities. Stage 04 does not automatically turn them into a candidate-scoring dimension.

## Project-specific Stage 05 spatial/transitivity-consistency metrics

```text
balanced23_anchor_score
combined23_anchor_score
M_X_down / M_X_up
D_24AS
D_24S
antisense_specific_directionality
F24_AS_down
F24_AS_up
delta_F24_AS
crosscorr_23_to_24(lag)
lag_asymmetry_AS_minus_S
```

---

# 27A. Validated Stage 05 canonical results

Stage 05 is an observational spatial analysis and remains separate from vdCHIBIN per-window ranking.

## `delta_F24_AS` — validated canonical pattern

The canonical sample-balanced analysis found a reproducible positive `delta_F24_AS` at 250–500 nt across both weighting modes and both 23-nt anchor definitions.

Key estimates:

| Weighting | Anchor | Window | Estimate | 95% CI | BH P |
|---|---|---:|---:|---:|---:|
| abundance | balanced23 | 250 | +0.006185 | [-0.000957, +0.015907] | 0.037493 |
| abundance | balanced23 | 500 | +0.016160 | [+0.007511, +0.025209] | 0.007199 |
| abundance | combined23 | 250 | +0.005753 | [+0.001319, +0.011406] | 0.007199 |
| abundance | combined23 | 500 | +0.013266 | [+0.006500, +0.019799] | 0.007199 |
| unique_sequence | balanced23 | 250 | +0.003698 | [+0.002126, +0.006149] | 0.011998 |
| unique_sequence | balanced23 | 500 | +0.003344 | [+0.000876, +0.006022] | 0.004799 |
| unique_sequence | combined23 | 250 | +0.001928 | [-0.001638, +0.003904] | 0.010398 |
| unique_sequence | combined23 | 500 | +0.002153 | [+0.000332, +0.005789] | 0.005999 |

Interpretation:

```text
positive delta_F24_AS
=
downstream antisense 23+24 population is relatively more 24-nt-rich
```

It does **not** mean that total downstream 24AS abundance increased by the same percentage.

All eight 250/500-nt effects remained positive in leave-one-virus-out sensitivity analysis.

## `antisense_specific_directionality` — validated canonical pattern

No canonical combination supports a positive antisense-specific 24-nt downstream directionality effect after multiple-testing correction:

```text
BH P = 1.0 for all 12 pre-specified combinations
```

Therefore Stage 05 supports a **relative compositional shift**, not an absolute antisense-specific downstream 24-nt wave.

## Historical regression provenance

```text
historical_effect_size_regression = PASS
historical_permutation_regression = NOT_EXACTLY_REPRODUCED
historical_source_package_status = unavailable
historical_rng_stream_status = unavailable
historical_raw_p_checkpoint_status = unavailable
```

Historical permutation/BH values are provenance/regression information only. Canonical sample-balanced permutation inference is used for current biological interpretation.

---

# 28. Current biological interpretation

The clean canonical pipeline has established, in sequence:

```text
Stage 01
24 nt is the dominant 15–35-nt viral small-RNA length class;
23 nt is the second most prominent;
both are strongly antisense-biased;
24 nt is more strongly antisense-biased.

Stage 02
23- and 24-nt terminal nucleotide-enrichment landscapes
are highly concordant, with recurring features including
3p1 T/U enrichment and depletion of 5p1 G, 3p1 A and 3p2 A.

Stage 03
official stepRNA reconstructs complementary passengers in all
four focal classes. Distance 0 is a prominent marginal end-distance,
especially in antisense populations, but the same-duplex spectrum
shows that fully blunt (0,0) duplexes are only a minority.

The pre-specified (+2,-2) geometry is also a minority same-duplex
feature. A marginal 3p -2 component is reproducibly represented,
especially in antisense populations, while the matching 5p +2
component is substantially weaker.

Stage 04
sample-aware aggregation confirms that the geometry observations are
broadly reproducible across the dataset. Sample-balanced fully blunt
(0,0) fractions are approximately 2.15% (23S), 9.31% (23AS), 2.27%
(24S), and 3.11% (24AS), with 0/54 runs in every focal class having
>50% (0,0).

Sample-balanced (+2,-2) fractions are approximately 1.53% (23S),
2.00% (23AS), 1.98% (24S), and 0.58% (24AS). In antisense
populations, (+2,-2) support is greater for 23 nt than for 24 nt,
but this is not sufficient to classify the two length classes as
primary/Dicer versus secondary/RdRP products.

Geometry-conditioned terminal sequence effects are not consistently
strong across abundance and unique-sequence weighting and partly
overlap the Stage 02 terminal-enrichment landscape. A separate
geometry/Dicer ranking feature is therefore not carried forward.

Stage 05
14 samples and 19 eligible positive-sense sample-virus units across
BMLV, VDV-5, and VDV-9 were analysed spatially.

At 250–500 nt downstream of predefined 23-nt hotspots,
delta_F24_AS is reproducibly positive across abundance/unique
weighting and balanced/combined anchor definitions. All eight
pre-specified 250/500-nt combinations are BH-significant in the
canonical spatial permutation analysis, and all remain positive
when any one virus is left out.

This means the downstream antisense 23+24 population is modestly
more 24-nt-rich in composition.

In contrast, antisense_specific_directionality has no positive
BH-significant result (BH P = 1.0 for all 12 combinations).
Therefore Stage 05 does not show an absolute antisense-specific
downstream wave of 24-nt signal.

The combined result is consistent with a spatial
amplification/transitivity-associated pattern, but does not prove
a 23→24 precursor-product relationship, RdRP-dependent secondary
siRNA production, or Dicer-versus-RdRP pathway assignment.
```

Therefore the canonical project must **not currently describe 23 nt as proven primary/Dicer products or 24 nt as proven secondary/RdRP products**.

The pipeline must also preserve these distinctions:

```text
marginal distance 0
    ≠
fully blunt same-duplex (0,0)

positive delta_F24_AS
    ≠
absolute increase in downstream 24AS abundance

spatial 23→24 association
    ≠
proof of biochemical precursor→product order
```

The pipeline must never simplify the evidence to:

```text
all 23-mers = Dicer/primary
all 24-mers = secondary/RdRP
24-mers = Dicer-independent
marginal distance 0 = fully blunt duplex
(+2,-2) = a universal Dicer definition
positive Stage 05 = proven RdRP transitivity
500 nt = a universal biological propagation distance
Stage 05 spatial effect = intrinsic vdCHIBIN window score
```

Natural viral infection creates structured complementary RNA substrates, and viral replication, subgenomic transcription, RNA stability, library preparation, passenger recovery, and local sequence/mappability can all contribute to the observed population.

Stages 03–05 therefore provide biological/pathway context while remaining deliberately separate from per-window vdCHIBIN ranking.

---


# 29. Metrics deliberately not yet defined

Do not create these until later analyses justify them:

```text
S23_primary
S24_secondary
overall vdCHIBIN window score
duplex-geometry compatibility ranking weight
transitivity kernel K(d)
construct-level secondary score
```

---

# 30. Rule for introducing a new metric

Every future metric must document:

1. name;
2. mathematical definition;
3. whether it is standard, published-method-derived, project-specific, historical, or canonical;
4. input data;
5. analysis unit;
6. deduplication/weighting rule;
7. normalization;
8. biological question;
9. interpretation of high/low values;
10. zero/undefined handling;
11. uncertainty method;
12. limitations;
13. relationship to existing metrics;
14. whether the definition was fixed before inspecting the result.

No important statistic should exist only inside a Python script.

---

# 31. Main methodological references

- Murcott B, Pawluk RJ, Protasio AV, Akinmusola RY, Lastik D, Hunt VL. 2022. *stepRNA: Identification of Dicer cleavage signatures and passenger strand lengths in small RNA sequences*. Frontiers in Bioinformatics 2:994871. DOI: 10.3389/fbinf.2022.994871.
- Benjamini Y, Hochberg Y. 1995. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. Journal of the Royal Statistical Society Series B 57:289–300.
- Phipson B, Smyth GK. 2010. *Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn*. Statistical Applications in Genetics and Molecular Biology 9:Article 39.
- Saravanan V, Berman GJ, Sober SJ. 2020. *Application of the hierarchical bootstrap to multi-level data in neuroscience*. General methodological support for clustered/nested resampling.
- The uploaded Varroa vsiRNA v1.4.1 strengthened-transitivity code and archived outputs are the exact historical source for the Stage-05 coordinate, anchor, window, weighting, null, and aggregation definitions reproduced here.
