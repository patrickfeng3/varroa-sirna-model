# Varroa vsiRNA Metric Dictionary

**Version:** 0.3  
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

For Stage 05, identity is exactly:

```text
virus × strand × length × sequence
```

and total sequence weight 1 is divided across all exact compatible loci for that sequence.

---

# 4. Basic 23/24 population metrics

## `count_23`

**Class:** Standard/descriptive

Number or weighted abundance of eligible 23-nt viral small RNAs under the stated weighting mode.

Always retain sample, virus, strand, and weighting mode.

## `count_24`

Equivalent quantity for eligible 24-nt viral small RNAs.

## `antisense_fraction`

**Class:** Standard/descriptive

```text
antisense_fraction = antisense / (sense + antisense)
```

Interpretation:

```text
0.5  = equal sense and antisense contribution
>0.5 = antisense-biased
<0.5 = sense-biased
```

If the denominator is zero, report `NA`.

## `sense_fraction`

```text
sense_fraction = sense / (sense + antisense)
```

When defined:

```text
sense_fraction + antisense_fraction = 1
```

---

# 5. Terminal nucleotide coordinates

Terminal positions are defined relative to the **physical sequenced RNA in its own 5′→3′ orientation**:

```text
5′ N1 N2 ................ N(n-1) Nn 3′
   ↑  ↑                    ↑     ↑
 5p1 5p2                  3p2   3p1
```

Antisense RNAs are therefore evaluated in antisense RNA orientation, not by simply reading the viral reference left-to-right.

---

# 6. Terminal nucleotide enrichment

## `observed_fraction(b,p)`

**Class:** Standard/descriptive

For nucleotide `b` at terminal position `p`:

```text
observed_fraction(b,p)
    = number/weight of eligible observed RNAs with b at p
      / total number/weight of eligible observed RNAs
```

## `expected_fraction(b,p)`

**Class:** Project-specific matched-background quantity

Frequency of nucleotide `b` at position `p` among all fully depth-supported windows of the **same RNA length** in the corresponding sample-specific viral background.

- Sense expectation uses reference orientation.
- Antisense expectation uses reverse-complement orientation.
- Combined expectation uses the observed strand mixture rather than forcing a 50:50 mixture.

Conceptually:

```text
expected_combined
    = wS × expected_sense
      + wAS × expected_antisense
```

with `wS + wAS = 1`.

## `enrichment_ratio`

**Class:** Project-specific empirical effect size

```text
enrichment_ratio = observed_fraction / expected_fraction
```

Interpretation:

```text
1   = observed as often as sequence availability predicts
>1  = enriched
<1  = depleted
```

If `expected_fraction = 0`, report `NA`.

This metric does not identify which molecular process created the enrichment.

## `pair_median_enrichment_ratio`

**Class:** Project-specific historical/cross-pair summary

Median `enrichment_ratio` across eligible sample-virus units.

This corresponds most closely to the historical design-facing `median_enrichment_ratio` used in the project.

## `sample_balanced_median_enrichment_ratio`

**Class:** Project-specific canonical summary

1. median enrichment across eligible viruses within each sample;
2. median across sample-level medians.

Use a sample-clustered bootstrap for its confidence interval.

Both pair- and sample-balanced fields should be exported so historical design references are not silently redefined.

## `spearman_rho_23_24`

**Class:** Standard statistical metric

Spearman rank correlation between matched 23- and 24-nt terminal enrichment landscapes.

```text
+1 = identical rank order
 0 = no monotonic association
-1 = opposite rank order
```

This measures similarity of enrichment patterns, not shared enzymatic origin.

---

# 7. stepRNA geometry metrics

The following quantities follow official stepRNA conventions.

## `steprna_5p_distance`

**Class:** Published-method-derived

Signed 5′ distance relative to the File-A reference RNA:

```text
negative = reference overhang
positive = reference underhang
0        = blunt
```

## `steprna_3p_distance`

**Class:** Published-method-derived

Equivalent signed distance at the 3′ end using the same sign convention.

## `passenger_length`

**Class:** Published-method-derived

Length of the complementary passenger sequence identified by stepRNA.

## `passenger_recovery_fraction`

**Class:** Descriptive quantity based on stepRNA output

```text
passenger_recovery_fraction
    = File-A references with at least one recovered passenger
      / all eligible File-A references
```

Low passenger recovery does not by itself mean Dicer processing is absent.

## `steprna_log_ratio` / installed stepRNA enrichment field

**Class:** Published-method-derived

Use the enrichment/log-ratio field produced by the installed official stepRNA version. Do not silently reimplement it with a different formula or rename it as a generic odds ratio unless the software output explicitly uses that definition.

The published method describes a log ratio comparing a distance-specific count with the mean end-distance count.

## `steprna_wald_z`

**Class:** Published-method-derived

Official stepRNA Wald Z-score for enrichment of an overhang/underhang distance.

Use the value produced by stepRNA rather than creating a different statistic under the same name.

---

# 8. Pre-specified Varroa joint Dicer-like geometry

The term “canonical Dicer geometry” is avoided as a universal label because Dicer-associated end geometry varies by pathway and organism.

## `varroa_2nt_joint_geometry`

**Class:** Project-specific, pre-specified pathway feature

Historical Varroa geometry of interest:

```text
5′ distance = +2
3′ distance = -2
```

under official stepRNA sign convention.

Equivalent historical label:

```text
5p_underhang_2__3p_overhang_2
```

This is a pre-specified Varroa feature of interest, not a universal definition of Dicer cleavage.

## `varroa_2nt_fraction_all_refs`

```text
number of all eligible File-A references supporting the joint geometry
/ all eligible File-A references
```

This combines passenger recoverability and geometry and must not be interpreted alone.

## `varroa_2nt_fraction_recovered`

```text
number of passenger-recovered references supporting the joint geometry
/ number of File-A references with at least one recovered passenger
```

This asks how common the geometry is **conditional on a passenger actually being recoverable**.

---

# 9. `delta_dicer`

Written `Δ_Dicer`.

**Class:** Project-specific secondary statistic

For pre-specified target distance `d*` and pre-specified comparison set `D0`:

```text
Δ_Dicer
    = support(d*)
      - mean[support(d) for d in D0]
```

`D0` must be defined before examining the result.

Interpretation:

```text
≈0 = target geometry does not stand out
>0 = target geometry exceeds the comparison distances
```

This is not an official stepRNA statistic and is not directly a candidate-window score.

---

# 10. Dicer-conditioned sequence metrics

## `E_Dicer_absolute(f)`

**Class:** Project-specific exploratory/candidate-development metric

```text
E_Dicer_absolute(f)
    = frequency of feature f among the pre-specified Dicer-supported subset
      / matched viral-sequence expected frequency of f
```

The background must be matched by sample, virus, length, and strand.

## `E_all_observed(f)`

General matched enrichment for feature `f` from Stage 02.

## `dicer_specific_log2_contrast(f)`

**Class:** Project-specific exploratory/candidate-development metric

```text
dicer_specific_log2_contrast(f)
    = log2(E_Dicer_absolute(f) / E_all_observed(f))
```

Interpretation:

```text
0  = Dicer-supported subset adds little beyond general enrichment
>0 = feature is more enriched in Dicer-supported RNAs
<0 = feature is less enriched in Dicer-supported RNAs
```

If the required ratio is undefined or non-positive, report `NA`; do not add an arbitrary pseudocount solely to force a finite logarithm.

## `dicer_general_correlation`

**Class:** Standard correlation applied to project metrics

Correlation between general terminal enrichment and Dicer-conditioned enrichment.

Purpose: detect redundancy before any Dicer-derived feature is considered for later sequence ranking.

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

It is a **primary-like anchor score**, not proof that the locus is biologically primary.

## `combined23_anchor_score`

Previously `combined23`.

**Class:** Project-specific

```text
combined23_anchor_score = 23S + 23AS
```

Measures total local 23-nt signal without requiring strand balance.

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

**Class:** Project-specific primary/secondary endpoint

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

For the primary sample-balanced analysis, one family is defined for each biological endpoint across:

```text
2 weighting modes × 2 anchor definitions × 3 windows = 12 tests
```

Separate 12-test families are used for:

```text
delta_F24_AS
antisense_specific_directionality
```

Robustness/descriptive analyses are not silently used as additional routes to a primary significance claim.

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

The two main biological endpoints answer different questions.

## `antisense_specific_directionality`

Asks:

> Is 24-AS more downstream-biased than the 24-S control?

This is about **absolute spatial directionality after strand control**.

## `delta_F24_AS`

Asks:

> Does the downstream antisense 23+24 population contain a larger fraction of 24-mers than the upstream population?

This is about **length composition**, not total abundance.

A dataset can show a positive composition shift while showing little or no antisense-specific absolute directionality. That pattern is exactly why both metrics must be kept separate.

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

## Descriptive population metrics

```text
count_23
count_24
sense_fraction
antisense_fraction
```

## Empirical sequence metrics

```text
observed_fraction
expected_fraction
enrichment_ratio
pair_median_enrichment_ratio
sample_balanced_median_enrichment_ratio
spearman_rho_23_24
```

## Published Dicer-analysis metrics

```text
steprna_5p_distance
steprna_3p_distance
passenger_length
passenger_recovery_fraction
installed stepRNA enrichment/log-ratio output
steprna_wald_z
```

## Project-specific Dicer summaries

```text
varroa_2nt_joint_geometry
varroa_2nt_fraction_all_refs
varroa_2nt_fraction_recovered
Δ_Dicer
```

## Potential Dicer-derived candidate features

```text
E_Dicer_absolute
dicer_specific_log2_contrast
dicer_general_correlation
```

These remain exploratory until reproducibility and non-redundancy are demonstrated.

## Project-specific spatial/transitivity metrics

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

# 28. Current biological interpretation

The working population-level model is:

```text
23 nt
→ more strongly associated with the primary/Dicer-like response

24 nt antisense
→ associated with a later/secondary-like population
→ while still showing evidence compatible with Dicer/Dicer-like processing
```

These are population associations, not identities.

The pipeline must never simplify this to:

```text
all 23-mers = primary
all 24-mers = secondary
24-mers = Dicer-independent
```

Natural viral infection also creates its own spatially structured RNA substrates, so Stage 05 provides evidence **consistent with** secondary/transitive biology rather than direct mechanistic proof.

---

# 29. Metrics deliberately not yet defined

Do not create these until later analyses justify them:

```text
S23_primary
S24_secondary
overall vdCHIBIN window score
Dicer-compatibility ranking weight
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
