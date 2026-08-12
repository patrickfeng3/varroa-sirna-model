# Varroa vsiRNA Metric Dictionary

**Version:** 0.2
**Scope:** Canonical viral pipeline through viral spatial/transitivity-consistency analysis

---

# 1. Purpose

This document defines the important quantities used by the canonical Varroa viral small-RNA pipeline.

Every metric is classified as one of:

* **Standard/descriptive** — conventional mathematical or statistical quantity.
* **Published-method-derived** — produced by or directly based on a published method.
* **Project-specific** — developed for this Varroa analysis for a defined biological question.
* **Provisional** — conceptually defined, but exact implementation must be locked before canonical coding.

A project-specific metric is not automatically a weak metric. It simply means the exact statistic was designed for this biological question rather than copied from a published software package.

---

# 2. Analysis units

## `sample`

One biological small-RNA sequencing library.

This is the main biological clustering level for uncertainty estimation.

---

## `sample_virus_unit`

One virus analysed within one sample.

Example:

```text
sample A × DWV-A
sample A × DWV-B
sample B × DWV-A
```

The first two units share the same biological library and therefore should not automatically be treated as fully independent replicates.

---

# 3. Weighting modes

## `abundance_weighted`

**Class:** Standard/descriptive

Repeated sequencing observations retain their abundance.

Example:

```text
sequence A = 1000 reads
sequence B =   10 reads
```

A contributes 100 times as much as B.

### Question answered

> What dominates the accumulated sequenced small-RNA population?

---

## `unique_sequence`

**Class:** Standard/descriptive

Each distinct sequence contributes once within a precisely defined analysis unit.

Example:

```text
sequence A = 1000 reads → 1 sequence
sequence B =   10 reads → 1 sequence
```

### Question answered

> Is the pattern represented across many distinct small-RNA sequences rather than being driven by a few abundant sequences?

### Required implementation rule

The deduplication unit must always be stated.

For example:

```text
sample × virus × length × strand
```

or, for spatial analysis, the exact corresponding spatial unit.

The earlier viral-transitivity implementation contained a bug where a nominal unique-sequence analysis remained effectively read-weighted. The canonical implementation must perform true sequence-level deduplication.

---

# 4. Basic 23/24 population metrics

## `count_23`

**Class:** Standard/descriptive

Number of eligible 23-nt viral small RNAs under the specified weighting mode.

Always report the associated:

* sample
* virus
* strand
* weighting mode

---

## `count_24`

Equivalent quantity for eligible 24-nt viral small RNAs.

---

## `antisense_fraction`

**Class:** Standard/descriptive

### Formula

```text
antisense_fraction
=
antisense
/
(sense + antisense)
```

### Range

```text
0 to 1
```

### Interpretation

```text
0.5  → balanced
>0.5 → antisense-biased
<0.5 → sense-biased
```

### Undefined case

If:

```text
sense + antisense = 0
```

report:

```text
NA
```

Do not introduce a pseudocount merely to create a value.

---

## `sense_fraction`

```text
sense_fraction
=
sense
/
(sense + antisense)
```

When both are defined:

```text
sense_fraction + antisense_fraction = 1
```

---

# 5. Terminal nucleotide coordinates

Terminal positions are defined relative to the **physical sequenced RNA in its own 5′→3′ orientation**.

```text
5′ N1 N2 ................ N(n-1) Nn 3′
   ↑  ↑                    ↑     ↑
 5p1 5p2                  3p2   3p1
```

Therefore an antisense RNA must be interpreted in its antisense RNA orientation rather than simply reading reference-genome coordinates.

---

# 6. Terminal nucleotide enrichment

## `observed_fraction`

**Class:** Standard/descriptive

For nucleotide (b) at terminal position (p):

```text
observed_fraction(b,p)
=
number of observed eligible RNAs with b at p
/
number of eligible observed RNAs
```

Example:

```text
40 of 100 observed 23-mers contain U at 5p1

observed_fraction(U,5p1) = 0.40
```

---

## `expected_fraction`

**Class:** Project-specific matched-background quantity

The expected frequency of that nucleotide among all **sequence-available, fully depth-supported windows of the same length** in the corresponding sample-specific viral consensus.

This corrects for viral sequence composition.

The existing validated pipeline generates depth-masked consensuses specifically so poorly supported reference positions do not dominate this background.

### Sense

Use the viral reference orientation.

### Antisense

Use the reverse-complement orientation.

### Combined strand

The expected background is weighted according to the observed strand mixture:

```text
expected_combined
=
wS × expected_sense
+
wAS × expected_antisense
```

where:

```text
wS + wAS = 1
```

Do not automatically assume a 50:50 strand mixture.

---

## `enrichment_ratio`

**Class:** Project-specific empirical effect size

### Formula

```text
enrichment_ratio
=
observed_fraction
/
expected_fraction
```

### Interpretation

```text
1    → exactly as frequent as expected
>1   → enriched
<1   → depleted
```

Example:

```text
observed 5p1-U = 0.40
expected 5p1-U = 0.20

enrichment_ratio = 2.0
```

### Undefined case

If:

```text
expected_fraction = 0
```

report:

```text
NA
```

rather than infinity or an arbitrary pseudocount.

### Biological interpretation

This is an empirical preference among sequenced Varroa viral small RNAs.

It does not independently identify whether enrichment arose from:

* Dicer processing
* Argonaute loading
* strand selection
* RNA degradation/stability
* sequencing/library effects
* another biological process

---

## `median_enrichment_ratio`

**Class:** Project-specific cross-dataset summary

Median `enrichment_ratio` across eligible biological units.

This remains the main empirical nucleotide statistic intended for later design work.

### Why median?

The median is less sensitive than a mean to one extreme sample-virus infection.

### Required accompanying information

Always retain:

```text
median
95% sample-aware bootstrap CI
number of contributing samples
number of contributing sample-virus units
```

---

# 7. 23-vs-24 similarity

## `spearman_rho_23_24`

**Class:** Standard statistical metric

Spearman rank correlation between matched 23- and 24-nt enrichment landscapes.

### Range

```text
+1 → identical rank order
 0 → no monotonic association
-1 → opposite rank order
```

### Purpose

Tests whether nucleotide features enriched among 23-mers tend also to be enriched among 24-mers.

This is an association statistic, not evidence that the same enzyme generated both populations.

---

# 8. stepRNA geometry quantities

stepRNA is the published primary method for Dicer-overhang analysis.

It directly aligns candidate reference small RNAs against potential passenger RNAs and determines 5′ and 3′ overhang/underhang geometry and passenger length. It uses exact matching by default.

---

## `steprna_5p_distance`

**Class:** Published-method-derived

Signed 5′ duplex-end distance reported according to the official stepRNA convention.

stepRNA defines:

```text
negative = overhang
positive = underhang
0        = blunt
```

relative to the File-A reference strand.

Downstream code must preserve this sign convention exactly.

---

## `steprna_3p_distance`

**Class:** Published-method-derived

Equivalent signed distance at the 3′ end.

Again:

```text
negative = overhang
positive = underhang
0        = blunt
```

according to stepRNA's reference-strand convention.

---

## `passenger_length`

**Class:** Published-method-derived

Length of the complementary passenger sequence identified by stepRNA.

stepRNA explicitly reports passenger-length distributions as one of its primary outputs.

---

# 9. Passenger recovery

This should be kept separate from Dicer geometry.

## `passenger_recovery_fraction`

**Class:** Published-method-derived descriptive quantity

### Formula

```text
number of File-A reference RNAs
with ≥1 predicted passenger
/
total number of eligible File-A reference RNAs
```

### Interpretation

Higher values mean complementary passenger sequences can be reconstructed for a greater proportion of the population.

### Important limitation

Low passenger recovery does not imply absence of Dicer processing.

Passenger strands can be:

* degraded
* depleted during RISC maturation
* poorly represented in sequencing
* rare because of strong strand bias

stepRNA itself demonstrates that a Dicer signature can still be detected when only a small fraction of non-collapsed reference RNAs have recoverable passengers.

---

# 10. Varroa pre-specified 2-nt geometry

The term **“canonical Dicer support” should no longer be used without qualification.**

Dicer frequently produces approximately 2-nt 3′ overhangs, but Dicer cleavage geometry can vary with substrate and biological pathway. Experimental and stepRNA analyses show enriched 0-, 1-, 2- or 3-nt geometries in different contexts.

Therefore use:

## `varroa_2nt_joint_geometry`

**Class:** Project-specific, pre-specified pathway feature

Historical Varroa representation:

```text
5p_underhang_2__3p_overhang_2
```

Under official stepRNA sign convention this corresponds conceptually to:

```text
5′ distance = +2
3′ distance = -2
```

relative to the selected File-A reference strand.

This is a **pre-specified Varroa geometry of interest**, not a universal definition of Dicer cleavage.

---

## `varroa_2nt_fraction_all_refs`

### Formula

```text
number of File-A references supporting
the pre-specified joint geometry
/
all eligible File-A references
```

### Purpose

Measures recoverable population-level support.

### Limitation

Strongly affected by passenger availability.

It should therefore never be interpreted alone.

---

## `varroa_2nt_fraction_recovered`

### Formula

```text
number of passenger-recovered reference RNAs
supporting the pre-specified joint geometry
/
number of reference RNAs with ≥1 recovered passenger
```

### Purpose

Asks:

> Among references for which a complementary partner was actually recoverable, how common is the Varroa 2-nt geometry?

This separates geometry from the overall passenger-recovery problem.

Both `varroa_2nt_fraction_all_refs` and `varroa_2nt_fraction_recovered` should be reported.

---

# 11. Official stepRNA enrichment metrics

## `steprna_log_odds`

**Class:** Published-method-derived

Use the value generated by official stepRNA.

stepRNA calculates a log ratio comparing the count at a particular end-distance with the mean count across end distances.

Do not silently recreate this statistic with a different formula.

---

## `steprna_wald_z`

**Class:** Published-method-derived

Official stepRNA significance statistic for overhang-distance enrichment.

stepRNA calculates Z-scores using a Wald-test framework.

### Interpretation

Larger positive Z:

> stronger enrichment of the tested geometry relative to stepRNA's internal background model.

### Use

This is the **primary published statistical evidence** for enriched Dicer-like duplex geometry.

---

# 12. `delta_dicer`

Written:

```text
Δ_Dicer
```

**Class:** Project-specific secondary statistic

This is **not an official stepRNA statistic**.

### Definition

For a pre-specified target end-distance (d^*) and a pre-specified comparison set (D_0):

```text
Δ_Dicer
=
support(d*)
-
mean[support(d) for d in D0]
```

### Critical requirement

The comparison-distance set `D0` must be fixed in configuration **before analysing the results**.

It cannot be chosen after inspecting which alternative distances are low.

### Interpretation

```text
≈0 → target geometry does not stand out
>0 → target geometry exceeds comparison distances
```

### Role

Secondary validation of our historical Varroa result.

Official stepRNA inference remains primary.

---

# 13. Dicer-conditioned sequence features

This analysis asks whether the subset with strong Dicer-like geometry has nucleotide properties beyond those already present in all observed siRNAs.

---

## `E_Dicer_absolute`

**Class:** Project-specific exploratory/candidate-development metric

For feature (f):

```text
E_Dicer_absolute(f)
=
frequency of f among the pre-specified
Dicer-supported subset
/
expected frequency of f from the matched
viral-sequence background
```

### Interpretation

```text
1    → as frequent as sequence availability predicts
>1   → enriched
<1   → depleted
```

The background must be matched by:

* sample
* virus
* length
* strand

in the same manner as Stage 02.

---

## `E_all_observed`

General matched enrichment for the same feature from Stage 02.

---

## `dicer_specific_log2_contrast`

**Class:** Project-specific exploratory/candidate-development metric

Preferred over the previous raw ratio-of-ratios because it is symmetric around zero.

### Formula

```text
dicer_specific_log2_contrast
=
log2(
    E_Dicer_absolute
    /
    E_all_observed
)
```

### Interpretation

```text
0  → Dicer subset adds no enrichment beyond the general population

>0 → feature is more enriched among Dicer-supported RNAs

<0 → feature is less enriched among Dicer-supported RNAs
```

### Undefined cases

If either required enrichment is mathematically undefined or zero in a way that makes the logarithm undefined:

```text
report NA
```

Do not introduce an arbitrary pseudocount merely to create a finite value.

A formal count-based model may later be preferable if sparse cells become important.

---

## `dicer_general_correlation`

**Class:** Standard statistical comparison applied to project metrics

Spearman correlation between:

```text
general terminal enrichment
```

and:

```text
Dicer-conditioned terminal enrichment
```

### Purpose

Determine whether the proposed Dicer-derived metric provides information independent of the existing nucleotide-enrichment score.

A high correlation would argue against automatically treating both as independent ranking dimensions.

---

# 14. Dicer metric interpretation

Dicer metrics fall into two fundamentally different groups:

### Pathway-level

```text
stepRNA distance distribution
steprna_wald_z
passenger_recovery_fraction
varroa_2nt_fraction_all_refs
varroa_2nt_fraction_recovered
Δ_Dicer
```

These describe **how a population appears to have been processed**.

They are not intrinsic scores for an untested candidate sequence.

### Candidate-feature-development

```text
E_Dicer_absolute
dicer_specific_log2_contrast
```

These may eventually become candidate-level features **only if they are reproducible and non-redundant**.

---

# 15. Viral spatial/transitivity metrics

Secondary-siRNA/transitivity literature supports examining small RNAs arising beyond or along the target transcript after primary silencing and demonstrates that secondary-siRNA production can spread directionally from an initiating event.

However, the exact metrics below are **Varroa project-specific**, not standard published transitivity statistics.

Our previous v1.4.1 analysis used `balanced23`, `combined23`, `D_24AS − D_24S`, `F24_AS`, true sequence deduplication, paired movement of 24S/24AS during permutation, and BH correction.

The exact spatial hotspot implementation still needs to be recovered before Stage 05 is frozen.

---

# 16. `balanced23_anchor_score`

Previously called:

```text
balanced23
```

**Class:** Project-specific, provisional until Stage-05 implementation lock

### Formula

```text
balanced23_anchor_score
=
sqrt(23_sense × 23_antisense)
```

This is the geometric mean of the two strand-specific local signals.

### Interpretation

High values require substantial signal from **both strands**.

Example:

```text
23S  = 100
23AS = 100

balanced23 = 100
```

If one strand approaches zero, the score approaches zero.

### Biological purpose

Identify regions more compatible with bidirectional processing of dsRNA than a purely one-sided hotspot.

### Important note

This is **not a published standard RNAi metric**.

It is a biologically motivated anchor score created for this analysis.

The exact underlying track unit—raw abundance, normalized abundance or unique-sequence signal—must be specified separately.

---

# 17. `combined23_anchor_score`

Previously:

```text
combined23
```

**Class:** Project-specific, provisional

### Formula

```text
combined23_anchor_score
=
23_sense + 23_antisense
```

### Purpose

Measures total local 23-nt activity without requiring strand balance.

Used as a looser alternative to `balanced23_anchor_score`.

---

# 18. `F24_AS`

**Class:** Standard proportion used in a project-specific biological comparison

### Formula

```text
F24_AS
=
24AS
/
(23AS + 24AS)
```

### Meaning

Among antisense 23+24-nt small RNAs in the region, what proportion are 24 nt?

### Range

```text
0 → entirely 23 nt
1 → entirely 24 nt
```

### Zero-signal rule

If:

```text
23AS + 24AS = 0
```

then:

```text
F24_AS = NA
```

There is no composition to estimate.

No pseudocount should manufacture a value.

---

# 19. `delta_F24_AS`

Written:

```text
ΔF24_AS
```

**Class:** Project-specific spatial contrast

### Formula

```text
ΔF24_AS
=
F24_AS_downstream
-
F24_AS_upstream
```

### Interpretation

```text
>0
downstream antisense population is relatively more 24-nt dominated

≈0
little change in 23/24 composition

<0
downstream population is relatively more 23-nt dominated
```

### Important distinction

A positive value does **not** necessarily mean more total small RNA is produced downstream.

It describes **composition**.

---

# 20. `D_24AS`

**Class:** Project-specific, provisional until Stage-05 implementation lock

Conceptually:

```text
D_24AS
=
24AS_downstream
-
24AS_upstream
```

However, the exact definition of the underlying downstream/upstream signal—particularly:

* normalization
* boundary handling
* hotspot aggregation
* positional multimapping

must be recovered from the final v1.4.1 implementation before this metric is considered fully canonical.

Do not code from this conceptual definition alone.

---

# 21. `D_24S`

Equivalent directional contrast for 24-nt sense RNA.

Also **provisional** until the exact Stage-05 spatial implementation is recovered.

---

# 22. `antisense_specific_directionality`

Historical notation:

```text
D_24AS − D_24S
```

**Class:** Project-specific difference-of-directionality contrast

### Conceptual formula

```text
antisense_specific_directionality
=
D_24AS
-
D_24S
```

### Purpose

Control for general downstream-versus-upstream effects that influence both 24-nt strands.

Question answered:

> Is downstream behaviour stronger specifically for the 24-nt antisense population than for the corresponding 24-nt sense population?

### Interpretation

```text
>0 → more antisense-specific downstream directionality
≈0 → no strong antisense-specific effect
<0 → opposite directional tendency
```

This remains provisional until `D_24AS` and `D_24S` are implementation-locked.

---

# 23. Canonical spatial distances

```text
100 nt
250 nt
500 nt
```

**Class:** Project analysis parameters

These are predefined comparison distances.

They are **not metrics** and must not be interpreted as established biological propagation distances.

---

# 24. Permutation P-value

## `p_permutation`

**Class:** Standard resampling inference

For (m) randomly generated null permutations and (b) null statistics at least as extreme as the observed statistic under the pre-specified tail:

```text
p_permutation
=
(b + 1)
/
(m + 1)
```

The +1 correction prevents reporting an impossible zero P-value when only a finite number of random permutations were sampled, following recommended practice for Monte-Carlo permutation tests.

### Requirement

Before analysis, specify whether the test is:

* upper-tailed
* lower-tailed
* two-sided

according to the biological hypothesis.

Do not choose the tail after seeing the result.

---

# 25. Benjamini-Hochberg adjusted value

## `p_BH`

or preferably:

```text
q_BH
```

**Class:** Standard multiple-testing procedure

Benjamini-Hochberg controls the false-discovery rate across a defined family of hypothesis tests.

### Requirement

The test family must be defined before examining significance.

Examples might involve the predefined combination of:

* endpoint
* distance
* anchor definition
* weighting mode

The exact Stage-05 family from v1.4.1 still needs to be recovered.

---

# 26. Sample-clustered bootstrap

## `sample_clustered_CI95`

**Class:** Standard hierarchical-resampling approach adapted to this dataset

### Procedure

The biological sample is the top-level resampling cluster.

When a sample is selected during bootstrap resampling, its relevant sample-virus observations are retained together.

### Why?

Multiple observations originating from the same sequencing library are correlated rather than independent.

Hierarchical/clustered bootstrap methods are specifically designed for nested data where lower-level observations share higher-level biological units.

### Output

Report:

```text
point estimate
2.5th bootstrap percentile
97.5th bootstrap percentile
number of bootstrap replicates
random seed
```

The exact CI method must be recorded.

---

# 27. Pair-balanced summary

**Class:** Project-specific descriptive/sensitivity aggregation

Each eligible sample-virus unit contributes comparably instead of weighting units by sequencing depth.

Useful for determining whether extremely deep infections dominate the result.

This is a sensitivity/descriptive analysis, not a substitute for respecting sample-level clustering.

---

# 28. Virus-balanced summary

**Class:** Project-specific sensitivity aggregation

Each biological virus contributes comparably to the higher-level summary.

### Purpose

Determine whether the overall conclusion is dominated by viruses represented in many samples or with unusually strong signal.

---

# 29. Leave-one-virus-out analysis

**Class:** Standard sensitivity-analysis principle applied to this project

Repeatedly calculate the result after excluding one biological virus.

### Interpretation

Stable result:

> stronger evidence that the conclusion is not driven by one virus.

Strongly changing result:

> conclusion depends substantially on particular virus biology.

---

# 30. Metric classes

## Descriptive population metrics

```text
count_23
count_24
sense_fraction
antisense_fraction
```

---

## Empirical sequence metrics

```text
observed_fraction
expected_fraction
enrichment_ratio
median_enrichment_ratio
spearman_rho_23_24
```

---

## Published Dicer-analysis metrics

```text
steprna_5p_distance
steprna_3p_distance
passenger_length
passenger_recovery_fraction
steprna_log_odds
steprna_wald_z
```

---

## Project-specific Dicer summaries

```text
varroa_2nt_joint_geometry
varroa_2nt_fraction_all_refs
varroa_2nt_fraction_recovered
Δ_Dicer
```

---

## Potential Dicer-derived candidate features

```text
E_Dicer_absolute
dicer_specific_log2_contrast
dicer_general_correlation
```

These remain exploratory until reproducibility and independence are demonstrated.

---

## Project-specific spatial/transitivity metrics

```text
balanced23_anchor_score
combined23_anchor_score
F24_AS
ΔF24_AS
D_24AS
D_24S
antisense_specific_directionality
```

`D_24AS`, `D_24S` and their contrast remain provisional until the final v1.4.1 spatial implementation is recovered.

---

# 31. Current biological interpretation

Current results support the cautious population-level interpretation:

```text
23 nt
→ more strongly primary/Dicer-associated

24 nt
→ more antisense-biased and secondary-associated
→ but also carries evidence of Dicer/Dicer-like processing
```

The previous Varroa analyses found significant Dicer-like geometry in both length populations, with stronger evidence for the 23-nt class, and a modest rather than dramatic viral spatial shift toward 24-nt antisense RNA.

This does not mean:

```text
all 23-mers are primary
all 24-mers are secondary
24-mers are Dicer-independent
```

---

# 32. Metrics deliberately not yet defined

Do not create these until the empirical evidence supports them:

```text
S23_primary
S24_secondary
overall vdCHIBIN score
Dicer-compatibility ranking weight
transitivity kernel K(d)
construct-level secondary score
```

The canonical viral pipeline should first establish the underlying biology.

---

# 33. Rule for introducing future metrics

Every new metric must document:

1. name
2. mathematical definition
3. whether it is standard, published-method-derived or project-specific
4. input data
5. analysis unit
6. deduplication/weighting rule
7. normalization
8. biological question
9. interpretation
10. undefined/zero cases
11. uncertainty method
12. limitations
13. relationship to existing metrics
14. whether the definition was fixed before examining results

No important statistic should exist only inside a Python script.

---

# 34. Primary methodological sources

The metric framework is grounded in:

* Murcott et al. (2022), **stepRNA**, for computational identification of Dicer cleavage geometry and passenger strands.
* Experimental and structural Dicer literature showing that Dicer commonly produces short 3′ overhangs but that exact cleavage geometry can vary with substrate.
* Hierarchical/bootstrap methodology for nested biological data.
* Benjamini-Hochberg FDR control for related multiple hypothesis tests.
* Monte-Carlo permutation-test methodology using non-zero corrected empirical P-values.
* Secondary-siRNA/transitivity literature for the biological concept of secondary small-RNA spread.

Project-specific spatial endpoints are explicitly identified as such rather than presented as published standard statistics.
