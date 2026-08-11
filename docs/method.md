# Method: Nascent residual and mechanism annotation

```{admonition} Scope of this page
:class: note

**Software mathematics** for the reference-corrected unspliced residual (and
optional sample-level residual permutation), as implemented in
{func}`~scatrans.active_score` and used by
{func}`~scatrans.partition_de_by_mechanism`.

This page is **not** a substitute for the manuscript Results (figure numbers
such as “Figure S1” in the paper refer to GR / design checks there, not to the
mode diagram below). Equation numbers here are local to this documentation
page; the residual core matches the manuscript Methods (reference-corrected
excess and Huber gene-structure residualization).

**Primary software path (package 0.10.9):** DE
defines gene-list membership; the residual annotates transcription vs.
stabilization support; program-level calls (competitive and optional
induction-matched) are preferred over hard per-gene claims; for **absolute**
program placement use
{func}`~scatrans.program_mechanism_permutation_calibrated` (observed mean minus
the same gene set’s expectation under shuffled condition labels — Methods
subsection *Permutation-calibrated program placement*); optional
`add_nascent_score` is a **detection** axis and never rewrites mechanism labels.
Gene-wise residual permutation FDR (`use_permutation=True` on `active_score`) is
**optional** (Mode A) and separate from program-level calibration. Product
scope: {doc}`faq`.

- Output columns and reporting: {doc}`statistical_guidance`
- Domain conventions: {doc}`domain_assumptions`
- Function signatures: {doc}`api/index`, {doc}`user_guide/index`
```

## Introduction

Conventional single-cell differential expression (DE) compares mature mRNA
abundance between conditions. Mature levels mix transcription, processing,
degradation, and cellular history, so the same fold-change can arise from
faster synthesis or slower decay.

Many scRNA-seq protocols also report unspliced (intron-containing) counts
alongside spliced mature mRNA. Unspliced molecules are a closer proxy for
recent transcription than mature mRNA, but they still depend on splicing
kinetics, gene structure, and capture efficiency. Under the standard kinetic
model, unspliced abundance $u$ and spliced abundance $s$ of a gene evolve as

$$\frac{du}{dt} = \alpha - \beta u,\quad\quad\frac{ds}{dt} = \beta u - \gamma s,$$

with transcription rate $\alpha$, splicing rate $\beta$, and mRNA degradation
rate $\gamma$; at steady state $u^{*}/s^{*} = \gamma/\beta$.

Two practical points follow. First, the observed steady-state
unspliced-to-spliced ratio identifies only the compound $\gamma/\beta$, not a
single kinetic rate, so a ratio fit from data is an empirical calibration
constant. Second, when conditions differ, unspliced abundance can deviate from
the value predicted by a reference unspliced-to-spliced relationship; that
*excess* is a static, condition-comparative signal for mechanism annotation
among DE genes, not a substitute discovery score.

Among genes DE already selected, the residual helps separate
transcription-weighted from stabilization-weighted changes, especially when
pooled at the program level. Early after a stimulus it can also carry nascent
signal that later decays as mature mRNA accumulates. Limits include intron
capture quality, kinetic regime, and power; product scope is in {doc}`faq`.

Cell-level DE can pseudoreplicate when cells from one biological sample are
treated as independent. Pseudobulk makes the replicate the unit of inference,
but raw unspliced counts remain confounded by gene length and intron number.

**scATrans** pairs a pluggable DE step (membership) with a reference-corrected,
bias-adjusted unspliced residual (annotation). Optional sample-level
permutation calibrates residual significance when biological replicates exist.
The recommended entry point is
{func}`~scatrans.partition_de_by_mechanism`.

## Methods

### Overview and default configuration

**Mechanism partition (recommended software path).**
{func}`~scatrans.partition_de_by_mechanism` (i) runs regime diagnosis, (ii)
selects genes by DE (builtin, external table, or backend kwargs; prefer
pseudobulk when `sample_col` provides replicates), (iii) computes the
reference-corrected unspliced excess and optional Huber gene-structure
residualization below, (iv) forms soft mechanism labels and optional
program-level / induction-matched tables, and (v) optionally appends a
**decoupled** nascent detection score (`add_nascent_score=True`). Gene-wise
residual permutation is **not** required for that path.

**Program absolute placement (recommended for claims that need a zero).**
{func}`~scatrans.program_mechanism_permutation_calibrated` places a gene set on
the mechanism axis against an **empirical** null: the same program’s mean
`transcription_support` under shuffled condition labels. Competitive /
induction-matched program tests are *relative* (background or matched
non-members); genome-wide standardization does not remove a gene set’s own
structural offset (length, intron content). See **Program-level permutation
calibration** below.

**Optional gene-wise residual permutation path (Mode A).** When biological
sample structure exists, the lower-level scorer can recompute the residual under
sample-level label permutation to obtain residual FDR (`unspliced_excess_fdr`).
That is a separate, gene-level calibration — not DE membership and not the
program-level absolute placement API.

For each gene the residual core is: (i) a shrunken reference
unspliced-to-spliced ratio, (ii) target-group unspliced excess, (iii) optional
gene-structure bias residualization. Alternative observational units, ratio
estimators, and DE backends share the same mathematical core (analysis routes /
Modes A–C below; package documentation, not manuscript figure numbers).

### Notation and pseudobulk construction

Cells are aggregated by biological sample, so that the observational unit $i$ is a **biological replicate**, and the spliced and unspliced layers are size-factor normalized after aggregation; samples with fewer than `min_cells` cells or fewer than `min_counts` counts are discarded. Gene $g$ indexes the gene of interest and $h$ is the running index in sums over the gene set $\mathcal{G}$ ($\left| \mathcal{G} \right| = G$). $T$ and $R$ denote the target and reference groups, $U_{i,g}$ and $S_{i,g}$ the unspliced and spliced abundances, and ${\bar{U}}_{T,g}$, ${\bar{S}}_{T,g}$, ${\bar{U}}_{R,g}$, ${\bar{S}}_{R,g}$ the corresponding group means. To avoid collision with the kinetic parameters of Eq. 1, the *empirical* unspliced-to-spliced ratio is written $\rho$ throughout; under steady state $\rho$ estimates $\gamma/\beta$ and must not be read as a kinetic rate.

$\mathcal{G}_{valid}$ denotes the genes retained for modeling: those with available length and intron annotation and with total counts $\sum_{i}^{}\left( U_{i,h} + S_{i,h} \right) \geq \kappa$ (`min_total_counts`, default 50).

### Step 1: Reference-corrected unspliced excess

The gene-wise reference ratio is shrunk toward a global background ratio $\rho_{0}$ by additive pseudo-count shrinkage, and the excess is the vertical deviation of the target group from the reference steady-state line:

$${\widehat{\rho}}_{R,g} = \frac{{\bar{U}}_{R,g} + \eta\,\rho_{0}}{{\bar{S}}_{R,g} + \eta},\quad\quad\rho_{0} = \frac{\sum_{i \in R}^{}{\sum_{h = 1}^{G}U_{i,h}} + \epsilon}{\sum_{i \in R}^{}{\sum_{h = 1}^{G}S_{i,h}} + \epsilon},$$

$$\Delta_{g} = {\bar{U}}_{T,g} - {\widehat{\rho}}_{R,g}\,{\bar{S}}_{T,g}.$$

Here $\epsilon = 10^{- 8}$ is a numerical constant and $\eta > 0$ is a prior weight (default $\eta = 5$). Eq. 2 is equivalent to adding $\eta$ pseudo-spliced units carrying an expected unspliced abundance $\eta\rho_{0}$: it stabilizes the ratio of low-coverage genes while leaving well-covered genes essentially unchanged, since ${\widehat{\rho}}_{R,g} \rightarrow {\bar{U}}_{R,g}/{\bar{S}}_{R,g}$ when ${\bar{S}}_{R,g} \gg \eta$. A positive $\Delta_{g}$ indicates more unspliced RNA in the target group than expected from its mature-mRNA level under the reference $U/S$ relationship. In the implementation this quantity is also stored under a legacy `velocity_delta_raw` alias, but it is a static, condition-comparative contrast and must not be interpreted as a dynamical RNA-velocity estimate.

### Step 2: Robust gene-structure bias correction

Gene length $L_{g}$ and intron number $I_{g}$ systematically influence unspliced-read recovery. With the design vector $\mathbf{x}_{g} = \left\lbrack \, 1,\ \log\left( 1 + L_{g} \right),\ \log\left( 1 + I_{g} \right)\, \right\rbrack^{\top}$, coefficients are fitted once, globally, by weighted Huber M-estimation and subtracted:

$$\widehat{\mathbf{\theta}} = \arg\min_{\mathbf{\theta}}\sum_{h \in \mathcal{G}_{valid}}^{}\omega_{h}\,\rho_{c}\left( \frac{\Delta_{h} - \mathbf{x}_{h}^{\top}\mathbf{\theta}}{\sigma} \right),\quad\quad R_{g} = \Delta_{g} - \mathbf{x}_{g}^{\top}\widehat{\mathbf{\theta}}.$$

$\rho_{c}$ is the Huber loss with threshold $c = 1.35$ (quadratic for $|e| \leq c$, linear beyond), $\sigma$ a jointly estimated scale, and the observation weights $\omega_{h}$ are the total unspliced-plus-spliced counts of each gene, winsorized at their 95th percentile so that a few very highly expressed genes cannot dominate the fit. The residual $R_{g}$ is the **bias-corrected unspliced excess**. Expressed genes lacking length/intron annotation are median-centered instead ($R_{g} = \Delta_{g} - {median}_{h}\Delta_{h}$); genes outside $\mathcal{G}_{valid}$ are assigned $R_{g} = 0$; and if the regression cannot be fitted (fewer than 30 usable genes, or numerical failure) the pipeline falls back to global median centering and records `bias_corrected = False` in the run diagnostics. Because $\Delta_{g}$ is defined on the normalized-count scale while the covariates are log-transformed, $\widehat{\mathbf{\theta}}$ should be read as an empirical bias-removal fit rather than a mechanistic model of capture efficiency.

### Step 3: Differential expression

Pseudobulk counts are modeled with PyDESeq2 (negative-binomial GLM, Wald test, Benjamini–Hochberg adjustment), yielding ${logFC}_{g}$ and the adjusted p-value $p_{g}^{adj}$; the latter is the DE backend’s own adjusted p-value and is distinct from the permutation-calibrated FDR of Step 4. The log-fold change and adjusted significance are reported alongside the bias-corrected residual $R_{g}$; the residual is the quantity carried into the permutation calibration of Step 4 and into the induction-normalized mechanism annotation. Alternative DE backends (Supplementary Methods S3) supply ${logFC}_{g}$ and $p_{g}^{adj}$ without otherwise changing the pipeline.

### Step 4 (optional): Gene-wise residual FDR by sample-level permutation

Condition labels are permuted $B$ times while preserving the target and reference group sizes. Permutation is performed **at the pseudobulk sample level, never at the cell level**, so that the null distribution respects the biological replicate as the unit of inference and pseudoreplication is not reintroduced through the null. For each permutation $b$ the residual pipeline — unspliced excess and bias correction — is recomputed, yielding null residual values $R_{g}^{(b)}$. With $B_{eff} \leq B$ the number of permutations that completed successfully, the one-sided empirical p-values and their BH-adjusted counterparts are

$$p_{g}^{perm,R} = \frac{1 + \sum_{b = 1}^{B_{eff}}\mathbf{1}\left( R_{g}^{(b)} \geq R_{g}^{obs} \right)}{B_{eff} + 1},\quad\quad\{ q_{h}^{perm,R}\}_{h \in \mathcal{G}_{valid}} = BH\left( \{ p_{h}^{perm,R}\}_{h \in \mathcal{G}_{valid}} \right).$$

The $+ 1$ correction keeps the p-values strictly positive and valid. The residual-based value $q_{g}^{perm,R}$ (`unspliced_excess_fdr`) is an optional residual FDR (also used by the lower-level built-in `significant` conjunction on `active_score` when permutation is enabled), with the one-sided requirement $R_{g} > 0$. It does **not** define DE membership on the primary `partition_de_by_mechanism` path. Because the attainable resolution of an empirical p-value is $1/\left( B_{eff} + 1 \right)$, FDR values are reported as usable only when $B_{eff} \geq 100$; otherwise the p-values are returned but flagged (`use_fdr = False`, reason `small_permutation_space`). The permutation null tests the exchangeability of condition labels.

### Program-level permutation calibration (absolute placement)

Matches the manuscript Methods subsection *Permutation-calibrated program
placement* and the package function
{func}`~scatrans.program_mechanism_permutation_calibrated` (0.10.9).

Competitive {func}`~scatrans.program_mechanism` and induction-matched
{func}`~scatrans.program_mechanism_induction_matched` are **relative**: they
rank a set against a background or matched non-members. Neither supplies a
zero — `transcription_support` is standardized genome-wide, while any gene set
has its own length / intron / expression composition.

**Construction.** Let $m$ be the mean `transcription_support` of the program’s
members. Condition labels in `obs[groupby]` are permuted across cells
(optionally **within** `block_col` strata such as donor/batch so that
condition-by-block imbalance is preserved in the null). Each replicate is
pushed through the **identical** chain as the observation: same frozen DE
table, residual estimator, gene-structure bias correction, robust-z, and gene
universe — yielding replicate means $m_b^\*$. Membership is taken from the
observed run and never recomputed (a precomputed `de=` is required; otherwise
the null absorbs DE-selection noise). Then

$$\mathrm{calibrated} = m - \overline{m^\*},\qquad
\mathrm{null\_sd} = \mathrm{sd}_b(m_b^\*).$$

`null_sd` is the SD of the null replicate means — how stable the structural
offset is across shuffles — **not** the sampling SE of the observed mean $m$.
Negative `calibrated` → stabilization-weighted; positive → transcription-
weighted.

**What to report.** The manuscript focuses on the displacement `calibrated`
(and shows that a transcriptional control returns near zero while a true
stabilization program does not). The package still returns several optional
summaries so users can choose what fits their claim:

| Column | Role |
|--------|------|
| `calibrated` | Absolute displacement (recommended for absolute placement claims) |
| `null_mean` / `null_sd` | Structural offset and its shuffle-to-shuffle stability |
| `p_perm` | Phipson–Smyth $(b+1)/(n+1)$ (never zero; floor $1/(n+1)$) |
| `z` | Optional `calibrated / null_sd` (descriptive standardization of the offset) |

**Permutation p-value.** Default implementation: **two-sided** count
$\lvert m_b^\* - \overline{m^\*} \rvert \ge \lvert\mathrm{calibrated}\rvert$
(centered on the null mean). Manuscript prose sometimes phrases extremity in
the observed direction; for programs far in one tail both conventions often
hit the same floor $p=1/(n+1)$ (e.g. $n=50$ → $0.020$). Either framing is
usable if stated clearly.

**Membership.** Frozen `de=` is required so the null does not re-run DE
selection. `restrict_to_selected=True` intersects with the observed DE list
(same universe as program figures); leave it `False` to use all finite-support
members of the supplied gene set.

**Locked GSE226488 illustration (S5, `restrict_to_selected` on the 803-gene
DE list):** ARE ($n=16$) stays strongly stabilization-displaced after
calibration (observed far negative; label-shuffle null is a modest structural
offset; calibrated remains large). Primary NF-κB returns near the empirical
zero after calibration (structural offset removed without erasing the ARE
effect). That pair is the verification the Discussion requires of a correct
offset correction. Report final null / calibrated figures from the production
$n_{\mathrm{perm}}$ run (package default 200; manuscript uses a finer grid).

### Limitations

Two limitations follow directly from the construction. First, residual magnitudes are data-adaptive within a run: use them for mechanism annotation and within-run visualization among DE-selected genes, not as a substitute DE discovery score or as absolute cross-run units. Optional residual permutation FDR (Step 4) calibrates residual magnitude under a conditional null; gene-list membership remains DE-defined. Program-level calibration removes a gene set’s structural offset under label exchangeability but still depends on capture regime and design (prefer biological replicates / `block_col` when available). Second, the replicate-aware guarantees above depend on the availability of biological replicates. When replicates are unavailable, scATrans provides cell-level backends (Supplementary Methods, Mode B), but the resulting p-values are subject to pseudoreplication and permutation at the cell level cannot correct it; such residual calibrations should be regarded as exploratory only.

## Supplementary Methods

### S1. Analysis routes implemented in scATrans

**Mode A** is the recommended replicate-aware residual route when biological
samples exist. The software also implements a cell-level route (**Mode B**) and
a moment-smoothed cell-level route (**Mode C**). All three share the residual
core (reference-corrected excess and Huber bias correction) and differ in
(i) the observational unit $i$ used for group means $\bar{U}$ and $\bar{S}$,
and (ii) the DE backend supplying ${logFC}_{g}$ and $p_{g}^{adj}$. The diagram
below is a **documentation schematic** (analysis routes), not manuscript
Figure S1.

```{figure} _static/method_routes_s1.png
:name: fig-method-routes
:width: 88%
:align: center

**Analysis routes in scATrans (Modes A–C).** Modes differ in the observational unit and the DE backend; the reference-corrected excess and Huber bias correction are shared. Only Mode A permits sample-level residual permutation for residual FDR.
```

#### Mode A — Pseudobulk (main text; recommended)

Observational unit $i$ = biological sample. Cells are aggregated per `sample_col`, layers are size-factor normalized, and DE is estimated with PyDESeq2. Equations: **2, 3, 4, 5** exactly as in the main text, with permutation applied at the sample level. This is the only route for which the permutation null is a valid replicate-level null.

#### Mode B — Cell-level

Observational unit $i$ = single cell; no aggregation. Group means in Eqs. 2–3 are taken over cells, and DE is estimated with a Scanpy cell-level test (`t-test_overestim_var` by default, or `wilcoxon`). Eqs. 2–4 are unchanged. Eq. 5 may still be evaluated, but labels are then permuted at the cell level and the resulting p-values do **not** control the false-positive rate arising from pseudoreplication. Use only when biological replicates are unavailable, and report results as exploratory.

#### Mode C — Moment-smoothed (advanced)

Observational unit $i$ = single cell, after $k$-nearest-neighbor smoothing in principal-component space. Raw layers are replaced by first-order moments,

$$M_{i,g}^{U} = \sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j}\, U_{j,g},\quad\quad M_{i,g}^{S} = \sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j}\, S_{j,g},\quad\quad\sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j} = 1,$$

where $\mathcal{N}_{k}(i)$ is the $k$-NN set of cell $i$ (including $i$). Eqs. 2 and 3 are then applied verbatim with $U \rightarrow M^{U}$ and $S \rightarrow M^{S}$:

$${\widehat{\rho}}_{R,g}^{M} = \frac{{\bar{M}}_{R,g}^{U} + \eta\,\rho_{0}^{M}}{{\bar{M}}_{R,g}^{S} + \eta},\quad\quad\rho_{0}^{M} = \frac{\sum_{i \in R}^{}{\sum_{h = 1}^{G}M_{i,h}^{U}} + \epsilon}{\sum_{i \in R}^{}{\sum_{h = 1}^{G}M_{i,h}^{S}} + \epsilon},$$

$$\Delta_{g}^{M} = {\bar{M}}_{T,g}^{U} - {\widehat{\rho}}_{R,g}^{M}\,{\bar{M}}_{T,g}^{S},\quad\quad{\bar{M}}_{T,g}^{U} = \frac{1}{|T|}\sum_{i \in T}^{}M_{i,g}^{U}\ \ \left( \text{and analogously for }{\bar{M}}^{S},\ R \right).$$

Eq. 4 then follows unchanged with $\Delta_{g}^{M}$ in place of $\Delta_{g}$. **Statistical caveat:** neighborhood averaging induces dependence between cells, so Mode C is intended for exploratory, cell-level visualization of nascent-RNA excess and is not combined with permutation inference by default. It should not be used to make replicate-level claims.

### S2. Optional reference-ratio estimators

Eq. 2 (`gamma_method="heuristic_shrink"`) is the default. Three alternatives replace Eq. 2 only; Eqs. 3–5 are unaffected.

**(a)** `robust_median`**.** The shrinkage anchor $\rho_{0}$ is replaced by the median of the per-gene reference ratios, computed over expressed genes only:

$$\rho_{0}^{med} = {median}_{\, h\,:\,{\bar{U}}_{R,h} + {\bar{S}}_{R,h} > 0}\left( \frac{{\bar{U}}_{R,h} + \epsilon}{{\bar{S}}_{R,h} + \epsilon} \right),\quad\quad{\widehat{\rho}}_{R,g} = \frac{{\bar{U}}_{R,g} + \eta\,\rho_{0}^{med}}{{\bar{S}}_{R,g} + \eta}.$$

This is a robust heuristic, not a Bayesian estimator. Zero-expression genes are excluded from the anchor because their ratio $\epsilon/\epsilon \approx 1$ would otherwise dominate on sparse data.

**(b)** `empirical_bayes`**.** Per-gene log-ratios $r_{g} = \log\left( \left( {\bar{U}}_{R,g} + \epsilon \right)/\left( {\bar{S}}_{R,g} + \epsilon \right) \right)$ are shrunk toward a robust prior $\left( \mu_{0},\tau^{2} \right)$ estimated by trimmed median and MAD across genes:

$${\widehat{r}}_{g} = w_{g}r_{g} + \left( 1 - w_{g} \right)\,\mu_{0},\quad w_{g} = \frac{\tau^{2}}{\tau^{2} + \sigma_{g}^{2}},\quad\sigma_{g}^{2} = \frac{1}{n_{R}{\bar{U}}_{R,g} + c} + \frac{1}{n_{R}{\bar{S}}_{R,g} + c},\quad{\widehat{\rho}}_{R,g} = e^{{\widehat{r}}_{g}},$$

where $\sigma_{g}^{2}$ is the delta-method sampling variance of the mean log-ratio, $n_{R}$ the number of reference units and $c$ a pseudo-count derived from `prior_weight`. The posterior log-scale $\sqrt{\tau^{2}\sigma_{g}^{2}/\left( \tau^{2} + \sigma_{g}^{2} \right)}$ is exported as a per-gene diagnostic. During permutation the prior hyper-parameters are held **fixed** at their observed-data values, so that the null is not re-tuned on shuffled labels. Recommended when the reference group is small.

**(c)** `raw`**.** No shrinkage: ${\widehat{\rho}}_{R,g} = \left( {\bar{U}}_{R,g} + \epsilon \right)/\left( {\bar{S}}_{R,g} + \epsilon \right)$ for expressed genes, with $\rho_{0}$ substituted for genes with no reference expression.

### S3. Optional DE backends

The DE backend supplies ${logFC}_{g}$ and $p_{g}^{adj}$ as reported outputs; the rest of the pipeline is unchanged. In addition to PyDESeq2 (Mode A default) and the Scanpy cell-level tests (Mode B default), scATrans implements a mixed linear model (fixed condition effect, random sample effect; Wald or LRT p-value) and Memento (method-of-moments cell-level testing). The residual permutation (Eq. 5) does not depend on the DE backend.

### S4. Bias-correction options

`bias_correction="huber_length_intron"` (default) applies Eq. 4. `bias_correction="none"` sets $R_{g} = \Delta_{g}$ for expressed genes, bypassing the regression entirely; Eq. 5 is otherwise unchanged. The median-centering fallback described in the main text is applied automatically whenever the Huber fit cannot be performed.

### S5. Summary of options

| Option                      | Values (default in bold)                                        | Affects                              | Notes                                           |
|:----------------------------|:----------------------------------------------------------------|:-------------------------------------|:------------------------------------------------|
| Observational unit          | **pseudobulk (Mode A)** / cell (B) / kNN-smoothed cell (C)      | $\bar{U},\bar{S}$ in Eqs. 2–3        | Only A supports replicate-level permutation     |
| DE backend                  | **PyDESeq2** / Scanpy t-test / Wilcoxon / mixed LM / Memento    | ${logFC}_{g},\ p_{g}^{adj}$ in Eq. 5 | Must be identical in observed and permuted runs |
| `gamma_method`              | **heuristic\_shrink** / robust\_median / empirical\_bayes / raw | Eq. 2 → S1, S2                       | EB recommended for small reference groups       |
| `prior_weight` $\eta$       | **5.0**                                                         | Eq. 2                                | Larger = stronger shrinkage to $\rho_{0}$       |
| `bias_correction`           | **huber\_length\_intron** / none                                | Eq. 4                                | Median-centering fallback is automatic           |
| `n_perm` $B$ (residual FDR) | **100**                                                         | Eq. 5 residual p / FDR               | Residual permutation in `active_score` / design; FDR usable only when $B_{eff} \geq 100$. Separate from program calibration (`program_mechanism_permutation_calibrated` defaults to **200**) |
| `min_total_counts` $\kappa$ | **50**                                                          | $\mathcal{G}_{valid}$                | Genes below threshold get $R_{g} = 0$           |

### S6. Equation index (for the software documentation)

| Eq.   | Quantity                                       | Field in output                                                                            |
|:------|:-----------------------------------------------|:-------------------------------------------------------------------------------------------|
| 1     | Kinetic model (context only)                   | —                                                                                          |
| 2     | Reference $U/S$ ratio ${\widehat{\rho}}_{R,g}$ | `effective_gamma` (opt-in)                                                                 |
| 3     | Raw unspliced excess $\Delta_{g}$              | `unspliced_excess_delta` (legacy: `velocity_delta_raw`)                                        |
| 4     | Huber bias fit and residual $R_{g}$            | `unspliced_excess_residual`; coefficients in `uns["scatrans"]["diagnostics"]`              |
| 5     | Permutation p-value and FDR                    | `unspliced_excess_pval` / `unspliced_excess_fdr`                                           |
| C1–C3 | Moment smoothing (Mode C)                      | layers `Mu`, `Ms`                                                                          |
| S1–S2 | Alternative ratio estimators                   | `gamma_info` diagnostics                                                                   |
