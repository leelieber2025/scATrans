# Method: Active Transcription Framework

```{admonition} Scope of this page
:class: note

Mathematical definition of the scATrans nascent residual computed by
{func}`~scatrans.active_score`: the default replicate-aware configuration,
optional routes and estimators, and the induction-normalized unspliced residual
used for mechanism annotation (adapted from the scATrans manuscript).

The recommended software entry point is
{func}`~scatrans.partition_de_by_mechanism`. Product scope and deprecated
paths: {doc}`faq`.

- Output columns and reporting: {doc}`statistical_guidance`
- Domain conventions: {doc}`domain_assumptions`
- Function signatures: {doc}`api/index`, {doc}`user_guide/index`
```

## Introduction

Single-cell RNA sequencing (scRNA-seq) resolves transcriptional states at high
resolution and links molecular changes to phenotype in development, disease, and
perturbation experiments. Conventional single-cell differential expression (DE)
quantifies differences in mature mRNA abundance. Mature mRNA levels reflect the
combined effects of transcription, processing, degradation, and cellular
history, and therefore may not fully capture ongoing transcriptional activation.

Many scRNA-seq protocols assign reads to two related RNA populations within the
same library: spliced mature mRNA and unspliced, intron-containing pre-mRNA.
Unspliced molecules largely represent recently transcribed species and are
therefore a closer proxy for recent transcriptional activity than mature mRNA,
although their abundance is also shaped by splicing kinetics, gene structure,
and technical capture efficiency. Under the standard kinetic model, the
unspliced abundance $u$ and spliced abundance $s$ of a gene evolve as

$$\frac{du}{dt} = \alpha - \beta u,\quad\quad\frac{ds}{dt} = \beta u - \gamma s,$$

with transcription rate $\alpha$, splicing rate $\beta$, and mRNA degradation
rate $\gamma$; at steady state $u^{*}/s^{*} = \gamma/\beta$.

Two consequences motivate scATrans. First, the *observable* steady-state
unspliced-to-spliced ratio identifies only the compound quantity
$\gamma/\beta$, not any individual kinetic rate; a ratio estimated from data
must be treated as an empirical calibration constant, not as a transcription or
degradation rate. Second, away from steady state—for example during cell-state
transitions, stimulus responses, or disease onset—the unspliced abundance can
deviate from the value predicted by the reference ratio, and this *excess* is an
early signal of transcriptional activation before total mRNA has fully changed.

Relative to steady-state DE alone, the nascent-versus-mature contrast can
(i) highlight early transcriptional responses, (ii) help separate transcriptional
from post-transcriptional regulation, and (iii) flag inducible genes before
mature mRNA abundance has fully changed. These properties are limited by capture
quality, kinetic regime, and power; product scope is summarized in {doc}`faq`.

Single-cell DE is also subject to pseudoreplication when cells from one
biological replicate are treated as independent observations. Pseudobulk
methods improve type-I error control by making the replicate the unit of
inference, yet raw unspliced counts remain confounded by gene length and intron
number.

**scATrans** combines replicate-aware DE evidence with a reference-corrected,
bias-adjusted unspliced residual; significance for the residual is obtained by
sample-level permutation when requested. In software, the recommended entry
point is mechanism partition of DE-selected genes
({func}`~scatrans.partition_de_by_mechanism`).

## Methods

### Overview and default configuration

All results reported here were generated with the recommended, replicate-aware configuration of scATrans, which is described in full below: cells are aggregated into pseudobulk profiles per biological sample, differential expression is estimated with PyDESeq2, the unspliced excess is calibrated against a shrunken reference unspliced-to-spliced ratio, gene-structure bias is removed by weighted Huber regression, and significance is obtained by permuting condition labels **at the sample level**. Alternative observational units, ratio estimators and DE backends are implemented in the software; they share the same mathematical core and are specified in the Supplementary Methods (Supplementary Figure S1, Supplementary Table S1).

For each gene, the pipeline (i) estimates a shrunken reference unspliced-to-spliced ratio, (ii) uses it to compute a target-group unspliced excess, (iii) regresses out gene-structure-associated technical bias to obtain the residual, and (iv) calibrates the residual against an empirical null obtained by sample-level label permutation.

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

### Step 4: Empirical calibration by sample-level permutation

Condition labels are permuted $B$ times while preserving the target and reference group sizes. Permutation is performed **at the pseudobulk sample level, never at the cell level**, so that the null distribution respects the biological replicate as the unit of inference and pseudoreplication is not reintroduced through the null. For each permutation $b$ the residual pipeline — unspliced excess and bias correction — is recomputed, yielding null residual values $R_{g}^{(b)}$. With $B_{eff} \leq B$ the number of permutations that completed successfully, the one-sided empirical p-values and their BH-adjusted counterparts are

$$p_{g}^{perm,R} = \frac{1 + \sum_{b = 1}^{B_{eff}}\mathbf{1}\left( R_{g}^{(b)} \geq R_{g}^{obs} \right)}{B_{eff} + 1},\quad\quad\{ q_{h}^{perm,R}\}_{h \in \mathcal{G}_{valid}} = BH\left( \{ p_{h}^{perm,R}\}_{h \in \mathcal{G}_{valid}} \right).$$

The $+ 1$ correction keeps the p-values strictly positive and valid. The residual-based value $q_{g}^{perm,R}$ (`unspliced_excess_fdr`) is the quantity used by the built-in significant-gene list, in conjunction with the one-sided requirement $R_{g} > 0$. Because the attainable resolution of an empirical p-value is $1/\left( B_{eff} + 1 \right)$, FDR values are reported as usable only when $B_{eff} \geq 100$; otherwise the p-values are returned but flagged (`use_fdr = False`, reason `small_permutation_space`). The permutation null tests the exchangeability of condition labels.

### Limitations

Two limitations follow directly from the construction. First, the residual scale parameters are data-adaptive, so residual magnitudes support *prioritization within an analysis* while formal significance is taken from the permutation FDR (Eq. 5) rather than from a residual threshold. Second, the replicate-aware guarantees above depend on the availability of biological replicates. When replicates are unavailable, scATrans provides cell-level backends (Supplementary Methods, Mode B), but the resulting p-values are subject to pseudoreplication and permutation at the cell level cannot correct it; such analyses should be regarded as exploratory prioritization only.

## Supplementary Methods

### S1. Analysis routes implemented in scATrans

The main text describes **Mode A**, the recommended replicate-aware route. The software additionally implements a cell-level route (**Mode B**) and a moment-smoothed cell-level route (**Mode C**). All three share the same mathematical core — Eqs. 2–5 of the main text — and differ only in (i) the definition of the observational unit $i$ used to form the group means $\bar{U}$ and $\bar{S}$, and (ii) the DE backend supplying ${logFC}_{g}$ and $p_{g}^{adj}$. Supplementary Figure S1 makes this shared structure explicit.

```{figure} _static/method_routes_s1.png
:name: fig-method-routes
:width: 88%
:align: center

**Supplementary Figure S1. Analysis routes in scATrans.** Modes A, B and C differ only in the observational unit and the DE backend; the reference-corrected excess (Eqs. 2–3) and the Huber bias correction (Eq. 4) are identical across routes. Only Mode A permits sample-level permutation and therefore replicate-aware significance (Eq. 5).
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
| `n_perm` $B$                | **100**                                                         | Eq. 5                                | FDR usable only when $B_{eff} \geq 100$         |
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
