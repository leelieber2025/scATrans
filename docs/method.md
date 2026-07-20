# Method: the Active Transcription framework

```{admonition} What this page is
:class: note

This page documents the **method** behind {func}`~scatrans.active_score` — the mathematical framework, its default replicate-aware configuration, and the optional routes and estimators — adapted from the scATrans manuscript. It is the reference for *how* the Active Score is defined and calibrated.

- For what each **output column** means and how to report it, see {doc}`statistical_guidance`.
- For domain conventions (upregulation-oriented scoring, residual vs. DE, cutoffs), see {doc}`domain_assumptions`.
- For the **functions and arguments** that implement each equation, see {doc}`api/index` and {doc}`user_guide/index`.
```

## Introduction

Single-cell RNA sequencing (scRNA-seq) has transformed the study of cellular heterogeneity, allowing transcriptional states to be resolved at high resolution and molecular changes to be linked to phenotypic outcomes in development, disease and perturbation experiments. Conventional single-cell differential expression (DE) analysis quantifies differences in mature mRNA abundance. Mature mRNA levels, however, reflect the combined effects of transcription, processing, degradation and cellular history, and therefore may not fully capture ongoing transcriptional activation.

Many modern scRNA-seq protocols assign reads to two related RNA populations within the same library: spliced mature mRNA and unspliced, intron-containing pre-mRNA. Because unspliced molecules largely represent recently transcribed species, they are a closer proxy for recent transcriptional activity than mature mRNA, although their abundance is also shaped by splicing kinetics, gene structure and technical capture efficiency. Under the standard kinetic model, the unspliced abundance $u$ and spliced abundance $s$ of a gene evolve as

$$\frac{du}{dt} = \alpha - \beta u,\quad\quad\frac{ds}{dt} = \beta u - \gamma s,$$

with transcription rate $\alpha$, splicing rate $\beta$ and mRNA degradation rate $\gamma$; at steady state $u^{*}/s^{*} = \gamma/\beta$.

Two consequences motivate scATrans. First, the *observable* steady-state unspliced-to-spliced ratio identifies only the compound quantity $\gamma/\beta$, not any individual kinetic rate; a ratio estimated from data must therefore be treated as an empirical calibration constant, not as a transcription or degradation rate. Second, away from steady state — precisely the regime of cell-state transitions, stimulus responses and disease onset — the unspliced abundance deviates from the value predicted by the reference ratio, and this *excess* provides an early signal of transcriptional activation before changes in total mRNA become detectable.

This nascent-versus-mature distinction offers several mechanistic advantages over steady-state DE: (i) it prioritizes early candidate regulators associated with phenotypic reprogramming; (ii) it helps separate transcriptional regulation from post-transcriptional buffering; (iii) it can reveal poised or rapidly inducible genes whose mature mRNA does not yet show differential abundance; and (iv) it carries information about near-future transcriptional states on short timescales.

Rigorously quantifying active transcription nevertheless presents analytical challenges. A well-documented limitation of conventional single-cell DE is pseudoreplication: treating thousands of cells from one biological replicate as independent observations inflates statistical power and generates false-positive DE genes. Pseudobulk approaches substantially improve type-I error control by making the biological replicate the unit of inference, yet raw unspliced counts remain confounded by gene length and intron number. Few existing pipelines simultaneously address pseudoreplication, transcriptional dynamics, technical bias correction and empirical, score-level calibration.

We therefore developed **scATrans**, a single-cell **A**ctive **Trans**cription analysis framework that exploits the nascent-versus-mature distinction to prioritize genes with active transcriptional signatures. scATrans combines replicate-aware DE evidence with a reference-corrected, bias-adjusted unspliced residual, and summarizes these signals in a permutation-calibrated **Composite Active Evidence Score** bounded between 0 and 100. The score is explicitly an interpretable heuristic composite rather than a likelihood-based test statistic; statistical significance is supplied separately, by permutation calibration of the bias-corrected unspliced excess at the level of biological replicates.

## Methods

### Overview and default configuration

All results reported here were generated with the recommended, replicate-aware configuration of scATrans, which is described in full below: cells are aggregated into pseudobulk profiles per biological sample, differential expression is estimated with PyDESeq2, the unspliced excess is calibrated against a shrunken reference unspliced-to-spliced ratio, gene-structure bias is removed by weighted Huber regression, and significance is obtained by permuting condition labels **at the sample level**. Alternative observational units, ratio estimators and DE backends are implemented in the software; they share the same mathematical core and are specified in the Supplementary Methods (Supplementary Figure S1, Supplementary Table S1).

For each gene, the pipeline (i) estimates a shrunken reference unspliced-to-spliced ratio, (ii) uses it to compute a target-group unspliced excess, (iii) regresses out gene-structure-associated technical bias, (iv) combines the resulting residual with log-fold change and adjusted significance into a bounded Active Score, and (v) calibrates the residual against an empirical null obtained by sample-level label permutation.

### Notation and pseudobulk construction

Cells are aggregated by biological sample, so that the observational unit $i$ is a **biological replicate**, and the spliced and unspliced layers are size-factor normalized after aggregation; samples with fewer than `min_cells` cells or fewer than `min_counts` counts are discarded. Gene $g$ indexes the gene of interest and $h$ is the running index in sums over the gene set $\mathcal{G}$ ($\left| \mathcal{G} \right| = G$). $T$ and $R$ denote the target and reference groups, $U_{i,g}$ and $S_{i,g}$ the unspliced and spliced abundances, and ${\bar{U}}_{T,g}$, ${\bar{S}}_{T,g}$, ${\bar{U}}_{R,g}$, ${\bar{S}}_{R,g}$ the corresponding group means. To avoid collision with the kinetic parameters of Eq. 1, the *empirical* unspliced-to-spliced ratio is written $\rho$ throughout; under steady state $\rho$ estimates $\gamma/\beta$ and must not be read as a kinetic rate.

$\mathcal{G}_{valid}$ denotes the genes retained for modeling: those with available length and intron annotation and with total counts $\sum_{i}^{}\left( U_{i,h} + S_{i,h} \right) \geq \kappa$ (`min_total_counts`, default 50).

### Step 1: Reference-corrected unspliced excess

The gene-wise reference ratio is shrunk towards a global background ratio $\rho_{0}$ by additive pseudo-count shrinkage, and the excess is the vertical deviation of the target group from the reference steady-state line:

$${\widehat{\rho}}_{R,g} = \frac{{\bar{U}}_{R,g} + \eta\,\rho_{0}}{{\bar{S}}_{R,g} + \eta},\quad\quad\rho_{0} = \frac{\sum_{i \in R}^{}{\sum_{h = 1}^{G}U_{i,h}} + \epsilon}{\sum_{i \in R}^{}{\sum_{h = 1}^{G}S_{i,h}} + \epsilon},$$

$$\Delta_{g} = {\bar{U}}_{T,g} - {\widehat{\rho}}_{R,g}\,{\bar{S}}_{T,g}.$$

Here $\epsilon = 10^{- 8}$ is a numerical constant and $\eta > 0$ is a prior weight (default $\eta = 5$). Eq. 2 is equivalent to adding $\eta$ pseudo-spliced units carrying an expected unspliced abundance $\eta\rho_{0}$: it stabilises the ratio of low-coverage genes while leaving well-covered genes essentially unchanged, since ${\widehat{\rho}}_{R,g} \rightarrow {\bar{U}}_{R,g}/{\bar{S}}_{R,g}$ when ${\bar{S}}_{R,g} \gg \eta$. A positive $\Delta_{g}$ indicates more unspliced RNA in the target group than expected from its mature-mRNA level under the reference $U/S$ relationship. In the implementation this quantity is also stored under a legacy `velocity_delta_raw` alias, but it is a static, condition-comparative contrast and must not be interpreted as a dynamical RNA-velocity estimate.

### Step 2: Robust gene-structure bias correction

Gene length $L_{g}$ and intron number $I_{g}$ systematically influence unspliced-read recovery. With the design vector $\mathbf{x}_{g} = \left\lbrack \, 1,\ \log\left( 1 + L_{g} \right),\ \log\left( 1 + I_{g} \right)\, \right\rbrack^{\top}$, coefficients are fitted once, globally, by weighted Huber M-estimation and subtracted:

$$\widehat{\mathbf{\theta}} = \arg\min_{\mathbf{\theta}}\sum_{h \in \mathcal{G}_{valid}}^{}\omega_{h}\,\rho_{c}\left( \frac{\Delta_{h} - \mathbf{x}_{h}^{\top}\mathbf{\theta}}{\sigma} \right),\quad\quad R_{g} = \Delta_{g} - \mathbf{x}_{g}^{\top}\widehat{\mathbf{\theta}}.$$

$\rho_{c}$ is the Huber loss with threshold $c = 1.35$ (quadratic for $|e| \leq c$, linear beyond), $\sigma$ a jointly estimated scale, and the observation weights $\omega_{h}$ are the total unspliced-plus-spliced counts of each gene, winsorized at their 95th percentile so that a few very highly expressed genes cannot dominate the fit. The residual $R_{g}$ is the **bias-corrected unspliced excess**. Expressed genes lacking length/intron annotation are median-centered instead ($R_{g} = \Delta_{g} - {median}_{h}\Delta_{h}$); genes outside $\mathcal{G}_{valid}$ are assigned $R_{g} = 0$; and if the regression cannot be fitted (fewer than 30 usable genes, or numerical failure) the pipeline falls back to global median centering and records `bias_corrected = False` in the run diagnostics. Because $\Delta_{g}$ is defined on the normalized-count scale while the covariates are log-transformed, $\widehat{\mathbf{\theta}}$ should be read as an empirical bias-removal fit rather than a mechanistic model of capture efficiency.

### Step 3: Composite Active Evidence Score

Pseudobulk counts are modelled with PyDESeq2 (negative-binomial GLM, Wald test, Benjamini–Hochberg adjustment), yielding ${logFC}_{g}$ and the adjusted p-value $p_{g}^{adj}$; the latter is the DE backend’s own adjusted p-value and is distinct from the permutation-calibrated FDR of Step 4. The significance evidence is $E_{g} = - \log_{10}\left( p_{g}^{adj} + 10^{- 300} \right)$.

The three raw components $z_{g}^{FC} = {logFC}_{g}$, $z_{g}^{V} = R_{g}$ and $z_{g}^{P} = E_{g}$ live on different, unbounded scales. Each is mapped to $\lbrack 0,1\rbrack$ by a monotone, one-sided saturating transform and the three legs are combined as a weighted mean, rescaled to $\lbrack 0,100\rbrack$:

$$s_{g}^{k} = \delta_{g}^{k}\left\lbrack \, 1 - \exp\left( - \frac{\max\left( z_{g}^{k},\, 0 \right)}{\lambda_{k}} \right) \right\rbrack,\quad\quad A_{g} = 100 \times \frac{\sum_{k}^{}w_{k}\, s_{g}^{k}}{\sum_{k}^{}w_{k}},\quad\quad k \in \{ FC,\, V,\, P\}.$$

Three elements of Eq. 5 require comment.

**One-sidedness.** The $\max( \cdot ,0)$ truncation makes every leg one-sided, so that only evidence for *upregulation* and for *positive* unspliced excess can contribute to the score.

**Direction gating.** $E_{g}$ is directionless — a strongly *down*regulated gene can also have a very small $p^{adj}$. The significance leg is therefore gated by the sign of the effect estimate $d_{g}$ that $p_{g}^{adj}$ actually tests (the PyDESeq2 log-fold change): $\delta_{g}^{P} = \mathbf{1}\left( d_{g} > 0 \right)$, while $\delta_{g}^{FC} = \delta_{g}^{V} = 1$.

**Adaptive scale.** Each $\lambda_{k}$ is estimated from the empirical distribution of its own component as $\lambda_{k} = median\{ z_{h}^{k}:z_{h}^{k} > 0\}/\ln 2$ (subject to a lower bound for numerical stability), so that a gene at the median of the positive values receives $s^{k} = 0.5$; for the significance leg, $\lambda_{P}$ is estimated on direction-positive genes only. Because $\lambda_{k}$ depends on which genes enter the analysis, $A_{g}$ is a **within-analysis relative quantity** and is not directly comparable across datasets or gene filters.

The weights $w_{FC},w_{V},w_{P}$ are non-negative and equal to 1 by default. Larger $A_{g}$ indicates stronger combined evidence for increased mature-mRNA abundance, excess unspliced RNA relative to the reference state, and DE significance. $A_{g}$ is an explicitly heuristic composite chosen for interpretability rather than optimality; it is not a likelihood-based test statistic, and it is the identical treatment of observed and permuted data in Step 4 that licenses its empirical calibration.

### Step 4: Empirical calibration by sample-level permutation

Condition labels are permuted $B$ times while preserving the target and reference group sizes. Permutation is performed **at the pseudobulk sample level, never at the cell level**, so that the null distribution respects the biological replicate as the unit of inference and pseudoreplication is not reintroduced through the null. For each permutation $b$ the entire pipeline — DE statistics, unspliced excess, bias correction, soft scaling with the same $\lambda_{k}$ rule, and score assembly — is recomputed, yielding null values $X_{g}^{(b)}$ for $X \in \{ R,A\}$. With $B_{eff} \leq B$ the number of permutations that completed successfully, the one-sided empirical p-values and their BH-adjusted counterparts are

$$p_{g}^{perm,X} = \frac{1 + \sum_{b = 1}^{B_{eff}}\mathbf{1}\left( X_{g}^{(b)} \geq X_{g}^{obs} \right)}{B_{eff} + 1},\quad\quad\{ q_{h}^{perm,X}\}_{h \in \mathcal{G}_{valid}} = BH\left( \{ p_{h}^{perm,X}\}_{h \in \mathcal{G}_{valid}} \right).$$

The $+ 1$ correction keeps the p-values strictly positive and valid. The residual-based value $q_{g}^{perm,V}$ (`unspliced_excess_fdr`) is the quantity used by the built-in significant-gene list, in conjunction with the one-sided requirement $R_{g} > 0$; the score-based $q_{g}^{perm,A}$ is reported for reference. Because the attainable resolution of an empirical p-value is $1/\left( B_{eff} + 1 \right)$, FDR values are reported as usable only when $B_{eff} \geq 100$; otherwise the p-values are returned but flagged (`use_fdr = False`, reason `small_permutation_space`). The permutation null tests the exchangeability of condition labels; it does not, by itself, validate the composite weighting scheme of Eq. 5.

### Limitations

Two limitations follow directly from the construction. First, $A_{g}$ is a heuristic composite whose scale parameters are data-adaptive; it should be used for *prioritization within an analysis*, and formal significance should be taken from Eq. 6 rather than from a score threshold. Second, the replicate-aware guarantees above depend on the availability of biological replicates. When replicates are unavailable, scATrans provides cell-level backends (Supplementary Methods, Mode B), but the resulting p-values are subject to pseudoreplication and permutation at the cell level cannot correct it; such analyses should be regarded as exploratory prioritization only.

## Supplementary Methods

### S1. Analysis routes implemented in scATrans

The main text describes **Mode A**, the recommended replicate-aware route. The software additionally implements a cell-level route (**Mode B**) and a moment-smoothed cell-level route (**Mode C**). All three share the same mathematical core — Eqs. 2–6 of the main text — and differ only in (i) the definition of the observational unit $i$ used to form the group means $\bar{U}$ and $\bar{S}$, and (ii) the DE backend supplying ${logFC}_{g}$ and $p_{g}^{adj}$. Supplementary Figure S1 makes this shared structure explicit.

```{figure} _static/method_routes_s1.png
:name: fig-method-routes
:width: 88%
:align: center

**Supplementary Figure S1. Analysis routes in scATrans.** Modes A, B and C differ only in the observational unit and the DE backend; the reference-corrected excess (Eqs. 2–3), the Huber bias correction (Eq. 4) and the composite score (Eq. 5) are identical across routes. Only Mode A permits sample-level permutation and therefore replicate-aware significance (Eq. 6).
```

#### Mode A — Pseudobulk (main text; recommended)

Observational unit $i$ = biological sample. Cells are aggregated per `sample_col`, layers are size-factor normalized, and DE is estimated with PyDESeq2. Equations: **2, 3, 4, 5, 6** exactly as in the main text, with permutation applied at the sample level. This is the only route for which the permutation null is a valid replicate-level null.

#### Mode B — Cell-level

Observational unit $i$ = single cell; no aggregation. Group means in Eqs. 2–3 are taken over cells, and DE is estimated with a Scanpy cell-level test (`t-test_overestim_var` by default, or `wilcoxon`). Eqs. 2–5 are unchanged. Eq. 6 may still be evaluated, but labels are then permuted at the cell level and the resulting p-values do **not** control the false-positive rate arising from pseudoreplication. Use only when biological replicates are unavailable, and report results as exploratory.

#### Mode C — Moment-smoothed (advanced)

Observational unit $i$ = single cell, after $k$-nearest-neighbor smoothing in principal-component space. Raw layers are replaced by first-order moments,

$$M_{i,g}^{U} = \sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j}\, U_{j,g},\quad\quad M_{i,g}^{S} = \sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j}\, S_{j,g},\quad\quad\sum_{j \in \mathcal{N}_{k}(i)}^{}a_{i,j} = 1,$$

where $\mathcal{N}_{k}(i)$ is the $k$-NN set of cell $i$ (including $i$). Eqs. 2 and 3 are then applied verbatim with $U \rightarrow M^{U}$ and $S \rightarrow M^{S}$:

$${\widehat{\rho}}_{R,g}^{M} = \frac{{\bar{M}}_{R,g}^{U} + \eta\,\rho_{0}^{M}}{{\bar{M}}_{R,g}^{S} + \eta},\quad\quad\rho_{0}^{M} = \frac{\sum_{i \in R}^{}{\sum_{h = 1}^{G}M_{i,h}^{U}} + \epsilon}{\sum_{i \in R}^{}{\sum_{h = 1}^{G}M_{i,h}^{S}} + \epsilon},$$

$$\Delta_{g}^{M} = {\bar{M}}_{T,g}^{U} - {\widehat{\rho}}_{R,g}^{M}\,{\bar{M}}_{T,g}^{S},\quad\quad{\bar{M}}_{T,g}^{U} = \frac{1}{|T|}\sum_{i \in T}^{}M_{i,g}^{U}\ \ \left( \text{and analogously for }{\bar{M}}^{S},\ R \right).$$

Eqs. 4 and 5 then follow unchanged with $\Delta_{g}^{M}$ in place of $\Delta_{g}$. **Statistical caveat:** neighborhood averaging induces dependence between cells, so Mode C is intended for exploratory, cell-level visualisation of nascent-RNA excess and is not combined with permutation inference by default. It should not be used to make replicate-level claims.

### S2. Optional reference-ratio estimators

Eq. 2 (`gamma_method="heuristic_shrink"`) is the default. Three alternatives replace Eq. 2 only; Eqs. 3–6 are unaffected.

**(a)** `robust_median`**.** The shrinkage anchor $\rho_{0}$ is replaced by the median of the per-gene reference ratios, computed over expressed genes only:

$$\rho_{0}^{med} = {median}_{\, h\,:\,{\bar{U}}_{R,h} + {\bar{S}}_{R,h} > 0}\left( \frac{{\bar{U}}_{R,h} + \epsilon}{{\bar{S}}_{R,h} + \epsilon} \right),\quad\quad{\widehat{\rho}}_{R,g} = \frac{{\bar{U}}_{R,g} + \eta\,\rho_{0}^{med}}{{\bar{S}}_{R,g} + \eta}.$$

This is a robust heuristic, not a Bayesian estimator. Zero-expression genes are excluded from the anchor because their ratio $\epsilon/\epsilon \approx 1$ would otherwise dominate on sparse data.

**(b)** `empirical_bayes`**.** Per-gene log-ratios $r_{g} = \log\left( \left( {\bar{U}}_{R,g} + \epsilon \right)/\left( {\bar{S}}_{R,g} + \epsilon \right) \right)$ are shrunk towards a robust prior $\left( \mu_{0},\tau^{2} \right)$ estimated by trimmed median and MAD across genes:

$${\widehat{r}}_{g} = w_{g}r_{g} + \left( 1 - w_{g} \right)\,\mu_{0},\quad w_{g} = \frac{\tau^{2}}{\tau^{2} + \sigma_{g}^{2}},\quad\sigma_{g}^{2} = \frac{1}{n_{R}{\bar{U}}_{R,g} + c} + \frac{1}{n_{R}{\bar{S}}_{R,g} + c},\quad{\widehat{\rho}}_{R,g} = e^{{\widehat{r}}_{g}},$$

where $\sigma_{g}^{2}$ is the delta-method sampling variance of the mean log-ratio, $n_{R}$ the number of reference units and $c$ a pseudo-count derived from `prior_weight`. The posterior log-scale $\sqrt{\tau^{2}\sigma_{g}^{2}/\left( \tau^{2} + \sigma_{g}^{2} \right)}$ is exported as a per-gene diagnostic. During permutation the prior hyper-parameters are held **fixed** at their observed-data values, so that the null is not re-tuned on shuffled labels. Recommended when the reference group is small.

**(c)** `raw`**.** No shrinkage: ${\widehat{\rho}}_{R,g} = \left( {\bar{U}}_{R,g} + \epsilon \right)/\left( {\bar{S}}_{R,g} + \epsilon \right)$ for expressed genes, with $\rho_{0}$ substituted for genes with no reference expression.

### S3. Optional DE backends

The DE backend supplies only ${logFC}_{g}$ and $p_{g}^{adj}$ to Eq. 5; the rest of the pipeline is unchanged. In addition to PyDESeq2 (Mode A default) and the Scanpy cell-level tests (Mode B default), scATrans implements a mixed linear model (fixed condition effect, random sample effect; Wald or LRT p-value, with the fixed-effect coefficient used as the direction-gating effect $d_{g}$ in Eq. 5) and Memento (method-of-moments cell-level testing). For speed, permutation shuffles default to a cheaper backend unless the exact backend is explicitly requested; when this is done, observed and permuted statistics must be generated by the same backend for the calibration in Eq. 6 to be interpretable.

### S4. Bias-correction options

`bias_correction="huber_length_intron"` (default) applies Eq. 4. `bias_correction="none"` sets $R_{g} = \Delta_{g}$ for expressed genes, bypassing the regression entirely; Eqs. 5 and 6 are otherwise unchanged. The median-centering fallback described in the main text is applied automatically whenever the Huber fit cannot be performed.

### S5. Summary of options

| Option                      | Values (default in bold)                                        | Affects                              | Notes                                           |
|:----------------------------|:----------------------------------------------------------------|:-------------------------------------|:------------------------------------------------|
| Observational unit          | **pseudobulk (Mode A)** / cell (B) / kNN-smoothed cell (C)      | $\bar{U},\bar{S}$ in Eqs. 2–3        | Only A supports replicate-level permutation     |
| DE backend                  | **PyDESeq2** / Scanpy t-test / Wilcoxon / mixed LM / Memento    | ${logFC}_{g},\ p_{g}^{adj}$ in Eq. 5 | Must be identical in observed and permuted runs |
| `gamma_method`              | **heuristic\_shrink** / robust\_median / empirical\_bayes / raw | Eq. 2 → S1, S2                       | EB recommended for small reference groups       |
| `prior_weight` $\eta$       | **5.0**                                                         | Eq. 2                                | Larger = stronger shrinkage to $\rho_{0}$       |
| `bias_correction`           | **huber\_length\_intron** / none                                | Eq. 4                                | Median-centering fallback is automatic           |
| $w_{FC},w_{V},w_{P}$        | **1, 1, 1**                                                     | Eq. 5                                | Sensitivity analysis recommended                |
| `ranking_mode`              | **composite** / nascent\_excess                                 | Ranking                              | `nascent_excess` ranks by $R_{g}$ alone         |
| `n_perm` $B$                | **100**                                                         | Eq. 6                                | FDR usable only when $B_{eff} \geq 100$         |
| `min_total_counts` $\kappa$ | **50**                                                          | $\mathcal{G}_{valid}$                | Genes below threshold get $R_{g} = 0$           |

### S6. Equation index (for the software documentation)

| Eq.   | Quantity                                       | Field in output                                                                            |
|:------|:-----------------------------------------------|:-------------------------------------------------------------------------------------------|
| 1     | Kinetic model (context only)                   | —                                                                                          |
| 2     | Reference $U/S$ ratio ${\widehat{\rho}}_{R,g}$ | `effective_gamma` (opt-in)                                                                 |
| 3     | Raw unspliced excess $\Delta_{g}$              | `unspliced_excess_delta` (legacy: `velocity_delta_raw`)                                        |
| 4     | Huber bias fit and residual $R_{g}$            | `unspliced_excess_residual`; coefficients in `uns["scatrans"]["diagnostics"]`              |
| 5     | Composite Active Evidence Score $A_{g}$        | `active_score`                                                                             |
| 6     | Permutation p-value and FDR                    | `unspliced_excess_pval` / `unspliced_excess_fdr`; `active_score_pval` / `active_score_fdr` |
| C1–C3 | Moment smoothing (Mode C)                      | layers `Mu`, `Ms`                                                                          |
| S1–S2 | Alternative ratio estimators                   | `gamma_info` diagnostics                                                                   |
