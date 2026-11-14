# looOmics

Leave-one-out sensitivity analysis for omics differential models.

`looOmics` revalidates differential-analysis pipelines on gene expression
and microbiome data by refitting the model once per left-out sample,
quantifying coefficient shifts, significance changes, and within-term
feature-rank changes for every feature and model term, flagging influential
samples, and profiling those samples against biological and
technical metadata covariates.

**This is a sensitivity analysis, not a bad-sample detector.** An
influential sample may reflect genuine biological heterogeneity; flagged
samples merit inspection, not automatic removal.

## Workflow

```r
library(looOmics)

# SummarizedExperiment is the primary (recommended) input interface.
# This example assumes a normalized, log-transformed logcounts assay.
res <- loo_sensitivity(se, ~ condition + batch, engine = "limma",
                       assay_name = "logcounts")

# Flag influential samples
res <- flag_influential(res, method = "mad")

# Profile flags against metadata (batch, depth, sex, age, ...)
res <- profile_flags(res, model_covariates = c("condition", "batch"))

# Inspect
show(res)
influence_by_term(res)   # per-term influence scores
profiles(res)            # biological / technical metadata associations
plot(res)                # per-sample influence scores

# Parallel execution is controlled by the user via future
future::plan(future::multisession, workers = 4)
res <- loo_sensitivity(se, ~ condition, engine = "limma",
                       assay_name = "logcounts")
future::plan(future::sequential)
```

## Engines

| Engine | Status | Notes |
|---|---|---|
| `lm` | Implemented | Exact match to `stats::lm` |
| `glm` | Implemented | Family-aware inference (t-based for gaussian/quasi, z-based for binomial/poisson with fixed dispersion) |
| `limma` | Implemented | Requires limma; select explicitly for normalized log-expression data |
| `deseq2` | Implemented | Formula-driven; requires untransformed non-negative integer counts; size factors and dispersions re-estimated for every fit |
| custom | Extensible | `register_engine()` |

`engine = "auto"` resolves deterministically to `lm`; it does not infer the
appropriate engine from the assay or depend on which optional packages happen
to be installed. For raw RNA-seq counts, select `engine = "deseq2"` explicitly.
For normalized, log-transformed expression values, select `engine = "limma"`.

## Installation

```r
# Install from the Bioconductor development branch after acceptance
BiocManager::install("looOmics")

# Development
remotes::install_github("<org>/looOmics")
```

## Parallelism

All heavy computation routes through `future.apply`; the package never
changes the plan itself. Results are identical across sequential,
multisession, and cluster plans under the same seed.

## Biological and technical metadata profiling

Metadata profiling is a core part of the workflow, not just an optional
outlier annotation. `profile_flags()` evaluates whether influential-sample
status is associated with available sample-level variables such as biological
subtype, tumor purity, treatment group, sequencing depth, extraction protocol,
batch, or other technical quality measures. Continuous and categorical
covariates receive appropriate exploratory univariate tests, with
Benjamini-Hochberg adjustment. An optional adjusted logistic model estimates
associations across covariates when the sample size supports it. Variables
used in the differential model are retained and labeled in the output.
These analyses are associational, not causal; influential samples are not
assumed to be erroneous.

## Scope

- Leave-one-sample-out only (feature knock-outs and Bayesian LOO are out
  of scope).
- The full pipeline — including normalization and model-specific
  estimation steps — is rerun for every leave-one-out fit where
  appropriate.
- Influence is quantified as deletion-based standardized coefficient
  shifts: for ordinary linear models these coincide with classical
  DFBETAS; for other engines they are the analogous quantity.

## License

MIT.
