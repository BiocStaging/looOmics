# looOmics 0.99.0

Bioconductor pre-submission development version.

## New features

- `loo_sensitivity()`: sample-wise leave-one-out revalidation of omics
  differential models with coefficient-shift tracking (delta, standardized
  deletion shifts, significance flips, and within-term feature-rank changes).
- `flag_influential()`: threshold / MAD / top-k / quantile flagging rules.
- `profile_flags()`: univariate (Fisher / chi-square / Wilcoxon) and
  multivariable (logistic regression, ORs) profiling of biological and
  technical sample metadata, including variables in the original model;
  BH adjustment, expected-count-based test selection, deterministic simulated
  Fisher tests, missing-data handling, rare-level collapsing, and separation
  and small-sample safeguards.
- `influence_by_term()`: per-term influence summaries.
- Engines: `lm`, `glm` (family-aware inference), `limma`, `deseq2`
  (formula-driven, raw counts, per-fit size-factor and dispersion
  re-estimation), plus `register_engine()` for custom engines.
- Input adapters: matrix/data.frame always; SummarizedExperiment,
  ExpressionSet, phyloseq when installed. `engine = "auto"` resolves
  deterministically to lm, independent of optional package availability.
- `future`-based parallel execution; results identical across plans under
  the same seed. Fail-soft per-sample error capture.
- Single pseudocount mechanism for log2/CLR transforms.

## Bug fixes (pre-release review)

- CLR transform: corrected column orientation; each sample (column) is now
  centered across features, preserving the features x samples layout.
- Pseudocount: applied exactly once (inside `apply_transform()`); adapters
  no longer add their own.
- Engine selection: `auto` now resolves deterministically to lm instead of
  choosing an engine based on installed optional packages; count-based RNA-seq
  and normalized expression workflows must select DESeq2 or limma explicitly.
- DESeq2 engine: redesigned around the stored model formula and raw integer
  counts instead of the design matrix; transformed, fractional, negative, or
  non-finite inputs are rejected rather than rounded; size factors and
  dispersions are re-estimated for every full and leave-one-out fit.
- Feature ranks: added full-data and leave-one-out ranks within each model
  term, with positive rank changes indicating a weaker numerical rank after
  sample deletion.
- Metadata profiling: original model variables are inferred from the formula,
  retained in the default profile, and labeled as in-model; simulated Fisher
  tests use a reproducible local seed without changing the caller RNG state.
- glm engine: family-aware dispersion (estimated for gaussian/quasi,
  fixed at 1 for binomial/poisson) and inference (t-based vs z-based).
- Metadata profiling: expected-count-based Fisher selection, reproducible
  simulated Fisher tests, NA handling, dropped-covariate reporting with
  reasons, retention and labeling of original model covariates, and
  separation and small-sample safeguards for the adjusted model.
- `n_features_shifted` counts unique features rather than feature-term
  combinations.
