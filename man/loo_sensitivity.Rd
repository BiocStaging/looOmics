\name{loo_sensitivity}
\alias{loo_sensitivity}
\title{Sample-wise leave-one-out sensitivity analysis for omics models}
\usage{
loo_sensitivity(x, formula, engine = "auto", parallel = TRUE,
                seed = TRUE, ...)
}
\arguments{
  \item{x}{Omics object: matrix/data.frame (features x samples, with
    \code{col_data}), or an OmicsInput. SummarizedExperiment, ExpressionSet
    and phyloseq are supported when those packages are installed.}
  \item{formula}{Model formula evaluated in the sample covariate table; the
    left-hand side is implicit (each feature is the response).}
  \item{engine}{Engine id: "auto" (default; resolves deterministically to
    lm), "lm", "glm", "limma", "deseq2", or a name registered via
    \code{\link{register_engine}}.}
  \item{parallel}{Logical, default TRUE. If FALSE, runs sequentially via
    lapply. The parallel strategy itself is owned by the user via
    \code{future::plan()}.}
  \item{seed}{Integer or TRUE, passed to future.apply::future_lapply.}
  \item{...}{Input-adapter arguments (such as \code{col_data},
    \code{assay_name}, \code{transform}, and \code{pseudocount}) and
    engine-specific arguments (such as \code{family} for glm). Adapter
    arguments are not forwarded to the model-fitting engine.}
}
\value{
  A \code{\linkS4class{LOOResult}} with full-sample coefficients, LOO
  betas/SEs, the long shift table, and per-sample influence scores.
}
\description{
  For each sample i, refits the model on the remaining n-1 samples and
  records coefficient shifts for every feature and model term, then computes
  influence metrics per sample. The entire relevant pipeline (normalization
  and model-specific estimation steps) is rerun for every leave-one-out fit
  where appropriate.
}
\details{
  The column \code{dfbetas} in the shift table is a deletion-based
  standardized coefficient shift, (beta_full - beta_loo) / se_loo. The
  shift table also reports full and leave-one-out p-values, BH-adjusted
  q-values, and within-term feature ranks. A positive rank_change means the
  feature moved to a numerically larger (weaker) rank after sample deletion. For
  ordinary linear models this coincides with classical DFBETAS; for other
  engines it is the analogous standardized deletion shift and should be
  interpreted as such. This is a sensitivity analysis for detecting
  influential samples, not an automatic bad-sample detector: an influential
  sample may reflect genuine biological heterogeneity.

  Failed leave-one-out refits (e.g. rank-deficient designs after sample
  deletion) are recorded in the \code{failures} slot and do not abort the
  run.
}
\examples{
set.seed(1)
mat <- matrix(rnorm(100 * 12), 100, 12,
              dimnames = list(paste0("g", 1:100), paste0("S", 1:12)))
cd <- data.frame(condition = factor(rep(c("A", "B"), 6)),
                 batch = factor(rep(1:2, 6)), row.names = colnames(mat))
res <- loo_sensitivity(mat, ~ condition, col_data = cd)
flag_influential(res)
}
