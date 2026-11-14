\name{profile_flags}
\alias{profile_flags}
\title{Profile flagged samples against biological and technical covariates}
\usage{
profile_flags(object, covariates = NULL, model_covariates = NULL,
              multivariable = TRUE, min_cell_size = 3,
              seed = 1L, fisher_B = 2000L)
}
\arguments{
  \item{object}{A \code{\linkS4class{LOOResult}} with flags computed.}
  \item{covariates}{Character vector of sample-metadata columns to profile.
    Defaults to all available metadata, including variables in the differential
    model.}
  \item{model_covariates}{Character vector of variables in the original
    differential model. Defaults to variables inferred from its formula.
    These variables remain profiled and are labeled as in-model.}
  \item{multivariable}{Logical; fit an exploratory adjusted logistic regression
    when the sample size and design support it.}
  \item{min_cell_size}{Minimum observations per categorical level before rare
    levels are collapsed into \code{Other}.}
  \item{seed}{Integer seed for simulated Fisher tests. The caller's random
    number generator state is restored after each simulated test.}
  \item{fisher_B}{Number of Monte Carlo replicates for simulated Fisher tests.}
}
\value{
  The input object with the populated \code{profiles} slot; retrieve it
  with \code{\link{profiles}}.
}
\description{
  Evaluates associations between influential-sample status and biological or
  technical sample metadata. Continuous covariates are assessed with
  Wilcoxon rank-sum tests; categorical covariates use chi-square tests or
  simulated Fisher tests when expected counts are sparse. Univariate
  p-values are adjusted using Benjamini-Hochberg. An optional multivariable
  logistic regression estimates adjusted associations across selected
  metadata covariates. Numeric covariates are standardized to one standard
  deviation in the adjusted model; factor coefficients are relative to their
  reference levels. Variables used in the differential model are retained
  and labeled because their association with influence can be biologically
  informative. This is exploratory association analysis, not causal inference.
  Covariates that cannot be tested are reported with reasons.
}

\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(
  condition = factor(rep(c("A", "B"), each = 6)),
  batch = factor(rep(c("batch1", "batch2"), 6)),
  age = seq(40, 62, length.out = 12),
  row.names = samples
)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
res <- flag_influential(res, method = "topk", k = 3)
res <- profile_flags(res, covariates = c("batch", "age"),
                     multivariable = FALSE)
profiles(res)
}

