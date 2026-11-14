\name{LOOResult-class}
\alias{LOOResult}
\alias{LOOResult-class}
\title{Leave-one-out sensitivity analysis result}
\description{
  S4 container for the output of \code{\link{loo_sensitivity}}. Holds the
  full-sample fit, leave-one-out betas/SEs, the long shift table (per
  feature, term, and left-out sample: beta, beta_loo, delta, dfbetas, p_loo,
  q_loo, sig_flip), per-sample influence aggregates, failed-refit records,
  flags (\code{\link{flag_influential}}), metadata profiles
  (\code{\link{profile_flags}}), and diagnostics (runtime, future plan,
  RNG seed, sample metadata, per-term influence).
}
\value{
  An S4 object of class \code{LOOResult}, containing the full-sample
  model results, leave-one-out estimates, influence statistics,
  failed-refit records, flags, metadata profiles, and diagnostics.
}
\examples{
set.seed(1)
samples <- paste0("S", seq_len(8))
x <- matrix(rnorm(10 * 8), nrow = 10,
            dimnames = list(paste0("gene", 1:10), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
                 row.names = samples)
result <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
is(result, "LOOResult")
}

\section{Slots}{
  \describe{
    \item{\code{call}}{The matched call.}
    \item{\code{engine}}{Engine id.}
    \item{\code{engine_info}}{Engine metadata (package versions, terms).}
    \item{\code{coefficients}}{Full-sample fit: feature, term, beta, se, p.}
    \item{\code{loo_betas}, \code{loo_se}}{LOO fits; rows are
      feature::term, columns are left-out samples.}
    \item{\code{shifts}}{Long shift table.}
    \item{\code{influence}}{Per-sample aggregate scores.}
    \item{\code{failures}}{Failed LOO refits (fail-soft record).}
    \item{\code{flags}}{Output of flag_influential().}
    \item{\code{profiles}}{Output of profile_flags().}
    \item{\code{diagnostics}}{Runtime, seed, future plan, metadata,
      per-term influence.}
    \item{\code{input_hash}}{Digest of the input assay and metadata.}
  }
}
