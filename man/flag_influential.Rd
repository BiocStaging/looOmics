\name{flag_influential}
\alias{flag_influential}
\title{Flag influential samples}
\usage{
flag_influential(object, method = c("mad", "threshold", "topk", "quantile"),
                 threshold = NULL, k = NULL, q = 0.9, metric = "max_dfbetas")
}
\arguments{
  \item{object}{A \code{\linkS4class{LOOResult}} with influence metrics.}
  \item{method}{One of "mad" (robust z on max |DFBETAS|), "threshold"
    (absolute cutoff, default 1), "topk", or "quantile".}
  \item{threshold}{Numeric; MAD cutoff (default 5) or absolute cutoff
    (default 1) depending on method.}
  \item{k}{Number of top samples for "topk".}
  \item{q}{Quantile for "quantile" (default 0.9).}
  \item{metric}{Aggregate to flag on: "max_dfbetas" (default),
    "n_sig_flips", or "median_abs_delta".}
}
\value{
  The input object with the populated \code{flags} slot (columns: sample,
  score, metric, method, flagged).
}
\description{
  Applies a flagging rule to the per-sample influence scores. Flagged
  samples merit inspection; they are not automatically invalid.
}
\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
                 row.names = samples)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
res <- flag_influential(res, method = "topk", k = 2)
res@flags
}
