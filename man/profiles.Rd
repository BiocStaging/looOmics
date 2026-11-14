\name{profiles}
\alias{profiles}
\title{Get the flag profile from a LOOResult}
\usage{
profiles(object)
}
\arguments{
  \item{object}{A \code{\linkS4class{LOOResult}}.}
}
\value{
  The FlagProfile computed by \code{\link{profile_flags}} (NULL if not yet
  computed).
}
\description{
  Returns the metadata-profiling results computed by
  \code{\link{profile_flags}}.
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

