\name{print.FlagProfile}
\alias{print.FlagProfile}
\title{Print a biological and technical metadata profile}
\usage{
\method{print}{FlagProfile}(x, ...)
}
\arguments{
  \item{x}{A \code{FlagProfile} object returned by \code{profiles()}.}
  \item{...}{Additional arguments (currently unused).}
}
\value{
  The input object, invisibly.
}
\description{
  Prints univariate associations between influential-sample status and sample
  metadata, optional multivariable logistic-regression results, and reasons
  why selected covariates could not be tested.
}

\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(
  condition = factor(rep(c("A", "B"), each = 6)),
  batch = factor(rep(c("batch1", "batch2"), 6)),
  row.names = samples
)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
res <- flag_influential(res, method = "topk", k = 3)
res <- profile_flags(res, multivariable = FALSE)
print(profiles(res))
}

