\name{influence_by_term}
\alias{influence_by_term}
\title{Per-term influence scores}
\usage{
influence_by_term(object)
}
\arguments{
  \item{object}{A \code{\linkS4class{LOOResult}}.}
}
\value{
  A data.frame with columns term, sample, max_abs_shift, n_sig_flips.
}
\description{
  Returns the per-sample influence aggregates computed separately for each
  model term. The overall \code{influence} slot pools all terms; this
  accessor shows which model terms drive a sample's score. Pass
  \code{terms} filtering downstream if the intercept is not of interest.
}

\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
                 row.names = samples)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
head(influence_by_term(res))
}

