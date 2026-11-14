\name{reexports}
\alias{autoplot}
\alias{tidy}
\title{Objects re-exported from other packages}
\value{
  `autoplot()` returns a ggplot object; `tidy()` returns a tibble.
}

\description{
  These objects are imported from other packages and re-exported by
  looOmics so that \code{\link[ggplot2]{autoplot}} and
  \code{\link[broom]{tidy}} can be used directly on
  \code{\linkS4class{LOOResult}} objects. See
  \code{\link{autoplot.LOOResult}} and \code{\link{tidy.LOOResult}} for
  the looOmics methods.
}

\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
                 row.names = samples)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
autoplot(res)
head(tidy(res))
}

