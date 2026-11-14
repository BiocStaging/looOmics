\name{autoplot.LOOResult}
\alias{autoplot.LOOResult}
\alias{plot,LOOResult,missing-method}
\title{Plot a LOOResult}
\usage{
\method{autoplot}{LOOResult}(object, type = c("shift", "scree"), ...)
\S4method{plot}{LOOResult,missing}(x, y, ...)
}
\arguments{
  \item{object, x}{A \code{\linkS4class{LOOResult}}.}
  \item{type}{One of "shift" (per-sample score barplot, flagged in red) or
    "scree" (ranked influence scores).}
  \item{y, ...}{Unused / passed to layers.}
}
\value{
  A ggplot object.
}
\description{
  Visualizes per-sample influence scores, highlighting flagged samples.
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
autoplot(res, type = "shift")
}

