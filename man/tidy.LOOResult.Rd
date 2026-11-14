\name{tidy.LOOResult}
\alias{tidy.LOOResult}
\title{Tidy a LOOResult into a tibble}
\usage{
\method{tidy}{LOOResult}(x, what = c("shifts", "influence", "flags"), ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{LOOResult}}.}
  \item{what}{One of "shifts" (default), "influence", or "flags".}
  \item{...}{Unused.}
}
\value{
  A tibble.
}
\description{
  Extracts the shift table, influence scores, or flags as a tibble.
}

\examples{
set.seed(1)
samples <- paste0("S", seq_len(12))
x <- matrix(rnorm(20 * 12), nrow = 20,
            dimnames = list(paste0("gene", 1:20), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
                 row.names = samples)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
head(tidy(res, what = "influence"))
}

