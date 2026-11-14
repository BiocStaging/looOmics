\name{run_engine}
\alias{run_engine}
\alias{run_engine,OmicsInput-method}
\title{Fit the full-sample model with a registered engine}
\usage{
run_engine(input, design, contrast = NULL, block = NULL, ...)
\S4method{run_engine}{OmicsInput}(input, design, contrast = NULL,
  block = NULL, ...)
}
\arguments{
  \item{input}{An \code{\linkS4class{OmicsInput}}.}
  \item{design}{Model matrix (rows are samples).}
  \item{contrast, block}{Optional contrast specification and blocking
    factor, passed through to the engine.}
  \item{...}{Engine-specific arguments.}
}
\value{
  A list with elements betas, se, p (features x terms matrices) and
  colnames (model terms).
}
\examples{
set.seed(1)
samples <- paste0("S", seq_len(8))
x <- matrix(rnorm(10 * 8), nrow = 10,
            dimnames = list(paste0("gene", 1:10), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
                 row.names = samples)
input <- as_omics_input(x, ~ condition, col_data = cd)
design <- stats::model.matrix(~ condition, data = cd)
fit <- run_engine(input, design)
dim(fit$betas)
}

\description{
  Internal dispatch used by \code{\link{loo_sensitivity}}; also the
  extension point for custom engines together with
  \code{\link{register_engine}}. The method dispatches on the engine id
  stored on the input: built-in engines first, then registered custom
  engines.
}
