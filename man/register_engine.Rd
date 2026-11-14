\name{register_engine}
\alias{register_engine}
\title{Register a custom fitting engine}
\usage{
register_engine(name, fit_fun)
}
\arguments{
  \item{name}{Engine id used in \code{loo_sensitivity(engine = ...)}.}
  \item{fit_fun}{function(input, design, contrast, block, ...) returning
    list(betas, se, p) matrices (features x terms).}
}
\examples{
example_fit <- function(input, design, contrast = NULL, block = NULL, ...) {
  features <- rownames(input@assay)
  terms <- colnames(design)
  estimates <- matrix(0, length(features), length(terms),
                      dimnames = list(features, terms))
  list(betas = estimates, se = estimates, p = estimates,
       colnames = terms)
}
register_engine("example_engine", example_fit)
rm(list = "example_engine",
   envir = getFromNamespace(".engine_registry", "looOmics"))
}

\value{
  Invisibly returns the registered engine name.
}

\description{
  Registers a custom fitting engine for use with
  \code{\link{loo_sensitivity}}. The function receives the OmicsInput, the
  design matrix, and optional contrast/block arguments, and must return
  coefficient, standard-error, and p-value matrices with features as rows
  and model terms as columns.
}
