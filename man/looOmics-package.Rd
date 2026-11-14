\name{looOmics-package}
\alias{looOmics-package}
\alias{looOmics}
\docType{package}
\title{Leave-One-Out Sensitivity Analysis for Omics Differential Models}
\description{
  looOmics performs sample-wise leave-one-out (LOO) revalidation of
  differential-analysis pipelines on gene expression and microbiome data. For
  each sample it refits the model on the remaining n-1 samples, records
  coefficient shifts for every feature and model term, computes per-sample
  influence metrics, flags influential samples, and profiles flagged samples
  against biological and technical metadata covariates.

  The analysis is a sensitivity analysis: an influential sample may reflect
  genuine biological heterogeneity and is not automatically a "bad" sample.
}
\details{
  The main entry point is \code{\link{loo_sensitivity}}. Inputs are built with
  \code{\link{as_omics_input}}, which dispatches on the input class
  (matrix, data.frame, and - when installed - SummarizedExperiment,
  ExpressionSet, phyloseq). Results are stored in a
  \code{\linkS4class{LOOResult}} object; influence summaries per model term
  are available through \code{\link{influence_by_term}}.

  Workflow:
  \enumerate{
    \item \code{res <- loo_sensitivity(se, ~ condition + batch)}
    \item \code{res <- flag_influential(res)}
    \item \code{res <- profile_flags(res)}
    \item \code{plot(res)}
  }
}
\keyword{package}
