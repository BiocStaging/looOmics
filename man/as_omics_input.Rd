\name{as_omics_input}
\alias{as_omics_input}
\alias{as_omics_input,matrix-method}
\alias{as_omics_input,data.frame-method}
\alias{as_omics_input,SummarizedExperiment-method}
\alias{as_omics_input,ExpressionSet-method}
\alias{as_omics_input,phyloseq-method}
\title{Coerce an omics object into the unified input container}
\usage{
as_omics_input(x, formula, engine = "auto", transform = "none",
               pseudocount = 1, ...)
}
\arguments{
  \item{x}{An omics object.}
  \item{formula}{Model formula evaluated in the sample covariate table;
    the left-hand side is implicit.}
  \item{engine}{Engine id (see \code{\link{loo_sensitivity}}). "auto"
    resolves deterministically to lm. Select limma or DESeq2 explicitly when
    appropriate for the assay type.}
  \item{transform}{Transformation applied to the assay: "none" (default),
    "log2", or "clr". The pseudocount is added exactly once, inside this
    transformation.}
  \item{pseudocount}{Pseudocount added before log2/clr transforms
    (default 1).}
  \item{...}{For the matrix method: \code{col_data} (sample-level
    covariates) and optionally \code{row_data}. For
    SummarizedExperiment: \code{assay_name}.}
}
\value{
  An \code{\linkS4class{OmicsInput}} object.
}
\description{
  Dispatches on the class of \code{x}. Matrix and data.frame methods are
  always available; SummarizedExperiment, ExpressionSet and phyloseq
  methods require the corresponding packages. SummarizedExperiment is the
  recommended primary interface for Bioconductor workflows. The model formula
  is stored for formula-driven engines such as DESeq2. Assay values must be
  numeric and finite; missing and infinite values should be handled upstream.
}
\examples{
mat <- matrix(1:12, 3, 4, dimnames = list(paste0("f", 1:3), paste0("S", 1:4)))
cd <- data.frame(grp = factor(c("A","A","B","B")), row.names = colnames(mat))
inp <- as_omics_input(mat, ~ grp, col_data = cd)
}
