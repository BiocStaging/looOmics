\name{OmicsInput-class}
\alias{OmicsInput}
\alias{OmicsInput-class}
\title{Unified omics input container}
\description{
  S4 container holding the assay (features x samples), feature metadata,
  sample-level covariates, engine id, applied transformation, provenance
  class, and the deparsed model formula (used by formula-driven engines
  such as DESeq2). Constructed via \code{\link{as_omics_input}}.
}
\value{
  An S4 object of class \code{OmicsInput}, containing the assay matrix,
  feature metadata, sample metadata, engine, transformation,
  provenance information, and model formula.
}
\examples{
mat <- matrix(
  seq_len(12), nrow = 3,
  dimnames = list(paste0("feature", 1:3), paste0("sample", 1:4))
)
metadata <- data.frame(
  condition = factor(c("control", "control", "treated", "treated")),
  row.names = colnames(mat)
)
input <- as_omics_input(mat, ~ condition, col_data = metadata)
is(input, "OmicsInput")
}

\section{Slots}{
  \describe{
    \item{\code{assay}}{Numeric matrix, features x samples.}
    \item{\code{row_data}}{Feature metadata (data.frame).}
    \item{\code{col_data}}{Sample-level covariates (data.frame).}
    \item{\code{engine}}{Engine id (character).}
    \item{\code{transform}}{Recorded transformation.}
    \item{\code{data_class}}{Provenance class of the original object.}
    \item{\code{formula}}{Deparsed model formula.}
  }
}
