# Generics and helpers ------------------------------------------------------

requireNamespace_or_stop <- function(pkg, fn) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required for ", fn,
             "(). Install it from CRAN or Bioconductor, as appropriate.",
             call. = FALSE)
    invisible(TRUE)
}

#' Coerce an omics object into looOmics' unified input container
#'
#' Dispatches on the class of \code{x}. Methods are registered for
#' \code{matrix}, \code{data.frame} (features x samples), and — when the
#' corresponding packages are installed — \code{SummarizedExperiment},
#' \code{ExpressionSet} and \code{phyloseq}.
#'
#' @param x An omics object.
#' @param engine Engine id: "auto" (default; resolves deterministically to
#'   "lm"), "lm", "glm", "limma", or "deseq2". DESeq2 requires raw
#'   non-negative integer counts and uses the stored model formula.
#' @param formula Model formula evaluated within the sample covariate table;
#'   the left-hand side is implicit (each feature is the response).
#' @param transform Transformation applied to the assay: "none" (default),
#'   "log2", or "clr".
#' @param ... Passed to methods. For the matrix method: \code{col_data}
#'   (sample-level covariates) and optionally \code{row_data}.
#'
#' @return An \code{\linkS4class{OmicsInput}} object.
#' @examples
#' mat <- matrix(
#'   seq_len(12), nrow = 3,
#'   dimnames = list(paste0("feature", 1:3), paste0("sample", 1:4))
#' )
#' metadata <- data.frame(
#'   condition = factor(c("control", "control", "treated", "treated")),
#'   row.names = colnames(mat)
#' )
#' input <- as_omics_input(mat, ~ condition, col_data = metadata)
#' assay(input)
#' @export
setGeneric("as_omics_input",
    function(x, formula, engine = "auto", transform = "none", pseudocount = 1, ...)
        standardGeneric("as_omics_input"))

# Resolve "auto" deterministically. Engine choice must not depend on which
# optional packages happen to be installed on a user's system.
resolve_engine <- function(engine) {
    if (!identical(engine, "auto")) return(engine)
    "lm"
}

