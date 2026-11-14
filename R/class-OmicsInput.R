#' @include generics.R
NULL

#' Unified omics input container
#'
#' S4 container holding an assay matrix (features by samples), feature
#' metadata, sample-level covariates, engine identifier, transformation,
#' provenance class, and model formula.
#'
#' @slot assay Numeric matrix with features in rows and samples in columns.
#' @slot row_data Feature-level metadata.
#' @slot col_data Sample-level covariates.
#' @slot engine Engine identifier.
#' @slot transform Transformation applied to the assay.
#' @slot data_class Class of the original input object.
#' @slot formula Deparsed model formula used by formula-driven engines.
#' @name OmicsInput-class
#' @rdname OmicsInput-class
#' @aliases OmicsInput
#' @docType class
#' @return An S4 object of class `OmicsInput`, containing the
#'   assay matrix, feature metadata, sample metadata, engine,
#'   transformation, provenance, and model formula.
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
#' is(input, "OmicsInput")
NULL

# OmicsInput: unified features x samples container -------------------------

setClass("OmicsInput",
    slots = c(
        assay      = "matrix",      # features (rows) x samples (cols)
        row_data   = "data.frame",  # feature metadata
        col_data   = "data.frame",  # sample-level covariates
        engine     = "character",
        transform  = "character",   # "", "log2", "clr" - recorded, applied by adapter
        data_class = "character",  # provenance
        formula    = "character"   # deparsed model formula (needed by formula-based engines)
    ),
    prototype = list(
        assay      = matrix(0, 0, 0),
        row_data   = data.frame(),
        col_data   = data.frame(),
        engine     = "lm",
        transform  = "",
        data_class = "matrix",
        formula    = ""
    )
)

setValidity("OmicsInput", function(object) {
    msg <- character(0)
    if (ncol(object@assay) != nrow(object@col_data))
        msg <- c(msg, "ncol(assay) must equal nrow(col_data)")
    if (nrow(object@assay) != nrow(object@row_data))
        msg <- c(msg, "nrow(assay) must equal nrow(row_data)")
    if (anyDuplicated(colnames(object@assay)))
        msg <- c(msg, "sample names must be unique")
    if (anyNA(object@assay))
        msg <- c(msg, "assay must not contain NA (filter or impute upstream)")
    if (length(msg)) msg else TRUE
})

# Constructor (not exported; use as_omics_input())
new_OmicsInput <- function(assay, row_data, col_data, engine, transform,
                           data_class, formula = "") {
    methods::new("OmicsInput",
        assay = assay, row_data = row_data, col_data = col_data,
        engine = engine, transform = transform, data_class = data_class,
        formula = formula)
}
