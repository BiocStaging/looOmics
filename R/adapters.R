#' @include generics.R
NULL
# Input adapters ------------------------------------------------------------
# All adapter methods are registered unconditionally: S4 dispatch is by
# class name, so methods for Suggests classes (SummarizedExperiment,
# ExpressionSet, phyloseq) are defined even when those packages are absent;
# the method bodies call requireNamespace_or_stop() for a clean runtime
# error if the class is used without the package installed.

#' @rdname as_omics_input
#' @export
setMethod("as_omics_input", "matrix",
          function(x, formula, engine = "auto", transform = "none",
                   pseudocount = 1, col_data, row_data = NULL, ...) {
            engine <- resolve_engine(engine)
            if (!is.numeric(x) || is.complex(x))
              stop("The assay matrix must contain real numeric values.", call. = FALSE)
            if (any(!is.finite(x)))
              stop("The assay matrix must contain only finite values; handle missing or infinite values before analysis.",
                   call. = FALSE)
            if (missing(col_data))
              stop("The matrix method requires 'col_data' (sample-level covariates).",
                   call. = FALSE)
            x <- apply_transform(x, transform, pseudocount = pseudocount)
            if (is.null(rownames(x))) rownames(x) <- paste0("feat", seq_len(nrow(x)))
            if (is.null(colnames(x))) colnames(x) <- paste0("S", seq_len(ncol(x)))
            if (is.null(rownames(col_data)))
              stop("'col_data' must have rownames matching the sample names.", call. = FALSE)
            if (!identical(rownames(col_data), colnames(x)))
              stop("'col_data' rownames must match colnames(x).", call. = FALSE)
            if (is.null(row_data)) row_data <- data.frame(row.names = rownames(x))
            new_OmicsInput(assay = x, row_data = row_data, col_data = as.data.frame(col_data),
                           engine = engine, transform = transform, data_class = "matrix",
                           formula = paste(deparse(formula), collapse = " "))
          })

#' @rdname as_omics_input
#' @export
setMethod("as_omics_input", "data.frame",
          function(x, formula, engine = "auto", transform = "none",
                   pseudocount = 1, col_data, row_data = NULL, ...) {
            as_omics_input(as.matrix(x), formula = formula, engine = engine,
                           transform = transform, col_data = col_data,
                           row_data = row_data, pseudocount = pseudocount, ...)
          })

# Single pseudocount mechanism: apply_transform() adds it; adapters must NOT
# add their own. CLR centers each sample (column) across features.
apply_transform <- function(x, transform, pseudocount = 1) {
  switch(transform,
         none = x,
         log2 = {
           if (any(x < 0)) stop("log2 transform requires non-negative values.", call. = FALSE)
           log2(x + pseudocount)
         },
         clr = {
           if (any(x < 0)) stop("clr transform requires non-negative values.", call. = FALSE)
           log_x <- log(x + pseudocount)
           sweep(log_x, 2, colMeans(log_x), FUN = "-")
         },
         stop("Unknown transform '", transform, "' (use none/log2/clr).", call. = FALSE))
}

# --- Bioconductor / phyloseq adapters ---------------------------------------

#' @rdname as_omics_input
#' @export
setMethod("as_omics_input", "SummarizedExperiment",
          function(x, formula, engine = "auto", transform = "none",
                   pseudocount = 1, assay_name = "counts", ...) {
            requireNamespace_or_stop("SummarizedExperiment", "as_omics_input")
            a <- as.matrix(SummarizedExperiment::assay(x, assay_name))
            cd <- as.data.frame(SummarizedExperiment::colData(x))
            rd <- as.data.frame(SummarizedExperiment::rowData(x))
            as_omics_input(a, formula = formula, engine = engine, transform = transform,
                           col_data = cd, row_data = rd, pseudocount = pseudocount, ...)
          })

#' @rdname as_omics_input
#' @export
setMethod("as_omics_input", "ExpressionSet",
          function(x, formula, engine = "auto", transform = "none", pseudocount = 1, ...) {
            requireNamespace_or_stop("Biobase", "as_omics_input")
            as_omics_input(as.matrix(Biobase::exprs(x)), formula = formula, engine = engine,
                           transform = transform,
                           col_data = as.data.frame(Biobase::pData(x)),
                           pseudocount = pseudocount, ...)
          })

#' @rdname as_omics_input
#' @export
setMethod("as_omics_input", "phyloseq",
          function(x, formula, engine = "auto", transform = "clr",
                   pseudocount = 1, ...) {
            requireNamespace_or_stop("phyloseq", "as_omics_input")
            otu <- phyloseq::otu_table(x)
            a <- if (phyloseq::taxa_are_rows(otu)) as(otu, "matrix") else t(as(otu, "matrix"))
            cd <- as.data.frame(phyloseq::sample_data(x))
            rd <- as.data.frame(phyloseq::tax_table(x))
            as_omics_input(a, formula = formula, engine = engine, transform = transform,
                           col_data = cd, row_data = rd, pseudocount = pseudocount, ...)
          })