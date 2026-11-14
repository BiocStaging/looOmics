#' @include generics.R
NULL

#' Leave-one-out sensitivity analysis result
#'
#' S4 container for the output of \code{loo_sensitivity()}. Stores the
#' full-sample fit, leave-one-out estimates, feature-level shifts, per-sample
#' influence aggregates, failed refits, flags, metadata profiles, and
#' diagnostics.
#'
#' @slot call Matched call to \code{loo_sensitivity()}.
#' @slot engine Engine identifier.
#' @slot engine_info Engine metadata and model terms.
#' @slot coefficients Full-sample feature-level coefficients, standard errors,
#'   and p-values.
#' @slot loo_betas Leave-one-out coefficient estimates.
#' @slot loo_se Leave-one-out standard errors.
#' @slot shifts Long table of feature-, term-, and sample-level shifts.
#' @slot influence Per-sample influence aggregates.
#' @slot failures Records of failed leave-one-out refits.
#' @slot flags Results from \code{flag_influential()}.
#' @slot profiles Results from \code{profile_flags()}.
#' @slot diagnostics Runtime, model, metadata, and per-term influence details.
#' @slot input_hash Digest of the input assay and sample metadata.
#' @name LOOResult-class
#' @rdname LOOResult-class
#' @aliases LOOResult
#' @docType class
#' @return An S4 object of class `LOOResult`, containing the
#'   full-sample model results, leave-one-out estimates,
#'   influence statistics, flags, profiles, and diagnostics.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(8))
#' x <- matrix(rnorm(10 * 8), nrow = 10,
#'             dimnames = list(paste0("gene", 1:10), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
#'                  row.names = samples)
#' result <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' is(result, "LOOResult")
NULL

# LOOResult: result container ----------------------------------------------

setClass("LOOResult",
    slots = c(
        call         = "ANY",
        engine       = "character",
        engine_info  = "list",
        coefficients = "data.frame",  # full-sample fit: feature, term, beta, se, p
        loo_betas    = "matrix",     # rows = feature:term, cols = left-out samples
        loo_se       = "matrix",
        shifts       = "data.frame", # long: feature, term, sample, delta, dfbetas, ...
        influence    = "data.frame", # per-sample aggregate scores
        failures     = "data.frame", # per-sample refit errors (fail-soft)
        flags        = "data.frame", # populated by flag_influential()
        profiles     = "list",      # populated by profile_flags()
        diagnostics  = "list",
        input_hash   = "character"
    ),
    prototype = list(
        call = NULL,
        engine = "lm",
        engine_info = list(),
        coefficients = data.frame(),
        loo_betas = matrix(NA_real_, 0, 0),
        loo_se = matrix(NA_real_, 0, 0),
        shifts = data.frame(),
        influence = data.frame(),
        failures = data.frame(),
        flags = data.frame(),
        profiles = list(),
        diagnostics = list(),
        input_hash = NA_character_
    )
)

setValidity("LOOResult", function(object) {
    msg <- character(0)
    if (ncol(object@loo_betas) != ncol(object@loo_se))
        msg <- c(msg, "loo_betas and loo_se must have identical dimensions")
    if (nrow(object@shifts) > 0L &&
        !all(c("feature", "term", "sample", "delta", "dfbetas") %in% colnames(object@shifts)))
        msg <- c(msg, "shifts is missing required columns")
    if (length(msg)) msg else TRUE
})
