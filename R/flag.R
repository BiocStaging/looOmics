# Flagging ------------------------------------------------------------------

#' Flag influential samples
#'
#' @param object A \code{\linkS4class{LOOResult}} with influence metrics.
#' @param method One of "mad" (robust z on max |DFBETAS|), "threshold"
#'   (absolute cutoff, default 1), "topk", or "quantile".
#' @param threshold Numeric; MAD cutoff (default 5) or absolute cutoff
#'   (default 1) depending on method.
#' @param k Number of top samples for "topk".
#' @param q Quantile for "quantile" (default 0.9).
#' @param metric Aggregate to flag on: "max_dfbetas" (default),
#'   "n_sig_flips", or "median_abs_delta".
#'
#' @return The input object with the populated \code{flags} slot
#'   (columns: sample, score, metric, method, flagged).
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' res <- flag_influential(res, method = "topk", k = 2)
#' res@flags
#' @export
flag_influential <- function(object, method = c("mad", "threshold", "topk",
                                                "quantile"),
                             threshold = NULL, k = NULL, q = 0.9,
                             metric = "max_dfbetas") {
    method <- match.arg(method)
    if (!is(object, "LOOResult"))
        stop("'object' must be a LOOResult.", call. = FALSE)
    if (nrow(object@influence) == 0)
        stop("No influence metrics found; run loo_sensitivity() first.", call. = FALSE)

    inf <- object@influence[!object@influence$error, , drop = FALSE]
    if (!metric %in% colnames(inf))
        stop("Unknown metric '", metric, "'.", call. = FALSE)
    score <- inf[[metric]]
    valid <- !is.na(score) & !is.nan(score)
    if (!any(valid))
        stop("The selected influence metric has no defined values.", call. = FALSE)

    finite_score <- score[is.finite(score)]
    flagged <- rep(FALSE, length(score))
    flagged[valid] <- switch(method,
        mad = {
            th <- if (is.null(threshold)) 5 else threshold
            if (length(th) != 1L || is.na(th) || th < 0)
                stop("'threshold' must be a non-negative number.", call. = FALSE)
            if (length(finite_score) == 0L) {
                rep(TRUE, sum(valid))
            } else {
                med <- stats::median(finite_score)
                madv <- stats::mad(finite_score, constant = 1.4826)
                if (is.na(madv) || madv == 0) {
                    sdv <- stats::sd(finite_score)
                    if (is.na(sdv) || sdv == 0) !is.finite(score[valid])
                    else abs(score[valid] - med) > th * sdv
                } else abs(score[valid] - med) / madv > th
            }
        },
        threshold = {
            th <- if (is.null(threshold)) 1 else threshold
            if (length(th) != 1L || is.na(th) || !is.finite(th))
                stop("'threshold' must be one finite number.", call. = FALSE)
            score[valid] > th
        },
        topk = {
            if (is.null(k)) k <- max(1L, ceiling(0.05 * sum(valid)))
            if (length(k) != 1L || is.na(k) || k < 1L || k != as.integer(k))
                stop("'k' must be a positive integer.", call. = FALSE)
            rank(-score[valid], ties.method = "first") <= min(k, sum(valid))
        },
        quantile = {
            if (length(q) != 1L || is.na(q) || q < 0 || q > 1)
                stop("'q' must be between 0 and 1.", call. = FALSE)
            score[valid] > stats::quantile(score[valid], probs = q, na.rm = TRUE)
        }
    )

    object@flags <- data.frame(
        sample = inf$sample, score = score, metric = metric,
        method = method, flagged = as.logical(flagged),
        stringsAsFactors = FALSE)
    validObject(object)
    object
}
