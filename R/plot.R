# Plotting ------------------------------------------------------------------

#' Plot a LOOResult
#'
#' @param object A \code{\linkS4class{LOOResult}}.
#' @param type One of "shift" (per-sample score barplot, flagged in red) or
#'   "scree" (ranked influence scores).
#' @param ... Passed to ggplot layers.
#' @return A ggplot object.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' res <- flag_influential(res, method = "topk", k = 2)
#' autoplot(res, type = "shift")
#' @export
autoplot.LOOResult <- function(object, type = c("shift", "scree"), ...) {
    type <- match.arg(type)
    inf <- object@influence[!object@influence$error, , drop = FALSE]
    if (nrow(inf) == 0) stop("No influence metrics to plot.", call. = FALSE)
    fl <- object@flags
    flagged <- if (nrow(fl) > 0) fl$sample[fl$flagged] else character(0)

    if (type == "shift") {
        inf$flag <- inf$sample %in% flagged
        ggplot2::ggplot(inf, ggplot2::aes(x = stats::reorder(.data$sample,
                                                              -.data$max_dfbetas),
                                          y = .data$max_dfbetas,
                                          fill = .data$flag)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = c(`FALSE` = "grey60",
                                                  `TRUE` = "firebrick"),
                                       guide = "none") +
            ggplot2::labs(x = "Left-out sample", y = "max |DFBETAS|",
                          title = "Per-sample influence score") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                hjust = 1))
    } else {
        s <- sort(inf$max_dfbetas, decreasing = TRUE)
        ggplot2::ggplot(data.frame(rank = seq_along(s), score = s),
                        ggplot2::aes(x = .data$rank, y = .data$score)) +
            ggplot2::geom_point() +
            ggplot2::labs(x = "Rank", y = "max |DFBETAS|",
                          title = "Ranked influence scores")
    }
}

#' @export
setMethod("plot", c(x = "LOOResult", y = "missing"),
          function(x, y, ...) autoplot(x, ...))
