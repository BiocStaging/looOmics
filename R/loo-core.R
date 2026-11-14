# Core driver ---------------------------------------------------------------

#' Sample-wise leave-one-out sensitivity analysis for omics models
#'
#' For each sample i, refits the model on the remaining n-1 samples and
#' records coefficient shifts for every feature and model term, then computes
#' influence metrics per sample.
#'
#' The column \code{dfbetas} in the shift table is a deletion-based
#' standardized coefficient shift, (beta_full - beta_loo) / se_loo. For
#' ordinary linear models this coincides with classical DFBETAS; for other
#' engines it is the analogous standardized deletion shift and should be
#' interpreted as such. This is a sensitivity analysis for detecting
#' influential samples, not an automatic bad-sample detector: an influential
#' sample may reflect genuine biological heterogeneity.
#'
#' @param x Omics object: matrix/data.frame (features x samples, with
#'   \code{col_data}), or an \code{OmicsInput}. SummarizedExperiment /
#'   ExpressionSet / phyloseq are supported when those packages are installed.
#' @param formula Model formula evaluated in the sample covariate table;
#'   the left-hand side is implicit (each feature is the response).
#' @param engine Engine id: "auto" (default; resolves to "lm"), "lm",
#'   "glm", "limma", "deseq2", or a name registered
#'   via \code{\link{register_engine}}.
#' @param parallel Logical, default TRUE. If FALSE, runs sequentially via
#'   lapply (useful for debugging). The parallel strategy itself is owned by
#'   the user via \code{future::plan()}.
#' @param seed Integer or TRUE, passed to
#'   \code{future.apply::future_lapply(future.seed =)}.
#' @param ... Passed to the engine (e.g. \code{family} for glm).
#'
#' @return A \code{\linkS4class{LOOResult}} with full-sample coefficients,
#'   LOO betas/SEs, the long shift table, and per-sample influence scores.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(100 * 12), 100, 12,
#'               dimnames = list(paste0("g", 1:100), paste0("S", 1:12)))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), 6)),
#'                  batch = factor(rep(1:2, 6)), row.names = colnames(mat))
#' res <- loo_sensitivity(mat, ~ condition, col_data = cd)
#' flag_influential(res)
#' @export
loo_sensitivity <- function(x, formula, engine = "auto", parallel = TRUE,
                            seed = TRUE, ...) {
    cl <- match.call()
    dots <- list(...)
    adapter_args <- c("col_data", "row_data", "transform", "pseudocount",
                      "assay_name")
    if (is(x, "OmicsInput")) {
        input <- x
        engine_args <- dots
    } else {
        input <- do.call(as_omics_input,
            c(list(x = x, formula = formula, engine = engine), dots))
        # Arguments consumed by the input adapter must not leak into engine
        # calls (e.g. col_data must not be forwarded to DESeq2::DESeq()).
        engine_args <- dots[setdiff(names(dots), adapter_args)]
    }
    engine <- input@engine
    cd <- input@col_data

    # Formula: allow RHS-only, and store the exact one-sided design formula
    # used by formula-driven engines such as DESeq2.
    f <- if (length(formula) == 2L) formula else formula[-2]
    input@formula <- paste(deparse(f), collapse = " ")
    design <- stats::model.matrix(f, data = cd)

    # Full-sample fit
    full <- do.call(run_engine, c(list(input = input, design = design),
                                  engine_args))
    full <- normalize_engine_fit(full, rownames(input@assay), colnames(design))

    n <- ncol(input@assay)
    samples <- colnames(input@assay)

    ft <- full$colnames
    feats <- rownames(input@assay)
    fit_one <- function(i) {
        tryCatch({
            sub <- input
            sub@assay <- input@assay[, -i, drop = FALSE]
            sub@col_data <- input@col_data[-i, , drop = FALSE]
            dsub <- design[-i, , drop = FALSE]
            if (qr(dsub)$rank < ncol(dsub))
                stop("rank-deficient design after dropping sample ", samples[i])
            fit_args <- c(list(input = sub, design = dsub), engine_args)
            normalize_engine_fit(do.call(run_engine, fit_args),
                                 rownames(input@assay), ft)
        }, error = function(e) list(error = conditionMessage(e)))
    }

    loo_fits <- if (parallel) {
        future.apply::future_lapply(seq_len(n), fit_one, future.seed = seed)
    } else {
        lapply(seq_len(n), fit_one)
    }

    failures <- data.frame(
        sample = samples[vapply(loo_fits, function(f) !is.null(f$error), logical(1))],
        error = vapply(loo_fits, function(f) if (is.null(f$error)) "" else f$error,
                       character(1))[vapply(loo_fits, function(f) !is.null(f$error), logical(1))],
        stringsAsFactors = FALSE)

    ok_idx <- which(vapply(loo_fits, function(f) is.null(f$error), logical(1)))

    feats <- rownames(input@assay)
    # rows = feature:term
    key <- expand.grid(feature = feats, term = ft, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)[, c("feature", "term")]
    key_label <- paste(key$feature, key$term, sep = "::")

    loo_betas <- loo_se <- matrix(NA_real_, nrow(key), n)
    dimnames(loo_betas) <- dimnames(loo_se) <- list(key_label, samples)
    for (i in ok_idx) {
        fi <- loo_fits[[i]]
        bm <- fi$betas[feats, ft, drop = FALSE]
        sm <- fi$se[feats, ft, drop = FALSE]
        loo_betas[, i] <- as.vector(bm)
        loo_se[, i] <- as.vector(sm)
    }

    full_betas <- full$betas[feats, ft, drop = FALSE]
    full_se <- full$se[feats, ft, drop = FALSE]
    full_p <- full$p[feats, ft, drop = FALSE]

    coefficients <- data.frame(
        feature = key$feature, term = key$term,
        beta = as.vector(full_betas), se = as.vector(full_se),
        p = as.vector(full_p), stringsAsFactors = FALSE)

    # ---- influence metrics (computed eagerly; see influence.R for helpers)
    shifts <- build_shifts(key, key_label, samples, full_betas, full_p,
                           loo_betas, loo_se, loo_fits, ok_idx, feats, ft,
                           q_threshold = 0.05)
    influence <- build_influence(shifts, samples)
    influence_by_term <- build_influence_by_term(shifts, samples)

    methods::new("LOOResult",
        call = cl, engine = engine,
        engine_info = list(engine = engine, terms = ft),
        coefficients = coefficients,
        loo_betas = loo_betas, loo_se = loo_se,
        shifts = shifts, influence = influence, failures = failures,
        diagnostics = list(
            n_samples = n, n_features = length(feats), terms = ft,
            parallel = parallel, seed = seed,
            future_plan = paste(class(future::plan()), collapse = "/"),
            runtime = Sys.time(), col_data = input@col_data,
            formula = deparse(f), model_variables = all.vars(f),
            influence_by_term = influence_by_term),
        input_hash = digest::digest(list(input@assay, input@col_data)))
}


# Ensure engine outputs carry feature/term dimnames (custom engines may omit them)
normalize_engine_fit <- function(fit, feats, ft) {
    if (is.null(fit$error)) {
        for (m in c("betas", "se", "p")) {
            mat <- fit[[m]]
            if (is.null(rownames(mat))) rownames(mat) <- feats
            if (is.null(colnames(mat))) colnames(mat) <- ft
            fit[[m]] <- mat
        }
        if (is.null(fit$colnames)) fit$colnames <- ft
    }
    fit
}
