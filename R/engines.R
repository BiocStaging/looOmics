# Engine dispatch -----------------------------------------------------------

#' Fit a model using the engine stored in an OmicsInput
#'
#' Dispatches to the built-in or user-registered fitting engine selected by
#' \code{input@engine}. Most users should call \code{loo_sensitivity()},
#' which invokes this generic for the full-sample and leave-one-out fits.
#'
#' @param input An \code{\linkS4class{OmicsInput}}.
#' @param design Numeric model matrix with samples in rows.
#' @param contrast Optional contrast specification passed to the engine.
#' @param block Optional blocking factor passed to the engine.
#' @param ... Additional engine-specific arguments.
#'
#' @return A list containing feature-by-term matrices named \code{betas},
#'   \code{se}, and \code{p}, plus \code{colnames} for the model terms.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(8))
#' x <- matrix(rnorm(10 * 8), nrow = 10,
#'             dimnames = list(paste0("gene", 1:10), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
#'                  row.names = samples)
#' input <- as_omics_input(x, ~ condition, col_data = cd)
#' design <- stats::model.matrix(~ condition, data = cd)
#' fit <- run_engine(input, design)
#' dim(fit$betas)
setGeneric("run_engine",
           function(input, design, contrast = NULL, block = NULL, ...)
             standardGeneric("run_engine"))

setMethod("run_engine", c(input = "OmicsInput"),
          function(input, design, contrast = NULL, block = NULL, ...) {
            name <- input@engine
            fun <- paste0("run_engine_", name)
            if (exists(name, envir = .engine_registry, inherits = FALSE)) {
              return(.engine_registry[[name]](input, design, contrast, block, ...))
            }
            if (!exists(fun, where = asNamespace("looOmics"), inherits = FALSE))
              stop("Unknown engine '", name, "'. See ?register_engine.", call. = FALSE)
            do.call(fun, list(input, design, contrast, block, ...))
          })

# Engines -------------------------------------------------------------------
# Every engine returns list(betas, se, p) as feature x term matrices with
# rownames = features, colnames = model terms.

# --- lm engine: per-feature OLS via .lm.fit (fast, no formula overhead) ----

run_engine_lm <- function(input, design, contrast = NULL, block = NULL, ...) {
  y <- t(input@assay)  # samples x features
  p <- ncol(design)
  n <- nrow(design)
  XtX <- crossprod(design)
  if (qr(design)$rank < p) stop("Design matrix is rank deficient.", call. = FALSE)
  XtX_inv <- solve(XtX)
  # betas = (X'X)^{-1} X'y ; crossprod(design, y) is p x g
  betas <- t(XtX_inv %*% crossprod(design, y))   # features x p
  resid <- y - design %*% t(betas)
  dof <- n - p
  if (dof <= 0L)
    stop("The design has no residual degrees of freedom; reduce the number of model terms or provide more samples.",
         call. = FALSE)
  sigma2 <- colSums(resid^2) / dof               # per-feature
  vif <- diag(XtX_inv)                            # diag of (X'X)^{-1}
  se <- sqrt(outer(sigma2, vif))                  # features x p
  # t statistics
  tstat <- betas / se
  pvals <- 2 * pt(-abs(tstat), df = dof)
  dimnames(betas) <- dimnames(se) <- dimnames(pvals) <-
    list(rownames(input@assay), colnames(design))
  list(betas = betas, se = se, p = pvals,
       colnames = colnames(design))
}

# --- glm engine: per-feature glm via stats::glm (binomial default) ---------

run_engine_glm <- function(input, design, contrast = NULL, block = NULL,
                           family = gaussian(), ...) {
  fam <- family
  fam_name <- fam$family
  estimate_disp <- fam_name == "gaussian" || grepl("^quasi", fam_name)
  y <- t(input@assay)
  g <- ncol(y)
  p <- ncol(design)
  betas <- se <- pvals <- matrix(NA_real_, g, p,
                                 dimnames = list(rownames(input@assay),
                                                 colnames(design)))
  for (j in seq_len(g)) {
    fit <- tryCatch(
      stats::glm.fit(x = design, y = y[, j], family = fam),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged) ||
        is.null(fit$qr) || fit$rank < p || fit$df.residual <= 0L ||
        anyNA(fit$coefficients)) next
    betas[j, ] <- fit$coefficients
    # cov.unscaled = (X'WX)^{-1} from the final IRLS QR. Restore the
    # original coefficient order after QR pivoting, as summary.glm does.
    cov_unsc <- tryCatch({
      rmat <- qr.R(fit$qr)[seq_len(p), seq_len(p), drop = FALSE]
      cov_pivot <- chol2inv(rmat)
      cov <- matrix(0, p, p)
      piv <- fit$qr$pivot[seq_len(p)]
      cov[piv, piv] <- cov_pivot
      cov
    }, error = function(e) NULL)
    if (is.null(cov_unsc) || any(!is.finite(diag(cov_unsc)))) next
    # Families with a free dispersion parameter (gaussian, quasi-*)
    # estimate it from the weighted residual sum of squares and use
    # t-based inference; binomial / Poisson fix dispersion at 1 and use
    # z-based inference.
    disp <- if (estimate_disp)
      sum(fit$weights * fit$residuals^2) / fit$df.residual else 1
    se[j, ] <- sqrt(diag(cov_unsc) * disp)
    tstat <- betas[j, ] / se[j, ]
    pvals[j, ] <- if (estimate_disp)
      2 * stats::pt(-abs(tstat), df = fit$df.residual) else
        2 * stats::pnorm(-abs(tstat))
  }
  list(betas = betas, se = se, p = pvals, colnames = colnames(design))
}

# --- Bioconductor engines (require Suggests packages) ----------------------

run_engine_limma <- function(input, design, contrast = NULL, block = NULL, ...) {
  requireNamespace_or_stop("limma", "the limma engine")
  fit <- limma::lmFit(input@assay, design)
  if (!is.null(contrast)) {
    con <- limma::makeContrasts(contrasts = contrast, levels = design)
    fit <- limma::contrasts.fit(fit, con)
  }
  fit <- limma::eBayes(fit, ...)
  list(betas = fit$coefficients,
       se = sqrt(fit$s2.post) * fit$stdev.unscaled,
       p = fit$p.value,
       colnames = colnames(fit$coefficients))
}

# DESeq2 is formula-driven: the engine uses the model formula stored on the
# OmicsInput (set by as_omics_input) rather than the precomputed design
# matrix, so DESeq2 receives raw integer counts and the original metadata.
# Every full-sample and leave-one-out call constructs a fresh DESeqDataSet,
# so size factors and dispersions are re-estimated for every fit.
run_engine_deseq2 <- function(input, design, contrast = NULL, block = NULL, ...) {
  requireNamespace_or_stop("DESeq2", "the DESeq2 engine")
  if (!identical(input@transform, "none"))
    stop("DESeq2 requires untransformed raw counts; use transform = \"none\".",
         call. = FALSE)
  if (any(!is.finite(input@assay)) || any(input@assay < 0) ||
      any(abs(input@assay - round(input@assay)) > 1e-8))
    stop("DESeq2 requires finite, non-negative integer raw counts; values are never rounded automatically.",
         call. = FALSE)
  if (any(input@assay > .Machine$integer.max))
    stop("DESeq2 count values exceed R's supported integer range.",
         call. = FALSE)
  counts <- input@assay
  storage.mode(counts) <- "integer"
  fml_txt <- input@formula
  if (is.null(fml_txt) || !nzchar(fml_txt))
    stop("DESeq2 engine requires the model formula stored on the ",
         "OmicsInput; construct inputs with as_omics_input().",
         call. = FALSE)
  if (!identical(rownames(input@col_data), colnames(counts)))
    stop("Sample metadata rows must be in the same order as count-matrix columns.",
         call. = FALSE)
  cd <- input@col_data
  for (nm in names(cd)) {
    if (is.character(cd[[nm]])) cd[[nm]] <- factor(cd[[nm]])
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts, colData = cd,
    design = stats::as.formula(fml_txt))
  dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)
  
  feats <- rownames(input@assay)
  .mat <- function(v, nm) matrix(v, nrow = length(feats),
                                 dimnames = list(feats, nm))
  
  if (!is.null(contrast)) {
    res <- DESeq2::results(dds, contrast = contrast)
    out <- list(betas = .mat(res$log2FoldChange, "contrast"),
                se = .mat(res$lfcSE, "contrast"),
                p = .mat(res$pvalue, "contrast"))
    out$colnames <- "contrast"
    return(out)
  }
  
  # Term-wise results: map model-matrix column names (e.g. "conditionB")
  # to DESeq2 result names (e.g. "condition_B_vs_A") or contrast vectors.
  terms <- colnames(design)[colnames(design) != "(Intercept)"]
  resn <- DESeq2::resultsNames(dds)
  betas <- se <- p <- matrix(NA_real_, length(feats), length(terms),
                             dimnames = list(feats, terms))
  for (tm in terms) {
    nm <- .deseq_result_name(tm, input@col_data, resn)
    got <- if (length(nm) == 1L && nm %in% resn) {
      DESeq2::results(dds, name = nm)
    } else {
      con <- .deseq_contrast(tm, input@col_data)
      if (is.null(con))
        stop("Could not map model-matrix term '", tm,
             "' to a DESeq2 coefficient or contrast. Use a simple main-effect design or supply an explicit contrast.",
             call. = FALSE)
      DESeq2::results(dds, contrast = con)
    }
    betas[, tm] <- got$log2FoldChange
    se[, tm] <- got$lfcSE
    p[, tm] <- got$pvalue
  }
  list(betas = betas, se = se, p = p, colnames = terms)
}

# Map a model.matrix term like "conditionB" to a DESeq2 resultsName like
# "condition_B_vs_A" (factor terms) or the variable name (numeric terms).
.deseq_result_name <- function(term, col_data, resn) {
  if (term %in% names(col_data) && is.numeric(col_data[[term]]))
    return(term)
  for (v in names(col_data)) {
    vals <- col_data[[v]]
    if (is.factor(vals) || is.character(vals)) {
      levs <- levels(factor(vals))
      for (l in levs) {
        if (term == paste0(v, l)) {
          # DESeq2 names factor coefficients
          # "<variable>_<level>_vs_<reference>", where the
          # reference is the FIRST factor level.
          ref <- levs[1L]
          cand <- paste0(v, "_", l, "_vs_", ref)
          if (!identical(l, ref) && cand %in% resn) return(cand)
        }
      }
    }
  }
  term
}

# Fallback contrast vector c(variable, numerator, denominator)
.deseq_contrast <- function(term, col_data) {
  for (v in names(col_data)) {
    vals <- col_data[[v]]
    if (is.factor(vals) || is.character(vals)) {
      levs <- levels(factor(vals))
      for (l in levs) {
        if (term == paste0(v, l)) {
          ref <- levs[1L]
          if (!identical(l, ref)) return(c(v, l, ref))
        }
      }
    }
  }
  NULL
}

# Engine registry for custom engines ----------------------------------------

.engine_registry <- new.env(parent = emptyenv())

#' Register a custom fitting engine
#'
#' @param name Engine id used in \code{loo_sensitivity(engine = ...)}.
#' @param fit_fun function(input, design, contrast, block, ...) returning
#'   \code{list(betas, se, p)} matrices (features x terms).
#' @return Invisibly returns the registered engine name.
#' @examples
#' example_fit <- function(input, design, contrast = NULL, block = NULL, ...) {
#'   features <- rownames(input@assay)
#'   terms <- colnames(design)
#'   estimates <- matrix(0, length(features), length(terms),
#'                       dimnames = list(features, terms))
#'   list(betas = estimates, se = estimates, p = estimates,
#'        colnames = terms)
#' }
#' register_engine("example_engine", example_fit)
#' rm(list = "example_engine",
#'    envir = getFromNamespace(".engine_registry", "looOmics"))
#' @export
register_engine <- function(name, fit_fun) {
  if (!is.function(fit_fun))
    stop("'fit_fun' must be a function.", call. = FALSE)
  .engine_registry[[name]] <- fit_fun
  invisible(name)
}