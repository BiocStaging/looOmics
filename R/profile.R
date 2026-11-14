# Metadata profiling --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Collapse rare factor levels while retaining a usable reference level.
.prepare_covariate <- function(v, name, min_cell_size) {
  if (is.numeric(v)) return(v)
  v <- droplevels(factor(v))
  tab <- table(v)
  keep <- names(tab)[tab >= min_cell_size]
  if (length(keep) < 2L) return(NULL)
  collapsed <- as.character(v)
  rare <- !is.na(collapsed) & !(collapsed %in% keep)
  collapsed[rare] <- "Other"
  out <- droplevels(factor(collapsed))
  if (nlevels(out) < 2L) return(NULL)
  out
}

# Simulated Fisher tests are seeded via withr::with_seed(), which restores
# the caller's RNG state afterwards, so the user's random number
# stream is left untouched.

#' Profile flagged samples against biological and technical covariates
#'
#' Evaluates associations between influential-sample status and sample-level
#' biological or technical metadata. Continuous covariates are assessed with
#' Wilcoxon rank-sum tests; categorical covariates use chi-square tests or
#' Fisher's exact test when expected cell counts are sparse. Univariate
#' p-values are adjusted using Benjamini-Hochberg. An optional multivariable
#' logistic regression estimates adjusted associations across eligible
#' metadata covariates. Covariates used in the differential model are retained
#' and labeled, because their association with influence can itself be
#' biologically informative; interpretation should account for their in-model
#' status. This is exploratory association analysis, not causal inference.
#'
#' @param object A \code{\linkS4class{LOOResult}} with flags computed.
#' @param covariates Character vector of metadata columns to profile. Defaults
#'   to all available sample metadata, including covariates in the differential
#'   model.
#' @param model_covariates Character vector of covariates used in the original
#'   differential model. Defaults to variables inferred from the stored model
#'   formula; these are retained and labeled in the output.
#' @param multivariable Logical; fit an exploratory adjusted logistic model.
#'   Numeric covariates are standardized to one standard deviation in this
#'   model; factor coefficients are relative to their reference level.
#' @param min_cell_size Minimum observations per categorical level before rare
#'   levels are collapsed into an \code{Other} category.
#' @param seed Integer seed for simulated Fisher tests.
#' @param fisher_B Number of Monte Carlo replicates for simulated Fisher tests.
#'
#' @return The input object with a \code{FlagProfile} in its profiles slot.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(
#'   condition = factor(rep(c("A", "B"), each = 6)),
#'   batch = factor(rep(c("batch1", "batch2"), 6)),
#'   age = seq(40, 62, length.out = 12),
#'   row.names = samples
#' )
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' res <- flag_influential(res, method = "topk", k = 3)
#' res <- profile_flags(res, covariates = c("batch", "age"),
#'                      multivariable = FALSE)
#' profiles(res)
#' @export
profile_flags <- function(object, covariates = NULL, model_covariates = NULL,
                          multivariable = TRUE, min_cell_size = 3,
                          seed = 1L, fisher_B = 2000L) {
  if (!is(object, "LOOResult"))
    stop("'object' must be a LOOResult.", call. = FALSE)
  if (nrow(object@flags) == 0L)
    stop("No flags found; run flag_influential() first.", call. = FALSE)
  if (length(min_cell_size) != 1L || is.na(min_cell_size) ||
      min_cell_size < 1 || min_cell_size != as.integer(min_cell_size))
    stop("'min_cell_size' must be a positive integer.", call. = FALSE)
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed))
    stop("'seed' must be one non-missing number.", call. = FALSE)
  if (length(fisher_B) != 1L || is.na(fisher_B) || fisher_B < 1 ||
      fisher_B != as.integer(fisher_B))
    stop("'fisher_B' must be a positive integer.", call. = FALSE)
  
  meta <- object@diagnostics$col_data
  if (is.null(meta))
    stop("Sample metadata are not available on this result object.", call. = FALSE)
  if (is.null(model_covariates))
    model_covariates <- object@diagnostics$model_variables %||% character()
  if (!is.character(model_covariates))
    stop("'model_covariates' must be a character vector.", call. = FALSE)
  if (is.null(covariates)) covariates <- names(meta)
  if (!is.character(covariates))
    stop("'covariates' must be a character vector.", call. = FALSE)
  
  unknown <- setdiff(unique(c(covariates, model_covariates)), names(meta))
  if (length(unknown))
    stop("Unknown metadata covariate(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  covariates <- unique(covariates)
  model_covariates <- unique(model_covariates)
  
  flags <- object@flags
  flag_by_sample <- setNames(flags$flagged, flags$sample)
  common <- intersect(rownames(meta), flags$sample)
  if (length(common) < 4L)
    stop("At least four matched samples are required for metadata profiling.",
         call. = FALSE)
  fl <- as.logical(flag_by_sample[common])
  meta <- meta[common, , drop = FALSE]
  if (anyNA(fl) || length(unique(fl)) < 2L)
    stop("Metadata profiling requires both flagged and unflagged samples.",
         call. = FALSE)
  
  rows <- list()
  dropped <- data.frame(covariate = character(), reason = character(),
                        stringsAsFactors = FALSE)
  .drop <- function(df, cv, reason) {
    rbind(df, data.frame(covariate = cv, reason = reason,
                         stringsAsFactors = FALSE))
  }
  
  for (cv in covariates) {
    v_all <- meta[[cv]]
    ok <- !is.na(v_all)
    if (sum(ok) < max(4L, ceiling(length(fl) * 0.5))) {
      dropped <- .drop(dropped, cv, paste0("too many missing values (", sum(!ok), ")"))
      next
    }
    v <- v_all[ok]
    fl_v <- fl[ok]
    if (length(unique(fl_v)) < 2L) {
      dropped <- .drop(dropped, cv, "complete cases contain only one flag-status group")
      next
    }
    if (is.numeric(v)) {
      if (length(unique(v)) < 2L) {
        dropped <- .drop(dropped, cv, "covariate has no variation")
        next
      }
      wt <- tryCatch(suppressWarnings(stats::wilcox.test(v ~ fl_v)),
                     error = function(e) NULL)
      if (is.null(wt)) {
        dropped <- .drop(dropped, cv, "Wilcoxon test could not be computed")
        next
      }
      eff <- suppressWarnings(stats::cor(v, as.numeric(fl_v),
                                         method = "spearman"))
      rows[[cv]] <- data.frame(
        covariate = cv, type = "continuous", test = "wilcox",
        statistic = unname(wt$statistic), p_value = wt$p.value,
        effect = unname(eff), effect_ci = NA_character_,
        in_model = cv %in% model_covariates,
        n_used = sum(ok), stringsAsFactors = FALSE)
    } else {
      v <- .prepare_covariate(v, cv, min_cell_size)
      if (is.null(v)) {
        dropped <- .drop(dropped, cv, "fewer than two categorical levels remain after rare-level handling")
        next
      }
      tab <- table(v, fl_v)
      if (nrow(tab) < 2L || ncol(tab) < 2L) {
        dropped <- .drop(dropped, cv, "contingency table has fewer than two levels in a dimension")
        next
      }
      expected <- tryCatch(suppressWarnings(stats::chisq.test(tab)$expected),
                           error = function(e) NULL)
      use_fisher <- is.null(expected) || any(expected < 5)
      tst <- if (use_fisher) {
        withr::with_seed(as.integer(seed),
                         suppressWarnings(stats::fisher.test(
                           tab, simulate.p.value = TRUE, B = as.integer(fisher_B))))
      } else suppressWarnings(stats::chisq.test(tab))
      cramers_v <- tryCatch({
        chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
        sqrt(as.numeric(chi$statistic) /
               (sum(tab) * (min(dim(tab)) - 1)))
      }, error = function(e) NA_real_)
      rows[[cv]] <- data.frame(
        covariate = cv, type = "categorical",
        test = if (use_fisher) "fisher_simulated" else "chisq",
        statistic = if (is.null(tst$statistic)) NA_real_
        else unname(as.numeric(tst$statistic)),
        p_value = tst$p.value, effect = unname(cramers_v),
        effect_ci = NA_character_, in_model = cv %in% model_covariates,
        n_used = sum(ok), stringsAsFactors = FALSE)
    }
  }
  
  if (length(rows) == 0L) {
    univariate <- data.frame(
      covariate = character(), type = character(), test = character(),
      statistic = numeric(), p_value = numeric(), effect = numeric(),
      effect_ci = character(), in_model = logical(), n_used = integer(),
      p_adj = numeric(), stringsAsFactors = FALSE)
  } else {
    univariate <- do.call(rbind, rows)
    univariate$p_adj <- stats::p.adjust(univariate$p_value, method = "BH")
    rownames(univariate) <- NULL
  }
  
  # The adjusted model includes both explicitly profiled covariates and
  # original model covariates, preserving biological context. It uses one
  # complete-case cohort so all adjusted coefficients refer to the same set
  # of samples. Rare categorical levels are collapsed consistently.
  mv <- NULL
  mv_note <- NULL
  mv_covariates <- unique(c(covariates, model_covariates))
  if (multivariable) {
    if (length(mv_covariates) == 0L) {
      mv_note <- "no metadata covariates were selected"
    } else if (nrow(univariate) == 0L) {
      # None of the explicitly profiled covariates could be tested
      # (all dropped for missing values or insufficient levels), so
      # there is nothing to adjust: report why instead of fitting a
      # model on model covariates alone.
      mv_note <- paste("none of the profiled covariates could be tested;",
                       "no covariates remained available for",
                       "multivariable profiling")
    } else if (sum(fl) < 5L || sum(!fl) < 5L) {
      mv_note <- "too few flagged or unflagged samples for adjusted logistic regression (at least five in each group are required)"
    } else {
      mv_meta <- meta[, mv_covariates, drop = FALSE]
      prepared <- list()
      for (cv in mv_covariates) {
        v <- mv_meta[[cv]]
        if (is.numeric(v)) {
          if (length(unique(v[!is.na(v)])) < 2L) {
            dropped <- .drop(dropped, cv, "excluded from multivariable model: no variation")
            next
          }
          prepared[[cv]] <- as.numeric(scale(v))
        } else {
          vv <- .prepare_covariate(v, cv, min_cell_size)
          if (is.null(vv)) {
            dropped <- .drop(dropped, cv, "excluded from multivariable model: insufficient categorical levels")
            next
          }
          prepared[[cv]] <- vv
        }
      }
      if (length(prepared) == 0L) {
        mv_note <- "no eligible covariates for the adjusted model"
      } else {
        dat <- data.frame(flagged = as.integer(fl),
                          prepared, check.names = FALSE)
        complete <- stats::complete.cases(dat)
        dat <- dat[complete, , drop = FALSE]
        if (sum(dat$flagged == 1L) < 5L || sum(dat$flagged == 0L) < 5L) {
          mv_note <- "complete-case data leave fewer than five samples in one flag-status group"
        } else {
          dat <- droplevels(dat)
          mm <- tryCatch(stats::model.matrix(flagged ~ ., data = dat),
                         error = function(e) NULL)
          if (is.null(mm)) {
            mv_note <- "adjusted model matrix could not be constructed after complete-case filtering"
          } else if (qr(mm)$rank < ncol(mm)) {
            mv_note <- "adjusted model design is rank deficient; remove redundant or confounded covariates"
          } else if (min(sum(dat$flagged == 1L), sum(dat$flagged == 0L)) <
                     5L * (ncol(mm) - 1L)) {
            mv_note <- "limited events per model parameter; adjusted odds ratios may be unstable"
          } else {
            warn_env <- new.env(parent = emptyenv())
            warn_env$msgs <- character()
            fit <- withCallingHandlers(
              tryCatch(stats::glm(flagged ~ ., data = dat,
                                  family = stats::binomial()),
                       error = function(e) NULL),
              warning = function(w) {
                warn_env$msgs <- c(warn_env$msgs, conditionMessage(w))
                invokeRestart("muffleWarning")
              })
            if (is.null(fit) || !isTRUE(fit$converged)) {
              mv_note <- "adjusted logistic regression failed to converge"
            } else {
              if (length(warn_env$msgs))
                mv_note <- paste(unique(warn_env$msgs), collapse = "; ")
              if (any(fit$fitted.values < 1e-6 |
                      fit$fitted.values > 1 - 1e-6))
                mv_note <- paste(c(mv_note,
                                   "possible separation or near-separation; odds ratios may be unreliable"),
                                 collapse = "; ")
              mv <- tryCatch(
                broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
                error = function(e) broom::tidy(fit,
                                                exponentiate = TRUE, conf.int = FALSE))
            }
          }
        }
      }
    }
  }
  
  prof <- structure(
    list(univariate = univariate, multivariable = mv,
         flag_counts = c(flagged = sum(fl), not_flagged = sum(!fl)),
         covariates = covariates, adjust = "BH",
         model_covariates = model_covariates,
         multivariable_covariates = mv_covariates,
         dropped = dropped, mv_note = mv_note,
         seed = as.integer(seed), fisher_B = as.integer(fisher_B)),
    class = "FlagProfile")
  object@profiles <- list(profile = prof)
  object
}

# FlagProfile S3 methods ----------------------------------------------------

#' Print a biological and technical metadata profile
#'
#' Displays univariate associations between influential-sample status and
#' sample metadata, optional multivariable logistic-regression results, and
#' reasons why selected covariates could not be tested.
#'
#' @param x A \code{FlagProfile} object returned by \code{profiles()}.
#' @param ... Additional arguments (currently unused).
#' @return The input \code{FlagProfile} object, invisibly.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(
#'   condition = factor(rep(c("A", "B"), each = 6)),
#'   batch = factor(rep(c("batch1", "batch2"), 6)),
#'   row.names = samples
#' )
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' res <- flag_influential(res, method = "topk", k = 3)
#' res <- profile_flags(res, multivariable = FALSE)
#' print(profiles(res))
#' @export
print.FlagProfile <- function(x, ...) {
  cat("FlagProfile: influential-sample metadata profiling\n")
  cat("  flagged:", x$flag_counts[["flagged"]],
      "| not flagged:", x$flag_counts[["not_flagged"]], "\n")
  if (nrow(x$univariate) > 0) {
    cat("\nUnivariate tests (BH-adjusted):\n")
    print(x$univariate[, c("covariate", "test", "p_value", "p_adj",
                           "effect", "in_model")], row.names = FALSE)
  } else cat("\nNo eligible covariates.\n")
  if (!is.null(x$mv_note))
    cat("\nMultivariable model note:", x$mv_note, "\n")
  if (!is.null(x$multivariable)) {
    cat("\nMultivariable logistic model (ORs):\n")
    print(x$multivariable[, intersect(c("term", "estimate", "conf.low",
                                        "conf.high", "p.value"),
                                      names(x$multivariable))],
          row.names = FALSE)
  }
  if (nrow(x$dropped) > 0) {
    cat("\nCovariates not tested:\n")
    print(x$dropped, row.names = FALSE)
  }
  invisible(x)
}