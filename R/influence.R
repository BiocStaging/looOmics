# Influence metrics ---------------------------------------------------------

build_shifts <- function(key, key_label, samples, full_betas, full_p,
                         loo_betas, loo_se, loo_fits, ok_idx, feats, ft,
                         q_threshold) {
    g <- nrow(key)
    n <- length(samples)

    # Adjust p-values and rank features independently within each model term.
    # Lower ranks correspond to smaller p-values (stronger evidence).
    full_q <- matrix(unlist(lapply(seq_len(ncol(full_p)), function(j)
        stats::p.adjust(full_p[, j], method = "BH"))),
        nrow = nrow(full_p), ncol = ncol(full_p))
    full_rank <- matrix(unlist(lapply(seq_len(ncol(full_p)), function(j)
        rank(full_p[, j], ties.method = "min", na.last = "keep"))),
        nrow = nrow(full_p), ncol = ncol(full_p))

    long <- data.frame(
        feature = rep(key$feature, times = n),
        term = rep(key$term, times = n),
        sample = rep(samples, each = g),
        beta = as.vector(full_betas),
        se = NA_real_, delta = NA_real_, dfbetas = NA_real_,
        p_full = as.vector(full_p), q_full = as.vector(full_q),
        rank_full = as.vector(full_rank),
        p_loo = NA_real_, q_loo = NA_real_, rank_loo = NA_real_,
        rank_change = NA_real_, sig_flip = NA,
        error = FALSE, stringsAsFactors = FALSE)

    for (i in ok_idx) {
        idx <- seq_len(g) + (i - 1L) * g
        loo_b <- loo_betas[, i]
        loo_s <- loo_se[, i]
        fi <- loo_fits[[i]]
        p_loo <- fi$p[feats, ft, drop = FALSE]
        q_loo <- matrix(unlist(lapply(seq_len(ncol(p_loo)), function(j)
            stats::p.adjust(p_loo[, j], method = "BH"))),
            nrow = nrow(p_loo), ncol = ncol(p_loo))
        rank_loo <- matrix(unlist(lapply(seq_len(ncol(p_loo)), function(j)
            rank(p_loo[, j], ties.method = "min", na.last = "keep"))),
            nrow = nrow(p_loo), ncol = ncol(p_loo))
        long$delta[idx] <- long$beta[idx] - loo_b
        long$dfbetas[idx] <- (long$beta[idx] - loo_b) / loo_s
        long$p_loo[idx] <- as.vector(p_loo)
        long$q_loo[idx] <- as.vector(q_loo)
        long$rank_loo[idx] <- as.vector(rank_loo)
        long$rank_change[idx] <- as.vector(rank_loo - full_rank)
        long$se[idx] <- loo_s
        long$sig_flip[idx] <- (as.vector(full_q) < q_threshold) !=
            (as.vector(q_loo) < q_threshold)
    }
    err_idx <- setdiff(seq_len(n), ok_idx)
    if (length(err_idx)) long$error[long$sample %in% samples[err_idx]] <- TRUE
    long
}

build_influence <- function(shifts, samples) {
    ok <- !shifts$error
    inf <- lapply(samples, function(s) {
        sub <- shifts[shifts$sample == s & ok, , drop = FALSE]
        if (nrow(sub) == 0)
            return(c(max_dfbetas = NA_real_, median_abs_delta = NA_real_,
                     n_sig_flips = NA_integer_, n_features_shifted = NA_integer_))
        .agg_row(sub)
    })
    out <- as.data.frame(do.call(rbind, inf), stringsAsFactors = FALSE)
    out$sample <- samples
    out$error <- samples %in% shifts$sample[shifts$error]
    out[!out$error, "n_sig_flips"] <- as.integer(out[!out$error, "n_sig_flips"])
    out[!out$error, "n_features_shifted"] <- as.integer(out[!out$error, "n_features_shifted"])
    out[, c("sample", "max_dfbetas", "median_abs_delta", "n_sig_flips",
            "n_features_shifted", "error")]
}

.agg_row <- function(sub) {
    d <- abs(sub$dfbetas)
    c(max_dfbetas = if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE),
      median_abs_delta = if (all(is.na(sub$delta))) NA_real_
                         else stats::median(abs(sub$delta), na.rm = TRUE),
      n_sig_flips = sum(sub$sig_flip, na.rm = TRUE),
      n_features_shifted = length(unique(sub$feature[d > 1])))
}

# Per-term influence: same aggregates as build_influence() but computed
# separately for each model term. Defaults to all terms (matching the pooled
# influence slot); pass terms = setdiff(terms, "(Intercept)") to exclude the
# intercept.
build_influence_by_term <- function(shifts, samples, terms = NULL) {
    ok <- shifts[!shifts$error, , drop = FALSE]
    if (is.null(terms)) terms <- unique(ok$term)
    out <- data.frame(stringsAsFactors = FALSE)
    for (tm in terms) {
        sub <- ok[ok$term == tm, , drop = FALSE]
        if (nrow(sub) == 0) next
        for (s in samples) {
            ss <- sub[sub$sample == s, , drop = FALSE]
            if (nrow(ss) == 0) next
            out <- rbind(out, data.frame(
                term = tm, sample = s,
                max_abs_shift = if (all(is.na(ss$dfbetas))) NA_real_
                                else max(abs(ss$dfbetas), na.rm = TRUE),
                n_sig_flips = sum(ss$sig_flip, na.rm = TRUE),
                stringsAsFactors = FALSE))
        }
    }
    out
}
