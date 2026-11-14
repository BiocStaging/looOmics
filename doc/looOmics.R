## ----knitr-setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5
)

## ----setup--------------------------------------------------------------------
library(looOmics)

## ----sim----------------------------------------------------------------------
set.seed(42)
n_feat <- 200
samples <- sprintf("S%02d", 1:16)
condition <- factor(rep(c("A", "B"), each = 8))
batch <- factor(rep(c("kit1", "kit2"), times = 8))

assay <- matrix(rnorm(n_feat * 16), n_feat, 16,
                dimnames = list(paste0("g", 1:n_feat), samples))
assay[1:30, condition == "B"] <- assay[1:30, condition == "B"] + 1.2

# planted influential sample: S09 (group B) has extreme values on g1-g20
assay[1:20, "S09"] <- assay[1:20, "S09"] + 6

col_data <- data.frame(
    condition = condition, batch = batch,
    depth = round(rnorm(16, mean = 5e6, sd = 4e5)),
    extraction = factor(rep(c("good", "good", "good", "bad"), times = 4)),
    row.names = samples
)

## ----loo----------------------------------------------------------------------
res <- loo_sensitivity(assay, ~ condition + batch, col_data = col_data)
res

## ----summarized-experiment----------------------------------------------------
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = assay),
    colData = S4Vectors::DataFrame(col_data)
)
res_se <- loo_sensitivity(
    se, ~ condition + batch, engine = "lm",
    assay_name = "logcounts", parallel = FALSE
)
res_se

## ----deseq2-summarized-experiment---------------------------------------------
set.seed(21)
count_samples <- sprintf("C%02d", 1:12)
count_condition <- factor(rep(c("control", "treated"), each = 6),
                          levels = c("control", "treated"))
count_meta <- data.frame(condition = count_condition,
                         row.names = count_samples)
count_matrix <- matrix(
    rnbinom(100 * 12, mu = 100, size = 10), nrow = 100,
    dimnames = list(paste0("gene", 1:100), count_samples)
)
count_matrix[1:20, count_condition == "treated"] <-
    count_matrix[1:20, count_condition == "treated"] + 80L
se_counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = count_matrix),
    colData = S4Vectors::DataFrame(count_meta)
)
res_deseq2 <- loo_sensitivity(
    se_counts, ~ condition, engine = "deseq2",
    assay_name = "counts", parallel = FALSE, fitType = "mean"
)
res_deseq2

## ----shifts-------------------------------------------------------------------
head(tidy(res, what = "shifts"))

## ----flag---------------------------------------------------------------------
res <- flag_influential(res, method = "topk", k = 2)
res

## ----topflags-----------------------------------------------------------------
flags <- tidy(res, what = "flags")
flags[order(-flags$score), ]

## ----byterm-------------------------------------------------------------------
head(influence_by_term(res))

## ----profile------------------------------------------------------------------
res <- profile_flags(res, model_covariates = c("condition", "batch"))
profiles(res)

## ----plot, fig.alt = "Bar plot of per-sample influence scores; the planted sample S09 stands out."----
plot(res)

## ----parallel, eval = FALSE---------------------------------------------------
# future::plan(future::multisession, workers = 4)
# res <- loo_sensitivity(se, ~ condition + batch, engine = "lm")
# future::plan(future::sequential)

## ----session------------------------------------------------------------------
sessionInfo()

