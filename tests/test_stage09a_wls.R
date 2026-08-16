#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "workflow", "scripts", "stage09a_wls.R"))

passed <- 0L
check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  passed <<- passed + 1L
}
close_enough <- function(observed, expected, tolerance = 1e-10) {
  max(abs(observed - expected)) <= tolerance
}

sample <- c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2")
group <- c("s1_v1", "s1_v1", "s1_v2", "s1_v2", "s2_v1", "s2_v1", "s2_v1", "s2_v1")
x <- cbind(
  terminal_C = c(0, 1, 0, 1, 0, 1, 0, 1),
  gc = c(0.1, 0.4, 0.6, 0.2, 0.3, 0.9, 0.5, 0.7)
)
y <- c(0.2, 1.0, 2.3, 2.8, -0.2, 1.4, 0.8, 2.1)
weights <- stage09a_sample_weights(sample, group)

check(abs(mean(weights) - 1) < 1e-14, "weights do not have mean one")
check(abs(sum(weights[sample == "s1"]) - sum(weights[sample == "s2"])) < 1e-14,
      "samples do not contribute equal total weight")
check(abs(sum(weights[group == "s1_v1"]) - sum(weights[group == "s1_v2"])) < 1e-14,
      "groups within sample do not contribute equal total weight")

within_fit <- stage09a_fit_wls(x, y, weights, group)
group_design <- stats::model.matrix(~ 0 + factor(group))
explicit_fit <- stats::lm.wfit(cbind(group_design, x), y, weights, singular.ok = FALSE)
explicit_sequence <- tail(explicit_fit$coefficients, ncol(x))
check(close_enough(within_fit$coefficients, explicit_sequence),
      "within-group WLS is not equivalent to explicit unpenalized nuisance intercepts")

prediction <- stage09a_predict_sequence_only(within_fit, x)
check(close_enough(prediction, as.numeric(x %*% within_fit$coefficients)),
      "candidate prediction contains something other than sequence coefficients")
check(within_fit$rank == ncol(x), "synthetic sequence design is rank deficient")
check(within_fit$n_groups == length(unique(group)), "group accounting is incorrect")

cat(sprintf("%d Stage 09A v0.20 WLS assertions passed\n", passed))
