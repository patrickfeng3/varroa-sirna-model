#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "workflow", "scripts", "stage09a_solver.R"))

passed <- 0L
check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  passed <<- passed + 1L
}
close_enough <- function(observed, expected, tolerance = 1e-7) {
  max(abs(observed - expected)) <= tolerance
}
soft_threshold <- function(value, threshold) {
  sign(value) * max(abs(value) - threshold, 0)
}

check(as.character(getRversion()) == "4.4.3", "unexpected R version")
check(as.character(packageVersion("glmnet")) == "4.1.10", "unexpected glmnet version")

# Gaussian solver regression: the two raw feature columns are weighted-
# orthogonal, so the specification's elastic-net solution has an independent
# closed form for each coefficient.
x <- cbind(
  x1 = c(1, -1, 0, 0),
  x2 = c(2, 1, 3, -2)
)
weights <- c(1, 2, 3, 4)
y <- c(4, 0, -1, 9)
spec_alpha <- 0.3
spec_l1_ratio <- 0.25
weight_sum <- sum(weights)

weighted_cross <- as.numeric(crossprod(x[, 1] * weights, x[, 2]))
check(abs(weighted_cross) < 1e-14, "synthetic Gaussian columns are not weighted-orthogonal")

weighted_a <- colSums(weights * x^2) / weight_sum
weighted_b <- colSums(weights * x * y) / weight_sum
expected_beta <- vapply(
  seq_along(weighted_a),
  function(index) {
    soft_threshold(weighted_b[[index]], spec_alpha * spec_l1_ratio) /
      (weighted_a[[index]] + spec_alpha * (1 - spec_l1_ratio))
  },
  numeric(1)
)

gaussian <- stage09a_fit_accumulation(
  x, y, weights, spec_alpha = spec_alpha, spec_l1_ratio = spec_l1_ratio
)
observed_beta <- as.numeric(gaussian$fit$beta[, 1])

check(gaussian$mapping$glmnet_lambda == spec_alpha, "Gaussian spec alpha did not map to glmnet lambda")
check(gaussian$mapping$glmnet_alpha == spec_l1_ratio, "Gaussian spec l1_ratio did not map to glmnet alpha")
check(close_enough(observed_beta, expected_beta), "Gaussian glmnet coefficients do not match the specified weighted objective")

# Unequal supplied weights, disabled internal x standardization, and disabled
# intercept each materially affect this fixture. These negative controls prove
# the canonical wrapper is not silently using glmnet defaults.
equal_weight_fit <- stage09a_fit_accumulation(
  x, y, rep(1, length(y)), spec_alpha, spec_l1_ratio
)$fit
standardized_fit <- glmnet::glmnet(
  x, y, family = stats::gaussian(), weights = weights,
  lambda = spec_alpha, alpha = spec_l1_ratio,
  standardize = TRUE, intercept = FALSE,
  thresh = 1e-14, maxit = 10000000L
)
intercept_fit <- glmnet::glmnet(
  x, y, family = stats::gaussian(), weights = weights,
  lambda = spec_alpha, alpha = spec_l1_ratio,
  standardize = FALSE, intercept = TRUE,
  thresh = 1e-14, maxit = 10000000L
)
check(max(abs(observed_beta - as.numeric(equal_weight_fit$beta[, 1]))) > 1e-3, "Gaussian supplied weights were not distinguishable")
check(max(abs(observed_beta - as.numeric(standardized_fit$beta[, 1]))) > 1e-3, "Gaussian x standardization negative control did not differ")
check(max(abs(observed_beta - as.numeric(intercept_fit$beta[, 1]))) > 1e-3, "Gaussian intercept negative control did not differ")
check(abs(as.numeric(gaussian$fit$a0[[1]])) < 1e-14, "Gaussian model fitted an intercept")

native_scaled_fit <- glmnet::glmnet(
  x, y, family = "gaussian", weights = weights,
  lambda = spec_alpha, alpha = spec_l1_ratio,
  standardize = FALSE, intercept = FALSE,
  thresh = 1e-14, maxit = 10000000L
)
check(
  max(abs(observed_beta - as.numeric(native_scaled_fit$beta[, 1]))) > 1e-3,
  "native Gaussian response-scaling guard did not distinguish the canonical family-object path"
)

# Representation solver regression. Two complete group-indicator columns are
# unpenalized nuisance intercepts and one sequence column is penalized.
group <- rep(c("g1", "g2"), each = 8)
sequence_x <- rep(c(-2, -1, 0, 1, 2, 3, -3, 0.5), 2)
represented <- c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0)
rep_weights <- c(1, 2, 1, 3, 2, 1, 2, 1, 2, 1, 3, 1, 2, 2, 1, 3)
group_design <- model.matrix(~ 0 + group)
design <- cbind(group_design, sequence_effect = sequence_x)
penalty_factor <- c(0, 0, 1)
rep_spec_alpha <- 0.12
rep_spec_l1_ratio <- 0.25

representation <- stage09a_fit_representation(
  design, represented, rep_weights, penalty_factor,
  spec_alpha = rep_spec_alpha, spec_l1_ratio = rep_spec_l1_ratio
)

check(representation$mapping$p_total == 3, "incorrect total representation-column count")
check(representation$mapping$p_penalized == 1, "incorrect penalized representation-column count")
check(representation$mapping$penalty_factor_internal_scale == 3, "incorrect glmnet penalty.factor internal scale")
check(representation$mapping$glmnet_lambda == rep_spec_alpha / 3, "representation lambda correction is incorrect")
check(representation$mapping$glmnet_alpha == rep_spec_l1_ratio, "representation l1_ratio mapping is incorrect")
check(abs(as.numeric(representation$fit$a0[[1]])) < 1e-14, "representation model fitted an extra intercept")

# Independently minimize the exact project-facing logistic objective with base
# R optim. This uses the specification's unscaled sequence penalty, not
# glmnet's internally rescaled penalty.factor representation.
log1pexp <- function(value) {
  ifelse(value > 0, value + log1p(exp(-value)), log1p(exp(value)))
}
representation_objective <- function(theta) {
  eta <- as.numeric(design %*% theta)
  weighted_loss <- sum(rep_weights * (log1pexp(eta) - represented * eta)) / sum(rep_weights)
  sequence_beta <- theta[[3]]
  penalty <- rep_spec_alpha * (
    rep_spec_l1_ratio * abs(sequence_beta) +
      0.5 * (1 - rep_spec_l1_ratio) * sequence_beta^2
  )
  weighted_loss + penalty
}
independent <- optim(
  par = c(-0.5, 0.2, 0.5), fn = representation_objective,
  method = "BFGS", control = list(reltol = 1e-14, maxit = 100000L)
)
glmnet_theta <- as.numeric(representation$fit$beta[, 1])
check(independent$convergence == 0, "independent logistic objective minimization did not converge")
check(close_enough(glmnet_theta, independent$par, tolerance = 2e-6), "corrected glmnet logistic fit does not match the specified objective")

uncorrected <- glmnet::glmnet(
  design, represented, family = "binomial", weights = rep_weights,
  lambda = rep_spec_alpha, alpha = rep_spec_l1_ratio,
  penalty.factor = penalty_factor, standardize = FALSE, intercept = FALSE,
  thresh = 1e-14, maxit = 10000000L
)
check(
  abs(as.numeric(uncorrected$beta["sequence_effect", 1]) - glmnet_theta[[3]]) > 1e-3,
  "uncorrected penalty.factor lambda unexpectedly matched the corrected fit"
)

cat(sprintf("%d Stage 09A solver-equivalence assertions passed\n", passed))
