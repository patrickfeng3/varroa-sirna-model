# Canonical v0.20 weighted fixed-effect linear regression helper.

stage09a_sample_weights <- function(sample, group) {
  stopifnot(length(sample) == length(group), length(sample) > 0L)
  group_sizes <- table(group)
  groups_per_sample <- vapply(split(group, sample), function(value) length(unique(value)), integer(1))
  raw <- 1 / (as.numeric(groups_per_sample[sample]) * as.numeric(group_sizes[group]))
  raw * length(raw) / sum(raw)
}


stage09a_within_group_transform <- function(x, y, weights, group) {
  stopifnot(
    is.matrix(x), nrow(x) == length(y), length(y) == length(weights),
    length(weights) == length(group), all(is.finite(x)), all(is.finite(y)),
    all(is.finite(weights)), all(weights > 0)
  )
  group_factor <- factor(group, levels = unique(group))
  denominator <- as.numeric(rowsum(weights, group_factor, reorder = FALSE))
  weighted_y <- as.numeric(rowsum(weights * y, group_factor, reorder = FALSE)) / denominator
  weighted_x <- rowsum(x * weights, group_factor, reorder = FALSE) / denominator
  group_index <- as.integer(group_factor)
  list(
    x = x - weighted_x[group_index, , drop = FALSE],
    y = y - weighted_y[group_index]
  )
}


stage09a_fit_wls <- function(x, y, weights, group) {
  transformed <- stage09a_within_group_transform(x, y, weights, group)
  fit <- stats::lm.wfit(
    x = transformed$x,
    y = transformed$y,
    w = weights,
    method = "qr",
    singular.ok = FALSE
  )
  coefficients <- as.numeric(fit$coefficients)
  names(coefficients) <- colnames(x)
  if (length(coefficients) != ncol(x) || any(!is.finite(coefficients))) {
    stop("Stage 09A weighted fixed-effect fit produced non-finite coefficients", call. = FALSE)
  }
  list(
    coefficients = coefficients,
    rank = fit$rank,
    n_observations = nrow(x),
    n_groups = length(unique(group)),
    weighted_residual_sum_squares = sum(weights * fit$residuals^2)
  )
}


stage09a_predict_sequence_only <- function(fit, x) {
  stopifnot(is.matrix(x), identical(colnames(x), names(fit$coefficients)))
  as.numeric(x %*% fit$coefficients)
}
