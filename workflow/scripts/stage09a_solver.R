# Canonical Stage 09A adapter from project hyperparameters to R glmnet.

stage09a_glmnet_mapping <- function(spec_alpha, spec_l1_ratio, penalty_factor) {
  stopifnot(
    length(spec_alpha) == 1L,
    is.finite(spec_alpha),
    spec_alpha >= 0,
    length(spec_l1_ratio) == 1L,
    is.finite(spec_l1_ratio),
    spec_l1_ratio >= 0,
    spec_l1_ratio <= 1,
    length(penalty_factor) > 0L,
    all(penalty_factor %in% c(0, 1)),
    sum(penalty_factor) > 0
  )

  # glmnet rescales penalty.factor to sum to nvars. For raw {0,1}
  # factors, the internal multiplier on every penalized coefficient is
  # p_total / p_penalized. Multiplying lambda by the inverse factor keeps
  # the effective sequence penalty equal to the specification's alpha.
  p_total <- length(penalty_factor)
  p_penalized <- sum(penalty_factor)

  list(
    spec_alpha = spec_alpha,
    spec_l1_ratio = spec_l1_ratio,
    glmnet_lambda = spec_alpha * p_penalized / p_total,
    glmnet_alpha = spec_l1_ratio,
    p_total = p_total,
    p_penalized = p_penalized,
    penalty_factor_internal_scale = p_total / p_penalized
  )
}


stage09a_fit_accumulation <- function(
  x,
  y_centered,
  sample_weights,
  spec_alpha,
  spec_l1_ratio,
  convergence_threshold = 1e-14,
  max_iterations = 10000000L
) {
  stopifnot(
    is.matrix(x),
    nrow(x) == length(y_centered),
    nrow(x) == length(sample_weights),
    all(is.finite(x)),
    all(is.finite(y_centered)),
    all(is.finite(sample_weights)),
    all(sample_weights > 0)
  )
  penalty_factor <- rep(1, ncol(x))
  mapping <- stage09a_glmnet_mapping(spec_alpha, spec_l1_ratio, penalty_factor)

  fit <- glmnet::glmnet(
    x = x,
    y = y_centered,
    # The character family "gaussian" internally scales the response even
    # when lambda is supplied. A family object uses glmnet's supported GLM
    # pathway and optimizes the unscaled Gaussian loss required by v0.19.1.
    family = stats::gaussian(),
    weights = sample_weights,
    lambda = mapping$glmnet_lambda,
    alpha = mapping$glmnet_alpha,
    penalty.factor = penalty_factor,
    standardize = FALSE,
    intercept = FALSE,
    thresh = convergence_threshold,
    maxit = max_iterations
  )
  list(fit = fit, mapping = mapping)
}


stage09a_fit_representation <- function(
  design,
  represented,
  sample_weights,
  penalty_factor,
  spec_alpha,
  spec_l1_ratio,
  convergence_threshold = 1e-14,
  max_iterations = 10000000L
) {
  stopifnot(
    is.matrix(design),
    nrow(design) == length(represented),
    nrow(design) == length(sample_weights),
    ncol(design) == length(penalty_factor),
    all(represented %in% c(0, 1)),
    all(is.finite(sample_weights)),
    all(sample_weights > 0)
  )
  mapping <- stage09a_glmnet_mapping(spec_alpha, spec_l1_ratio, penalty_factor)

  fit <- glmnet::glmnet(
    x = design,
    y = represented,
    family = "binomial",
    weights = sample_weights,
    lambda = mapping$glmnet_lambda,
    alpha = mapping$glmnet_alpha,
    penalty.factor = penalty_factor,
    standardize = FALSE,
    intercept = FALSE,
    thresh = convergence_threshold,
    maxit = max_iterations
  )
  list(fit = fit, mapping = mapping)
}
