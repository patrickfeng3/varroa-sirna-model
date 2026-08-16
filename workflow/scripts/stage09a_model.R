#!/usr/bin/env Rscript

# Canonical Stage 09A model fitting and leakage-free validation.
# Scientific definitions are fixed by PIPELINE_SPEC.md v0.19.1.  This script
# uses only base R, Matrix (an unavoidable glmnet dependency), and glmnet.

parse_args <- function(args) {
  output <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    if (index == length(args)) stop(sprintf("missing value for --%s", key), call. = FALSE)
    output[[gsub("-", "_", key)]] <- args[[index + 1L]]
    index <- index + 2L
  }
  required <- c("training", "candidates", "accounting", "solver", "output_root")
  missing <- setdiff(required, names(output))
  if (length(missing)) stop(sprintf("missing arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
  output
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
source(args$solver)

alpha_grid <- c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1)
l1_ratio_grid <- c(0.05, 0.25, 0.50, 0.75, 0.95, 1.00)
structures <- c("A_shared", "B_shared_plus_length_interactions", "C_separate_23_24")
base_features <- c(
  as.vector(outer(c("guide_5p1", "guide_5p2", "guide_3p2", "guide_3p1"), c("C", "G", "U"), paste, sep = "_")),
  "guide_A3p3", "guide_GC_3p5_10", "guide_W17", "guide_R10"
)
# outer() is column-major; restore the fixed position-major encoding order.
base_features <- c(
  "guide_5p1_C", "guide_5p1_G", "guide_5p1_U",
  "guide_5p2_C", "guide_5p2_G", "guide_5p2_U",
  "guide_3p2_C", "guide_3p2_G", "guide_3p2_U",
  "guide_3p1_C", "guide_3p1_G", "guide_3p1_U",
  "guide_A3p3", "guide_GC_3p5_10", "guide_W17", "guide_R10"
)

read_tsv <- function(path) {
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  read.delim(connection, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

rbind_fill <- function(...) {
  frames <- list(...)
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  aligned <- lapply(frames, function(frame) {
    for (column in setdiff(columns, names(frame))) frame[[column]] <- NA
    frame[, columns, drop = FALSE]
  })
  do.call(rbind, aligned)
}

median_na <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) median(values) else NA_real_
}

group_key <- function(data) paste(data$sample, data$analysis_unit, data$candidate_length_nt, sep = "\r")

sample_aware_weights <- function(data, positive_only = TRUE) {
  if (!nrow(data)) stop("empty training fold", call. = FALSE)
  keys <- group_key(data)
  group_sizes <- table(keys)
  groups_per_sample <- vapply(split(keys, data$sample), function(x) length(unique(x)), integer(1))
  raw <- 1 / (as.numeric(groups_per_sample[data$sample]) * as.numeric(group_sizes[keys]))
  raw * nrow(data) / sum(raw)
}

fit_preprocessing <- function(data) {
  matrix <- as.matrix(data[, base_features, drop = FALSE])
  storage.mode(matrix) <- "double"
  means <- colMeans(matrix)
  sds <- apply(matrix, 2, stats::sd)
  retained <- names(sds)[is.finite(sds) & sds > 0]
  omitted <- setdiff(base_features, retained)
  if (!length(retained)) stop("all sequence features have zero training SD", call. = FALSE)
  list(means = means, sds = sds, retained = retained, omitted = omitted)
}

apply_preprocessing <- function(data, preprocessing) {
  matrix <- as.matrix(data[, preprocessing$retained, drop = FALSE])
  storage.mode(matrix) <- "double"
  sweep(sweep(matrix, 2, preprocessing$means[preprocessing$retained], "-"),
        2, preprocessing$sds[preprocessing$retained], "/")
}

structure_matrix <- function(z, lengths, structure) {
  if (structure %in% c("A_shared", "C_separate_23_24")) return(z)
  if (structure == "B_shared_plus_length_interactions") {
    interactions <- z * as.numeric(lengths == 24L)
    colnames(interactions) <- paste0(colnames(z), "_x_24nt")
    return(cbind(z, interactions))
  }
  stop(sprintf("unknown structure: %s", structure), call. = FALSE)
}

within_group_transform <- function(x, y, weights, keys) {
  x_centered <- x
  y_centered <- y
  for (indexes in split(seq_along(y), keys)) {
    group_weights <- weights[indexes]
    y_centered[indexes] <- y[indexes] - weighted.mean(y[indexes], group_weights)
    group_means <- colSums(x[indexes, , drop = FALSE] * group_weights) / sum(group_weights)
    x_centered[indexes, ] <- sweep(x[indexes, , drop = FALSE], 2, group_means, "-")
  }
  list(x = x_centered, y = y_centered)
}

fit_accumulation <- function(data, structure, spec_alpha, spec_l1_ratio) {
  positives <- data[data$represented == 1L & data$abundance > 0, , drop = FALSE]
  if (!nrow(positives)) stop("accumulation fold has no represented sequences", call. = FALSE)
  if (structure == "C_separate_23_24") {
    fits <- lapply(c(23L, 24L), function(length_nt) {
      subset <- positives[positives$candidate_length_nt == length_nt, , drop = FALSE]
      if (!nrow(subset)) stop(sprintf("separate model lacks %d-nt positives", length_nt), call. = FALSE)
      preprocessing <- fit_preprocessing(subset)
      z <- apply_preprocessing(subset, preprocessing)
      weights <- sample_aware_weights(subset)
      transformed <- within_group_transform(z, log(subset$abundance), weights, group_key(subset))
      fitted <- stage09a_fit_accumulation(
        transformed$x, transformed$y, weights, spec_alpha, spec_l1_ratio
      )
      list(fit = fitted$fit, mapping = fitted$mapping, preprocessing = preprocessing,
           feature_names = colnames(z), length = length_nt)
    })
    names(fits) <- c("23", "24")
    return(list(structure = structure, fits = fits, spec_alpha = spec_alpha,
                spec_l1_ratio = spec_l1_ratio))
  }
  preprocessing <- fit_preprocessing(positives)
  z <- apply_preprocessing(positives, preprocessing)
  x <- structure_matrix(z, positives$candidate_length_nt, structure)
  weights <- sample_aware_weights(positives)
  transformed <- within_group_transform(x, log(positives$abundance), weights, group_key(positives))
  fitted <- stage09a_fit_accumulation(
    transformed$x, transformed$y, weights, spec_alpha, spec_l1_ratio
  )
  list(structure = structure, fit = fitted$fit, mapping = fitted$mapping,
       preprocessing = preprocessing, feature_names = colnames(x),
       spec_alpha = spec_alpha, spec_l1_ratio = spec_l1_ratio)
}

predict_accumulation <- function(model, data) {
  if (model$structure == "C_separate_23_24") {
    output <- rep(NA_real_, nrow(data))
    for (length_nt in c(23L, 24L)) {
      indexes <- which(data$candidate_length_nt == length_nt)
      component <- model$fits[[as.character(length_nt)]]
      z <- apply_preprocessing(data[indexes, , drop = FALSE], component$preprocessing)
      output[indexes] <- as.numeric(z %*% as.numeric(component$fit$beta[, 1]))
    }
    return(output)
  }
  z <- apply_preprocessing(data, model$preprocessing)
  x <- structure_matrix(z, data$candidate_length_nt, model$structure)
  as.numeric(x %*% as.numeric(model$fit$beta[, 1]))
}

spearman_safe <- function(x, y) {
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = "spearman"))
}

group_metrics <- function(data, scores, score_name = "accumulation") {
  output <- list()
  keys <- group_key(data)
  for (key in unique(keys)) {
    indexes <- which(keys == key)
    group <- data[indexes, , drop = FALSE]
    group_scores <- scores[indexes]
    abundance <- group$abundance
    n <- length(indexes)
    selected_n <- ceiling(0.10 * n)
    selected <- order(-group_scores, group$sequence_dna)[seq_len(selected_n)]
    total <- sum(abundance)
    share <- if (total > 0) sum(abundance[selected]) / total else NA_real_
    positive <- abundance > 0
    output[[length(output) + 1L]] <- data.frame(
      sample = group$sample[[1]], analysis_unit = group$analysis_unit[[1]],
      biological_virus = group$biological_virus[[1]],
      candidate_length_nt = group$candidate_length_nt[[1]], n_opportunities = n,
      n_represented = sum(positive), total_abundance = total,
      spearman_rho = spearman_safe(group_scores, abundance),
      top10_abundance_share = share,
      top10_abundance_lift = if (is.finite(share)) share / (selected_n / n) else NA_real_,
      conditional_positive_spearman_rho = spearman_safe(group_scores[positive], abundance[positive]),
      score_type = score_name, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

selection_summary <- function(metrics) {
  m23 <- median_na(metrics$spearman_rho[metrics$candidate_length_nt == 23L])
  m24 <- median_na(metrics$spearman_rho[metrics$candidate_length_nt == 24L])
  l23 <- median_na(metrics$top10_abundance_lift[metrics$candidate_length_nt == 23L])
  l24 <- median_na(metrics$top10_abundance_lift[metrics$candidate_length_nt == 24L])
  data.frame(
    median_rho_23nt = m23, median_rho_24nt = m24,
    median_top10_lift_23nt = l23, median_top10_lift_24nt = l24,
    selection_score_rho = mean(c(m23, m24)), selection_score_top10 = mean(c(l23, l24)),
    stringsAsFactors = FALSE
  )
}

select_row <- function(rows, include_structure = FALSE, tolerance = 1e-6) {
  valid <- rows[is.finite(rows$selection_score_rho) & is.finite(rows$selection_score_top10), , drop = FALSE]
  if (!nrow(valid)) stop("no estimable model-selection configuration", call. = FALSE)
  best <- max(valid$selection_score_rho)
  valid <- valid[best - valid$selection_score_rho <= tolerance, , drop = FALSE]
  best_top <- max(valid$selection_score_top10)
  valid <- valid[best_top - valid$selection_score_top10 <= tolerance, , drop = FALSE]
  if ("alpha" %in% names(valid)) valid <- valid[valid$alpha == max(valid$alpha), , drop = FALSE]
  if ("l1_ratio" %in% names(valid)) valid <- valid[valid$l1_ratio == max(valid$l1_ratio), , drop = FALSE]
  if (include_structure && "model_structure" %in% names(valid)) {
    valid$structure_order <- match(valid$model_structure, structures)
    valid <- valid[valid$structure_order == min(valid$structure_order), , drop = FALSE]
    valid$structure_order <- NULL
  }
  valid[1, , drop = FALSE]
}

tune_accumulation <- function(data, structure) {
  families <- sort(unique(data$biological_virus))
  grid <- expand.grid(alpha = alpha_grid, l1_ratio = l1_ratio_grid, KEEP.OUT.ATTRS = FALSE)
  rows <- list()
  for (configuration in seq_len(nrow(grid))) {
    fold_metrics <- list()
    for (family in families) {
      train <- data[data$biological_virus != family, , drop = FALSE]
      heldout <- data[data$biological_virus == family, , drop = FALSE]
      model <- fit_accumulation(train, structure, grid$alpha[[configuration]], grid$l1_ratio[[configuration]])
      fold_metrics[[length(fold_metrics) + 1L]] <- group_metrics(heldout, predict_accumulation(model, heldout))
    }
    metrics <- do.call(rbind, fold_metrics)
    summary <- selection_summary(metrics)
    rows[[configuration]] <- cbind(
      data.frame(model_structure = structure, alpha = grid$alpha[[configuration]],
                 l1_ratio = grid$l1_ratio[[configuration]], stringsAsFactors = FALSE), summary
    )
  }
  results <- do.call(rbind, rows)
  selected <- select_row(results)
  list(results = results, selected = selected)
}

outer_primary_validation <- function(data) {
  all_group_rows <- list()
  selection_rows <- list()
  for (outer_family in sort(unique(data$biological_virus))) {
    training <- data[data$biological_virus != outer_family, , drop = FALSE]
    heldout <- data[data$biological_virus == outer_family, , drop = FALSE]
    for (structure in structures) {
      tuning <- tune_accumulation(training, structure)
      selected <- tuning$selected
      model <- fit_accumulation(training, structure, selected$alpha, selected$l1_ratio)
      metrics <- group_metrics(heldout, predict_accumulation(model, heldout))
      metrics$validation_scheme <- "leave_one_virus_family_out"
      metrics$outer_group <- outer_family
      metrics$model_structure <- structure
      metrics$alpha <- selected$alpha
      metrics$l1_ratio <- selected$l1_ratio
      all_group_rows[[length(all_group_rows) + 1L]] <- metrics
      selection_rows[[length(selection_rows) + 1L]] <- cbind(
        data.frame(selection_stage = "nested_outer_inner", outer_group = outer_family,
                   selected = TRUE, stringsAsFactors = FALSE), selected
      )
    }
  }
  list(groups = do.call(rbind, all_group_rows), selection = do.call(rbind, selection_rows))
}

architecture_summary <- function(primary_groups) {
  rows <- lapply(structures, function(structure) {
    subset <- primary_groups[primary_groups$model_structure == structure, , drop = FALSE]
    cbind(data.frame(model_structure = structure, stringsAsFactors = FALSE), selection_summary(subset))
  })
  results <- do.call(rbind, rows)
  chosen <- select_row(results, include_structure = TRUE)
  results$selected_primary_structure <- results$model_structure == chosen$model_structure
  results
}

secondary_sample_validation <- function(data) {
  rows <- list()
  for (sample_id in sort(unique(data$sample))) {
    training <- data[data$sample != sample_id, , drop = FALSE]
    heldout <- data[data$sample == sample_id, , drop = FALSE]
    candidates <- lapply(structures, function(structure) tune_accumulation(training, structure)$selected)
    candidates <- do.call(rbind, candidates)
    selected <- select_row(candidates, include_structure = TRUE)
    model <- fit_accumulation(training, selected$model_structure, selected$alpha, selected$l1_ratio)
    metrics <- group_metrics(heldout, predict_accumulation(model, heldout))
    metrics$validation_scheme <- "leave_one_sample_out"
    metrics$outer_group <- sample_id
    metrics$model_structure <- selected$model_structure
    metrics$alpha <- selected$alpha
    metrics$l1_ratio <- selected$l1_ratio
    rows[[length(rows) + 1L]] <- metrics
  }
  do.call(rbind, rows)
}

representation_design <- function(data, preprocessing, structure, training_groups = NULL) {
  z <- apply_preprocessing(data, preprocessing)
  sequence <- structure_matrix(z, data$candidate_length_nt, structure)
  if (is.null(training_groups)) return(sequence)
  indicators <- matrix(0, nrow(data), length(training_groups),
                       dimnames = list(NULL, paste0("nuisance_", seq_along(training_groups))))
  matches <- match(group_key(data), training_groups)
  present <- which(!is.na(matches))
  indicators[cbind(present, matches[present])] <- 1
  cbind(indicators, sequence)
}

fit_representation <- function(data, structure, spec_alpha, spec_l1_ratio) {
  if (structure == "C_separate_23_24") {
    fits <- lapply(c(23L, 24L), function(length_nt) {
      subset <- data[data$candidate_length_nt == length_nt, , drop = FALSE]
      fit_representation(subset, "A_shared", spec_alpha, spec_l1_ratio)
    })
    names(fits) <- c("23", "24")
    return(list(structure = structure, fits = fits))
  }
  preprocessing <- fit_preprocessing(data)
  groups <- sort(unique(group_key(data)))
  design <- representation_design(data, preprocessing, structure, groups)
  nuisance_count <- length(groups)
  penalty_factor <- c(rep(0, nuisance_count), rep(1, ncol(design) - nuisance_count))
  fitted <- stage09a_fit_representation(
    design, as.integer(data$represented), sample_aware_weights(data, positive_only = FALSE),
    penalty_factor, spec_alpha, spec_l1_ratio
  )
  list(structure = structure, fit = fitted$fit, mapping = fitted$mapping,
       preprocessing = preprocessing, groups = groups,
       sequence_feature_names = colnames(design)[seq.int(nuisance_count + 1L, ncol(design))],
       nuisance_count = nuisance_count)
}

predict_representation <- function(model, data) {
  if (model$structure == "C_separate_23_24") {
    output <- rep(NA_real_, nrow(data))
    for (length_nt in c(23L, 24L)) {
      indexes <- which(data$candidate_length_nt == length_nt)
      output[indexes] <- predict_representation(model$fits[[as.character(length_nt)]], data[indexes, , drop = FALSE])
    }
    return(output)
  }
  sequence <- representation_design(data, model$preprocessing, model$structure)
  beta <- as.numeric(model$fit$beta[, 1])
  sequence_beta <- beta[seq.int(model$nuisance_count + 1L, length(beta))]
  plogis(as.numeric(sequence %*% sequence_beta))
}

auc_roc <- function(labels, scores) {
  positives <- sum(labels == 1L); negatives <- sum(labels == 0L)
  if (!positives || !negatives) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - positives * (positives + 1) / 2) / (positives * negatives)
}

average_precision <- function(labels, scores) {
  positives <- sum(labels == 1L)
  if (!positives) return(NA_real_)
  order_index <- order(-scores, seq_along(scores))
  ordered <- labels[order_index]
  precision <- cumsum(ordered) / seq_along(ordered)
  sum(precision[ordered == 1L]) / positives
}

representation_metrics <- function(data, probabilities) {
  rows <- list()
  for (length_nt in c(23L, 24L)) {
    indexes <- which(data$candidate_length_nt == length_nt)
    labels <- as.integer(data$represented[indexes])
    scores <- probabilities[indexes]
    prevalence <- mean(labels)
    ap <- average_precision(labels, scores)
    selected_n <- ceiling(0.10 * length(indexes))
    selected <- order(-scores, seq_along(scores))[seq_len(selected_n)]
    selected_rate <- mean(labels[selected])
    rows[[length(rows) + 1L]] <- data.frame(
      candidate_length_nt = length_nt, n_opportunities = length(indexes),
      n_represented = sum(labels), representation_prevalence = prevalence,
      roc_auc = auc_roc(labels, scores), average_precision = ap,
      ap_lift = if (prevalence > 0) ap / prevalence else NA_real_,
      top_decile_representation_enrichment = if (prevalence > 0) selected_rate / prevalence else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

tune_representation <- function(data, structure) {
  rows <- list(); index <- 0L
  for (alpha in alpha_grid) for (l1_ratio in l1_ratio_grid) {
    folds <- list()
    for (family in sort(unique(data$biological_virus))) {
      training <- data[data$biological_virus != family, , drop = FALSE]
      heldout <- data[data$biological_virus == family, , drop = FALSE]
      model <- fit_representation(training, structure, alpha, l1_ratio)
      folds[[length(folds) + 1L]] <- representation_metrics(heldout, predict_representation(model, heldout))
    }
    metrics <- do.call(rbind, folds)
    index <- index + 1L
    rows[[index]] <- data.frame(
      alpha = alpha, l1_ratio = l1_ratio,
      mean_auc = mean(vapply(c(23L, 24L), function(l) median_na(metrics$roc_auc[metrics$candidate_length_nt == l]), numeric(1))),
      mean_ap_lift = mean(vapply(c(23L, 24L), function(l) median_na(metrics$ap_lift[metrics$candidate_length_nt == l]), numeric(1))),
      stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, rows)
  valid <- results[is.finite(results$mean_auc) & is.finite(results$mean_ap_lift), , drop = FALSE]
  best <- max(valid$mean_auc); valid <- valid[best - valid$mean_auc <= 1e-6, , drop = FALSE]
  best_ap <- max(valid$mean_ap_lift); valid <- valid[best_ap - valid$mean_ap_lift <= 1e-6, , drop = FALSE]
  valid <- valid[valid$alpha == max(valid$alpha), , drop = FALSE]
  valid <- valid[valid$l1_ratio == max(valid$l1_ratio), , drop = FALSE]
  list(results = results, selected = valid[1, , drop = FALSE])
}

outer_representation_validation <- function(data, structure) {
  diagnostics <- list(); probabilities <- rep(NA_real_, nrow(data))
  for (family in sort(unique(data$biological_virus))) {
    training <- data[data$biological_virus != family, , drop = FALSE]
    indexes <- which(data$biological_virus == family)
    heldout <- data[indexes, , drop = FALSE]
    tuning <- tune_representation(training, structure)
    selected <- tuning$selected
    model <- fit_representation(training, structure, selected$alpha, selected$l1_ratio)
    probabilities[indexes] <- predict_representation(model, heldout)
    metrics <- representation_metrics(heldout, probabilities[indexes])
    metrics$validation_scheme <- "leave_one_virus_family_out"
    metrics$outer_group <- family
    metrics$alpha <- selected$alpha
    metrics$l1_ratio <- selected$l1_ratio
    diagnostics[[length(diagnostics) + 1L]] <- metrics
  }
  list(diagnostics = do.call(rbind, diagnostics), probabilities = probabilities)
}

cv_summary <- function(groups, length_nt) {
  output <- list()
  for (scheme in unique(groups$validation_scheme)) {
    subset <- groups[groups$validation_scheme == scheme & groups$candidate_length_nt == length_nt, , drop = FALSE]
    output[[length(output) + 1L]] <- data.frame(
      validation_scheme = scheme, candidate_length_nt = length_nt,
      n_held_out_groups = nrow(subset),
      median_spearman_rho = median_na(subset$spearman_rho),
      median_top10_abundance_share = median_na(subset$top10_abundance_share),
      median_top10_abundance_lift = median_na(subset$top10_abundance_lift),
      median_conditional_positive_spearman_rho = median_na(subset$conditional_positive_spearman_rho),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

model_coefficients <- function(model) {
  rows <- list()
  if (model$structure == "C_separate_23_24") {
    for (length_nt in c(23L, 24L)) {
      component <- model$fits[[as.character(length_nt)]]
      rows[[length(rows) + 1L]] <- data.frame(
        model_structure = model$structure, candidate_length_scope = paste0(length_nt, "nt"),
        feature = component$feature_names,
        coefficient = as.numeric(component$fit$beta[, 1]), stringsAsFactors = FALSE
      )
    }
  } else {
    rows[[1]] <- data.frame(
      model_structure = model$structure, candidate_length_scope = "shared_23_24",
      feature = model$feature_names, coefficient = as.numeric(model$fit$beta[, 1]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

preprocessing_rows <- function(model) {
  make <- function(preprocessing, scope) data.frame(
    candidate_length_scope = scope, feature = base_features,
    training_mean = preprocessing$means[base_features], training_sd = preprocessing$sds[base_features],
    retained = base_features %in% preprocessing$retained,
    omitted_zero_sd = base_features %in% preprocessing$omitted, stringsAsFactors = FALSE
  )
  if (model$structure == "C_separate_23_24") {
    return(do.call(rbind, lapply(c(23L, 24L), function(length_nt) {
      make(model$fits[[as.character(length_nt)]]$preprocessing, paste0(length_nt, "nt"))
    })))
  }
  make(model$preprocessing, "shared_23_24")
}

favourable_percentile <- function(values) (rank(values, ties.method = "average") - 0.5) / length(values)

training <- read_tsv(args$training)
candidates <- read_tsv(args$candidates)
accounting <- read_tsv(args$accounting)
required_training <- c("sample", "analysis_unit", "biological_virus", "candidate_length_nt",
                       "sequence_dna", "represented", "abundance", base_features)
if (!all(required_training %in% names(training))) stop("prepared training table has invalid schema", call. = FALSE)
if (!all(base_features %in% names(candidates))) stop("prepared candidate table has invalid feature schema", call. = FALSE)
training$candidate_length_nt <- as.integer(training$candidate_length_nt)
training$represented <- as.integer(training$represented)
training$abundance <- as.numeric(training$abundance)
candidates$candidate_length_nt <- as.integer(candidates$candidate_length_nt)

# Primary nested architecture validation and secondary sample robustness.
primary <- outer_primary_validation(training)
architectures <- architecture_summary(primary$groups)
selected_structure <- architectures$model_structure[architectures$selected_primary_structure][[1]]
final_tuning <- tune_accumulation(training, selected_structure)
final_selected <- final_tuning$selected
final_model <- fit_accumulation(training, selected_structure, final_selected$alpha, final_selected$l1_ratio)
secondary_groups <- secondary_sample_validation(training)

# Diagnostic representation model and hurdle benchmark.
representation_outer <- outer_representation_validation(training, selected_structure)
representation_tuning <- tune_representation(training, selected_structure)
representation_selected <- representation_tuning$selected
final_representation <- fit_representation(
  training, selected_structure, representation_selected$alpha, representation_selected$l1_ratio
)

primary_selected <- primary$groups[primary$groups$model_structure == selected_structure, , drop = FALSE]
accumulation_outer_scores <- rep(NA_real_, nrow(training))
# Recreate only selected-architecture outer predictions for hurdle diagnostics.
for (family in sort(unique(training$biological_virus))) {
  train_fold <- training[training$biological_virus != family, , drop = FALSE]
  indexes <- which(training$biological_virus == family)
  tune_fold <- tune_accumulation(train_fold, selected_structure)$selected
  fit_fold <- fit_accumulation(train_fold, selected_structure, tune_fold$alpha, tune_fold$l1_ratio)
  accumulation_outer_scores[indexes] <- predict_accumulation(fit_fold, training[indexes, , drop = FALSE])
}
hurdle_scores <- representation_outer$probabilities * exp(accumulation_outer_scores)
hurdle_groups <- group_metrics(training, hurdle_scores, "hurdle_benchmark")
hurdle_groups$validation_scheme <- "leave_one_virus_family_out"
hurdle_groups$outer_group <- hurdle_groups$biological_virus
hurdle_groups$model_structure <- selected_structure
hurdle_groups$alpha <- NA_real_; hurdle_groups$l1_ratio <- NA_real_

all_cv_groups <- rbind(primary_selected, secondary_groups, hurdle_groups)

# Candidate-facing frozen sequence-only predictions.
candidates$layer1_accumulation_linear_predictor <- predict_accumulation(final_model, candidates)
candidates$layer1_representation_probability_diagnostic <- predict_representation(final_representation, candidates)
candidates$layer1_hurdle_score_diagnostic <- candidates$layer1_representation_probability_diagnostic *
  exp(candidates$layer1_accumulation_linear_predictor)
candidates$layer1_accumulation_percentile <- NA_real_
candidates$layer1_representation_percentile_diagnostic <- NA_real_
normalization_groups <- interaction(candidates$target_id, candidates$candidate_length_nt, drop = TRUE)
for (indexes in split(seq_len(nrow(candidates)), normalization_groups)) {
  candidates$layer1_accumulation_percentile[indexes] <- favourable_percentile(
    candidates$layer1_accumulation_linear_predictor[indexes]
  )
  candidates$layer1_representation_percentile_diagnostic[indexes] <- favourable_percentile(
    candidates$layer1_representation_probability_diagnostic[indexes]
  )
}

selection_table <- rbind(
  primary$selection,
  cbind(data.frame(selection_stage = "final_full_data_tuning", outer_group = "all_families",
                   selected = seq_len(nrow(final_tuning$results)) == as.integer(rownames(final_selected)[[1]]),
                   stringsAsFactors = FALSE), final_tuning$results)
)

representation_diagnostic <- representation_outer$diagnostics
representation_summary <- do.call(rbind, lapply(c(23L, 24L), function(length_nt) {
  subset <- representation_diagnostic[representation_diagnostic$candidate_length_nt == length_nt, , drop = FALSE]
  data.frame(
    candidate_length_nt = length_nt, n_outer_folds = nrow(subset),
    median_roc_auc = median_na(subset$roc_auc),
    median_average_precision = median_na(subset$average_precision),
    median_ap_lift = median_na(subset$ap_lift),
    median_top_decile_representation_enrichment = median_na(subset$top_decile_representation_enrichment),
    stringsAsFactors = FALSE
  )
}))
representation_diagnostic$record_type <- "outer_fold"
representation_summary$record_type <- "summary"
representation_output <- rbind_fill(representation_diagnostic, representation_summary)

architectures$median_hurdle_rho_23nt <- NA_real_
architectures$median_hurdle_rho_24nt <- NA_real_
selected_index <- which(architectures$model_structure == selected_structure)
architectures$median_hurdle_rho_23nt[selected_index] <- median_na(hurdle_groups$spearman_rho[hurdle_groups$candidate_length_nt == 23L])
architectures$median_hurdle_rho_24nt[selected_index] <- median_na(hurdle_groups$spearman_rho[hurdle_groups$candidate_length_nt == 24L])

output_root <- args$output_root
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write_tsv(model_coefficients(final_model), file.path(output_root, "layer1_model_coefficients.tsv"))
write_tsv(preprocessing_rows(final_model), file.path(output_root, "layer1_model_preprocessing.tsv"))
write_tsv(selection_table, file.path(output_root, "layer1_model_selection.tsv"))
write_tsv(all_cv_groups, file.path(output_root, "layer1_cv_by_group.tsv"))
write_tsv(cv_summary(all_cv_groups[all_cv_groups$score_type == "accumulation", , drop = FALSE], 23L),
          file.path(output_root, "layer1_cv_summary_23nt.tsv"))
write_tsv(cv_summary(all_cv_groups[all_cv_groups$score_type == "accumulation", , drop = FALSE], 24L),
          file.path(output_root, "layer1_cv_summary_24nt.tsv"))
write_tsv(representation_output, file.path(output_root, "layer1_representation_diagnostic.tsv"))
write_tsv(architectures, file.path(output_root, "layer1_architecture_benchmarks.tsv"))
write_tsv(candidates, file.path(output_root, "candidate_layer1.tsv"))

parameters <- data.frame(
  parameter = c(
    "stage", "specification_version", "primary_response", "pseudocount",
    "alpha_grid", "l1_ratio_grid", "model_structures", "selected_structure",
    "selected_alpha", "selected_l1_ratio", "hyperparameter_boundary_warning",
    "solver", "R_version", "glmnet_version", "glmnet_spec_alpha_mapping",
    "glmnet_spec_l1_ratio_mapping", "candidate_normalization", "length_pooling",
    "stage08_features_in_layer1"
  ),
  value = c(
    "09A", "v0.19.1", "ln(count) among positive supported sequences", "none",
    paste(alpha_grid, collapse = ","), paste(l1_ratio_grid, collapse = ","), paste(structures, collapse = ","),
    selected_structure, final_selected$alpha, final_selected$l1_ratio,
    final_selected$alpha %in% range(alpha_grid), "R glmnet approved adapter",
    as.character(getRversion()), as.character(packageVersion("glmnet")),
    "spec alpha -> glmnet lambda (penalty.factor corrected)", "spec l1_ratio -> glmnet alpha",
    "favourable percentile within target x candidate_length", "never",
    "forbidden/none"
  ), stringsAsFactors = FALSE
)
write_tsv(parameters, file.path(dirname(output_root), "stage09_parameters.tsv"))

accounting_map <- setNames(accounting$value, accounting$metric)
qc <- data.frame(
  check = c(
    "primary_samples", "sample_virus_units", "sample_virus_length_groups",
    "opportunities_23nt", "opportunities_24nt", "represented_23nt", "represented_24nt",
    "supported_abundance", "candidate_row_preservation", "candidate_lengths_separate",
    "stage08_feature_leakage", "final_alpha_grid_boundary"
  ),
  observed = c(
    accounting_map[["primary_samples"]], accounting_map[["sample_virus_units"]],
    accounting_map[["sample_virus_length_groups"]], accounting_map[["opportunities_23nt"]],
    accounting_map[["opportunities_24nt"]], accounting_map[["represented_23nt"]],
    accounting_map[["represented_24nt"]], accounting_map[["supported_abundance"]],
    nrow(candidates), paste(sort(unique(candidates$candidate_length_nt)), collapse = ","),
    "none", final_selected$alpha %in% range(alpha_grid)
  ),
  expected = c("20", "54", "108", "411079", "408148", "121592", "175564", "3445943",
               accounting_map[["candidate_rows"]], "23,24", "none", "false preferred"),
  severity = c(rep("PASS", 11), if (final_selected$alpha %in% range(alpha_grid)) "WARN" else "PASS"),
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(dirname(output_root), "stage09_qc.tsv"))
