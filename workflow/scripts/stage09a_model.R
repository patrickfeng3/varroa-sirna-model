#!/usr/bin/env Rscript

# Canonical Stage 09A v0.20: two transparent weighted fixed-effect models.

parse_args <- function(args) {
  output <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[[index]]))
    if (index == length(args)) stop(sprintf("missing value for --%s", key), call. = FALSE)
    output[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  required <- c("training", "candidates", "accounting", "wls_helper", "output_root")
  missing <- setdiff(required, names(output))
  if (length(missing)) stop(sprintf("missing arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
  output
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
source(args$wls_helper)

feature_names <- c(
  "guide_5p1_C", "guide_5p1_G", "guide_5p1_U",
  "guide_5p2_C", "guide_5p2_G", "guide_5p2_U",
  "guide_3p2_C", "guide_3p2_G", "guide_3p2_U",
  "guide_3p1_C", "guide_3p1_G", "guide_3p1_U",
  "guide_A3p3", "guide_GC_3p5_10", "guide_W17", "guide_R10"
)
lengths <- c(23L, 24L)

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

median_finite <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) stats::median(values) else NA_real_
}

range_finite <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) range(values) else c(NA_real_, NA_real_)
}

sequence_matrix <- function(data) {
  matrix <- as.matrix(data[, feature_names, drop = FALSE])
  storage.mode(matrix) <- "double"
  if (any(!is.finite(matrix))) stop("non-finite Stage 09A predictor", call. = FALSE)
  matrix
}

positive_data <- function(data, length_nt) {
  data[data$candidate_length_nt == length_nt & data$represented == 1L & data$abundance > 0, , drop = FALSE]
}

fit_length_model <- function(data, length_nt) {
  positives <- positive_data(data, length_nt)
  if (!nrow(positives)) stop(sprintf("no positive %d-nt training opportunities", length_nt), call. = FALSE)
  group <- paste(positives$sample, positives$analysis_unit, sep = "\r")
  weights <- stage09a_sample_weights(positives$sample, group)
  fit <- stage09a_fit_wls(sequence_matrix(positives), log(positives$abundance), weights, group)
  fit$length_nt <- length_nt
  fit
}

spearman_safe <- function(x, y) {
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = "spearman"))
}

evaluate_groups <- function(data, scores, heldout_family, length_nt) {
  output <- list()
  keys <- paste(data$sample, data$analysis_unit, sep = "\r")
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
    output[[length(output) + 1L]] <- data.frame(
      heldout_biological_virus = heldout_family,
      sample = group$sample[[1]], analysis_unit = group$analysis_unit[[1]],
      candidate_length_nt = length_nt, n_opportunities = n,
      n_represented = sum(group$represented == 1L), observed_abundance = total,
      spearman_rho = spearman_safe(group_scores, abundance),
      top10_abundance_share = share,
      top10_abundance_lift = if (is.finite(share)) share / (selected_n / n) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

favourable_percentile <- function(values) {
  (rank(values, ties.method = "average") - 0.5) / length(values)
}

training <- read_tsv(args$training)
candidates <- read_tsv(args$candidates)
accounting <- read_tsv(args$accounting)
required <- c("sample", "analysis_unit", "biological_virus", "candidate_length_nt",
              "sequence_dna", "represented", "abundance", feature_names)
if (!all(required %in% names(training))) stop("invalid prepared Stage 09A training schema", call. = FALSE)
if (!all(feature_names %in% names(candidates))) stop("invalid prepared Stage 09A candidate schema", call. = FALSE)
training$candidate_length_nt <- as.integer(training$candidate_length_nt)
training$represented <- as.integer(training$represented)
training$abundance <- as.numeric(training$abundance)
candidates$candidate_length_nt <- as.integer(candidates$candidate_length_nt)
families <- sort(unique(training$biological_virus))
if (length(families) != 5L) stop(sprintf("expected five virus/family holdouts, observed %d", length(families)), call. = FALSE)

# Exactly five leave-one-virus/family-out fits per length.
holdout_fits <- list()
group_metrics <- list()
holdout_summaries <- list()
for (length_nt in lengths) {
  for (family in families) {
    training_fold <- training[training$biological_virus != family, , drop = FALSE]
    heldout <- training[training$biological_virus == family & training$candidate_length_nt == length_nt, , drop = FALSE]
    fit <- fit_length_model(training_fold, length_nt)
    holdout_fits[[paste(length_nt, family, sep = "\r")]] <- fit
    metrics <- evaluate_groups(heldout, stage09a_predict_sequence_only(fit, sequence_matrix(heldout)), family, length_nt)
    group_metrics[[length(group_metrics) + 1L]] <- metrics
    holdout_summaries[[length(holdout_summaries) + 1L]] <- data.frame(
      heldout_biological_virus = family, candidate_length_nt = length_nt,
      record_type = "virus_holdout_summary", sample = NA_character_, analysis_unit = NA_character_,
      n_opportunities = sum(metrics$n_opportunities), n_represented = sum(metrics$n_represented),
      observed_abundance = sum(metrics$observed_abundance),
      spearman_rho = median_finite(metrics$spearman_rho), top10_abundance_share = NA_real_,
      top10_abundance_lift = median_finite(metrics$top10_abundance_lift),
      stringsAsFactors = FALSE
    )
  }
}
group_output <- do.call(rbind, group_metrics)
group_output$record_type <- "sample_virus_group"
group_output <- group_output[, names(holdout_summaries[[1]]), drop = FALSE]
holdout_output <- rbind(group_output, do.call(rbind, holdout_summaries))

# Exactly one final full-data model per length.
final_models <- setNames(lapply(lengths, function(length_nt) fit_length_model(training, length_nt)), as.character(lengths))

coefficient_tables <- list()
stability_rows <- list()
for (length_nt in lengths) {
  final <- final_models[[as.character(length_nt)]]$coefficients
  holdout_matrix <- do.call(cbind, lapply(families, function(family) {
    holdout_fits[[paste(length_nt, family, sep = "\r")]]$coefficients
  }))
  colnames(holdout_matrix) <- families
  coefficient_tables[[as.character(length_nt)]] <- data.frame(
    candidate_length_nt = length_nt, coefficient = names(final), estimate = as.numeric(final),
    model_family = "weighted_fixed_effect_linear_regression",
    reference_category = ifelse(grepl("guide_(5p1|5p2|3p2|3p1)_", names(final)), "A", NA),
    stringsAsFactors = FALSE
  )
  stability_rows[[length(stability_rows) + 1L]] <- data.frame(
    candidate_length_nt = length_nt, coefficient = names(final),
    final_full_data_coefficient = as.numeric(final),
    median_holdout_fit_coefficient = apply(holdout_matrix, 1, median_finite),
    minimum_holdout_fit_coefficient = apply(holdout_matrix, 1, min),
    maximum_holdout_fit_coefficient = apply(holdout_matrix, 1, max),
    n_holdout_fits_same_sign_as_final = vapply(seq_along(final), function(index) {
      sum(sign(holdout_matrix[index, ]) == sign(final[[index]]))
    }, integer(1)), n_holdout_fits = ncol(holdout_matrix), stringsAsFactors = FALSE
  )
}

summary_for_length <- function(length_nt) {
  rows <- do.call(rbind, holdout_summaries)
  rows <- rows[rows$candidate_length_nt == length_nt, , drop = FALSE]
  rho_range <- range_finite(rows$spearman_rho)
  lift_range <- range_finite(rows$top10_abundance_lift)
  data.frame(
    candidate_length_nt = length_nt, n_virus_holdouts = nrow(rows),
    crossvirus_median_rho = median_finite(rows$spearman_rho),
    minimum_virus_holdout_rho = rho_range[[1]], maximum_virus_holdout_rho = rho_range[[2]],
    positive_rho_holdouts = sum(rows$spearman_rho > 0, na.rm = TRUE),
    valid_rho_holdouts = sum(is.finite(rows$spearman_rho)),
    crossvirus_median_top10_lift = median_finite(rows$top10_abundance_lift),
    minimum_virus_holdout_top10_lift = lift_range[[1]],
    maximum_virus_holdout_top10_lift = lift_range[[2]],
    top10_lift_gt1_holdouts = sum(rows$top10_abundance_lift > 1, na.rm = TRUE),
    valid_top10_lift_holdouts = sum(is.finite(rows$top10_abundance_lift)),
    stringsAsFactors = FALSE
  )
}

candidates$layer1_accumulation_linear_predictor <- NA_real_
for (length_nt in lengths) {
  indexes <- which(candidates$candidate_length_nt == length_nt)
  candidates$layer1_accumulation_linear_predictor[indexes] <- stage09a_predict_sequence_only(
    final_models[[as.character(length_nt)]], sequence_matrix(candidates[indexes, , drop = FALSE])
  )
}
candidates$layer1_accumulation_percentile <- NA_real_
normalization_group <- interaction(candidates$target_id, candidates$candidate_length_nt, drop = TRUE)
for (indexes in split(seq_len(nrow(candidates)), normalization_group)) {
  candidates$layer1_accumulation_percentile[indexes] <- favourable_percentile(
    candidates$layer1_accumulation_linear_predictor[indexes]
  )
}

output_root <- args$output_root
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write_tsv(accounting, file.path(output_root, "layer1_training_accounting.tsv"))
write_tsv(coefficient_tables[["23"]], file.path(output_root, "layer1_coefficients_23nt.tsv"))
write_tsv(coefficient_tables[["24"]], file.path(output_root, "layer1_coefficients_24nt.tsv"))
write_tsv(holdout_output, file.path(output_root, "layer1_leave_one_virus_out.tsv"))
write_tsv(summary_for_length(23L), file.path(output_root, "layer1_cv_summary_23nt.tsv"))
write_tsv(summary_for_length(24L), file.path(output_root, "layer1_cv_summary_24nt.tsv"))
write_tsv(do.call(rbind, stability_rows), file.path(output_root, "layer1_coefficient_stability.tsv"))
write_tsv(candidates, file.path(output_root, "candidate_layer1.tsv"))

provenance <- data.frame(
  parameter = c(
    "stage", "specification_version", "model_family", "response", "penalty",
    "hyperparameter_tuning", "model_structure_search", "length_models",
    "virus_family_holdouts", "primary_accumulation_fit_count", "R_version",
    "candidate_predictor", "nuisance_intercept_transfer", "representation_model",
    "hurdle_model", "stage08_features"
  ),
  value = c(
    "09A", "v0.20", "weighted_fixed_effect_linear_regression", "ln(count)", "none",
    "none", "none", "separate_23nt_and_24nt", length(families), 12L,
    as.character(getRversion()), "sequence_coefficients_only", "none", "not_fitted",
    "not_fitted", "excluded"
  ), stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(output_root, "layer1_model_provenance.tsv"))
write_tsv(provenance, file.path(dirname(output_root), "stage09_parameters.tsv"))

accounting_map <- setNames(accounting$value, accounting$metric)
qc <- data.frame(
  check = c(
    "primary_samples", "sample_virus_units", "opportunities_23nt", "opportunities_24nt",
    "represented_23nt", "represented_24nt", "supported_abundance",
    "outside_background_species", "outside_background_abundance", "virus_family_holdouts",
    "primary_accumulation_fit_count", "candidate_rows", "candidate_lengths_separate",
    "candidate_scores_finite", "candidate_percentiles_finite", "stage08_feature_leakage",
    "representation_model_absent", "hurdle_model_absent"
  ),
  observed = c(
    accounting_map[["primary_samples"]], accounting_map[["sample_virus_units"]],
    accounting_map[["opportunities_23nt"]], accounting_map[["opportunities_24nt"]],
    accounting_map[["represented_23nt"]], accounting_map[["represented_24nt"]],
    accounting_map[["supported_abundance"]], accounting_map[["outside_background_species"]],
    accounting_map[["outside_background_abundance"]], length(families), 12L, nrow(candidates),
    paste(sort(unique(candidates$candidate_length_nt)), collapse = ","),
    all(is.finite(candidates$layer1_accumulation_linear_predictor)),
    all(is.finite(candidates$layer1_accumulation_percentile)), "none", "absent", "absent"
  ),
  expected = c(
    "20", "54", "411079", "408148", "121592", "175564", "3445943", "3616", "3973",
    "5", "12", accounting_map[["candidate_rows"]], "23,24", "TRUE", "TRUE", "none", "absent", "absent"
  ),
  severity = "PASS", stringsAsFactors = FALSE
)
write_tsv(qc, file.path(dirname(output_root), "stage09_qc.tsv"))
