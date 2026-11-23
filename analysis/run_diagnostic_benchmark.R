#!/usr/bin/env Rscript

# Focused diagnostic simulations for
# "Nonparametric Change-Point Detection for Modern Data".
#
# The study is intentionally organized by the representation--search--
# calibration pipeline.  It is not a leaderboard of named software packages.
# All thresholds used for the representation study are learned from an
# independent null sample; evaluation nulls and alternatives are disjoint.

suppressPackageStartupMessages({
  library(data.table)
  library(Rcpp)
  library(ggplot2)
  library(gridExtra)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(root_dir)
plot_only <- identical(Sys.getenv("CPD_PLOT_ONLY", "0"), "1")
if (!plot_only) {
  sourceCpp(file.path("analysis", "benchmark_core.cpp"))
}

mode <- tolower(Sys.getenv("CPD_MODE", "pilot"))
if (!mode %in% c("pilot", "full")) stop("CPD_MODE must be pilot or full")
cores_requested <- as.integer(Sys.getenv("CPD_CORES", "4"))
if (is.na(cores_requested)) cores_requested <- 1L
cores <- max(1L, cores_requested)

cfg <- if (mode == "pilot") {
  list(
    amoc_cal = 35L, amoc_eval = 30L, amoc_alt = 30L,
    object_cal = 30L, object_eval = 25L, object_alt = 25L,
    dependence_null = 30L, dependence_alt = 25L, dependence_B = 49L,
    runtime_reps = 5L
  )
} else {
  list(
    amoc_cal = 10000L, amoc_eval = 2000L, amoc_alt = 500L,
    object_cal = 2000L, object_eval = 2000L, object_alt = 500L,
    dependence_null = 1000L, dependence_alt = 500L, dependence_B = 199L,
    runtime_reps = 30L
  )
}

dir.create(file.path("analysis", "results"), recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
suffix <- if (mode == "pilot") "_pilot" else ""
alpha <- 0.05
master_seed <- 20261118L
start_time <- Sys.time()

message(sprintf("Mode: %s; workers: %d", mode, cores))

wilson_interval <- function(successes, total, level = 0.95) {
  if (total == 0L) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level) / 2)
  phat <- successes / total
  denom <- 1 + z^2 / total
  center <- (phat + z^2 / (2 * total)) / denom
  half <- z * sqrt(phat * (1 - phat) / total + z^2 / (4 * total^2)) / denom
  c(max(0, center - half), min(1, center + half))
}

summarize_binary <- function(dt, group_cols, value_col = "reject") {
  dt[, {
    x <- sum(get(value_col))
    ci <- wilson_interval(x, .N)
    list(estimate = x / .N, ci_low = ci[1], ci_high = ci[2], n_rep = .N)
  }, by = group_cols]
}

pairwise_sqdist <- function(x) {
  squared_norm <- rowSums(x * x)
  d2 <- outer(squared_norm, squared_norm, "+") - 2 * tcrossprod(x)
  pmax(d2, 0)
}

cusum_paths <- function(x, cuts) {
  n <- nrow(x)
  p <- ncol(x)
  cumulative <- apply(x, 2L, cumsum)
  if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, ncol = 1L)
  left_sum <- cumulative[cuts, , drop = FALSE]
  total <- matrix(cumulative[n, ], nrow = length(cuts), ncol = p, byrow = TRUE)
  left_mean <- left_sum / cuts
  right_mean <- (total - left_sum) / (n - cuts)
  balance <- sqrt(cuts * (n - cuts) / n)
  sweep(left_mean - right_mean, 1L, balance, "*")
}

knn_graph <- function(d2, k = 5L) {
  n <- nrow(d2)
  k <- min(k, n - 1L)
  adjacency <- matrix(0, n, n)
  for (i in seq_len(n)) {
    candidates <- seq_len(n)[-i]
    neighbors <- candidates[order(d2[i, candidates])[seq_len(k)]]
    adjacency[i, neighbors] <- 1
  }
  1 * ((adjacency + t(adjacency)) > 0)
}

extract_scan <- function(result, cuts, fixed_cut) {
  fixed_index <- which.min(abs(cuts - fixed_cut))
  c(
    statistic = unname(result$statistic),
    location = unname(result$location),
    fixed = unname(result$values[fixed_index])
  )
}

compute_pipeline_statistics <- function(x, cuts = NULL, graph_k = 5L) {
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(cuts)) cuts <- seq.int(ceiling(0.2 * n), floor(0.8 * n))
  fixed_cut <- floor(n / 2)

  cpath <- cusum_paths(x, cuts)
  max_values <- apply(abs(cpath), 1L, max)
  dense_values <- abs(rowSums(cpath) / sqrt(p))
  max_index <- which.max(max_values)
  dense_index <- which.max(dense_values)

  d2 <- pairwise_sqdist(x)
  distance <- sqrt(d2)
  positive_distances <- d2[upper.tri(d2) & d2 > 0]
  bandwidth2 <- if (length(positive_distances)) median(positive_distances) else 1
  if (!is.finite(bandwidth2) || bandwidth2 <= 0) bandwidth2 <- 1
  kernel <- exp(-d2 / (2 * bandwidth2))
  graph <- knn_graph(d2, k = graph_k)

  energy <- extract_scan(scan_matrix_cpp(distance, cuts, 1L), cuts, fixed_cut)
  mmd <- extract_scan(scan_matrix_cpp(kernel, cuts, 2L), cuts, fixed_cut)
  graph_scan <- extract_scan(scan_matrix_cpp(graph, cuts, 3L), cuts, fixed_cut)
  fixed_index <- which.min(abs(cuts - fixed_cut))

  data.table(
    method = c("Max-CUSUM", "Dense projection", "Energy", "Gaussian MMD", "5-NN graph"),
    statistic = c(max_values[max_index], dense_values[dense_index],
                  energy["statistic"], mmd["statistic"], graph_scan["statistic"]),
    location = as.integer(c(cuts[max_index], cuts[dense_index],
                            energy["location"], mmd["location"], graph_scan["location"])),
    fixed = c(max_values[fixed_index], dense_values[fixed_index],
              energy["fixed"], mmd["fixed"], graph_scan["fixed"])
  )
}

generate_amoc <- function(n = 120L, p, scenario, tau = 60L) {
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  post <- seq.int(tau + 1L, n)

  if (scenario == "null") {
    return(x)
  } else if (scenario == "sparse") {
    x[post, 1L] <- x[post, 1L] + 0.80
  } else if (scenario == "dense") {
    x[post, ] <- x[post, ] + 0.65 / sqrt(p)
  } else if (scenario == "scale") {
    x[post, ] <- 1.25 * x[post, ]
  } else if (scenario == "correlation") {
    q <- min(5L, p)
    sigma <- diag(p)
    sigma[seq_len(q), seq_len(q)] <- 0.60
    diag(sigma) <- 1
    x[post, ] <- matrix(rnorm(length(post) * p), nrow = length(post)) %*% chol(sigma)
  } else if (scenario == "radial mixture") {
    low <- 0.65
    high <- sqrt((1 - 0.80 * low^2) / 0.20)
    scales <- ifelse(runif(length(post)) < 0.80, low, high)
    x[post, ] <- x[post, ] * scales
  } else {
    stop("unknown AMOC scenario: ", scenario)
  }
  x
}

parallel_replicates <- function(n_rep, seed_offset, worker) {
  indices <- seq_len(n_rep)
  pieces <- parallel::mclapply(
    indices,
    function(i) {
      set.seed(master_seed + seed_offset + i)
      ans <- worker(i)
      ans[, replicate := i]
      ans
    },
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  rbindlist(pieces, use.names = TRUE, fill = TRUE)
}

run_amoc_study <- function() {
  message("Running Euclidean AMOC representation study")
  dimensions <- c(5L, 50L, 200L)
  alternatives <- c("sparse", "dense", "scale", "correlation", "radial mixture")
  all_raw <- list()
  threshold_rows <- list()
  seed_cursor <- 0L

  for (p in dimensions) {
    message("  p = ", p, ": calibration null")
    cal <- parallel_replicates(cfg$amoc_cal, seed_cursor, function(i) {
      compute_pipeline_statistics(generate_amoc(p = p, scenario = "null"))
    })
    seed_cursor <- seed_cursor + 100000L
    cal[, `:=`(phase = "calibration", scenario = "null", p = p)]

    thresholds <- cal[, .(
      scan_threshold = as.numeric(quantile(statistic, 1 - alpha, type = 8)),
      fixed_threshold = as.numeric(quantile(fixed, 1 - alpha, type = 8))
    ), by = method]
    thresholds[, p := p]
    threshold_rows[[length(threshold_rows) + 1L]] <- thresholds
    all_raw[[length(all_raw) + 1L]] <- cal

    message("  p = ", p, ": evaluation null")
    evaluation <- parallel_replicates(cfg$amoc_eval, seed_cursor, function(i) {
      compute_pipeline_statistics(generate_amoc(p = p, scenario = "null"))
    })
    seed_cursor <- seed_cursor + 100000L
    evaluation[, `:=`(phase = "evaluation", scenario = "null", p = p)]
    all_raw[[length(all_raw) + 1L]] <- evaluation

    for (scenario_name in alternatives) {
      message("  p = ", p, ": ", scenario_name)
      alternative <- parallel_replicates(cfg$amoc_alt, seed_cursor, function(i) {
        x <- generate_amoc(p = p, scenario = scenario_name)
        statistics <- compute_pipeline_statistics(x)
        statistics[, scenario := scenario_name]
        if (p == 200L && scenario_name == "dense") {
          # A common reflection preserves all distances, the graph, and the
          # maximum absolute coordinate contrast on every replicate.
          x[, 101:200] <- -x[, 101:200]
          control <- compute_pipeline_statistics(x)
          control[, scenario := "orthogonal dense"]
          statistics <- rbind(statistics, control)
        }
        statistics
      })
      seed_cursor <- seed_cursor + 100000L
      alternative[, `:=`(phase = "alternative", p = p)]
      all_raw[[length(all_raw) + 1L]] <- alternative
    }
  }

  raw <- rbindlist(all_raw, use.names = TRUE)
  thresholds <- rbindlist(threshold_rows)
  setkey(thresholds, p, method)
  raw[thresholds, `:=`(
    scan_threshold = i.scan_threshold,
    fixed_threshold = i.fixed_threshold
  ), on = .(p, method)]
  raw[, reject := statistic > scan_threshold]
  raw[, abs_error := abs(location - 60L)]
  raw[, capped_error := fifelse(reject, abs_error, 120)]
  raw[scenario == "null", `:=`(abs_error = NA_integer_, capped_error = NA_real_)]

  summary <- raw[phase != "calibration", {
    successes <- sum(reject)
    ci <- wilson_interval(successes, .N)
    list(
      rejection = successes / .N,
      ci_low = ci[1],
      ci_high = ci[2],
      n_rep = .N,
      miss_rate = mean(!reject),
      mean_capped_error = mean(capped_error),
      median_error_detected = if (successes) as.numeric(median(abs_error[reject])) else NA_real_,
      mean_error_detected = if (successes) mean(abs_error[reject]) else NA_real_
    )
  }, by = .(scenario, p, method)]
  summary[scenario == "null", miss_rate := NA_real_]

  calibration <- raw[phase == "evaluation", {
    full_reject <- statistic > scan_threshold
    reused_reject <- statistic > fixed_threshold
    full_ci <- wilson_interval(sum(full_reject), .N)
    reused_ci <- wilson_interval(sum(reused_reject), .N)
    rbind(
      data.table(calibration = "Full-scan threshold", estimate = mean(full_reject),
                 ci_low = full_ci[1], ci_high = full_ci[2], n_rep = .N),
      data.table(calibration = "Fixed-split threshold reused", estimate = mean(reused_reject),
                 ci_low = reused_ci[1], ci_high = reused_ci[2], n_rep = .N)
    )
  }, by = .(p, method)]

  list(raw = raw, thresholds = thresholds, summary = summary, calibration = calibration)
}

scan_univariate <- function(x, cuts = seq.int(24L, 96L)) {
  n <- length(x)
  cumulative <- cumsum(x)
  contrast <- sqrt(cuts * (n - cuts) / n) *
    (cumulative[cuts] / cuts - (cumulative[n] - cumulative[cuts]) / (n - cuts))
  index <- which.max(abs(contrast))
  c(statistic = abs(contrast[index]), location = cuts[index])
}

generate_ar1 <- function(n = 120L, rho = 0.5, change = 0) {
  burn <- 200L
  innovations <- rnorm(n + burn, sd = sqrt(1 - rho^2))
  series <- numeric(n + burn)
  for (i in 2:length(series)) series[i] <- rho * series[i - 1L] + innovations[i]
  x <- tail(series, n)
  if (change != 0) x[61L:n] <- x[61L:n] + change
  x
}

block_permute <- function(x, block_length = 10L) {
  n <- length(x)
  if (n %% block_length != 0L) stop("n must be divisible by block_length")
  block_id <- rep(seq_len(n / block_length), each = block_length)
  blocks <- split(seq_len(n), block_id)
  x[unlist(blocks[sample.int(length(blocks))], use.names = FALSE)]
}

resampling_pvalues <- function(x, B, block_lengths = c(5L, 10L, 15L, 20L)) {
  observed <- scan_univariate(x)
  # Reordering the same observations can change summation roundoff. Count
  # numerical ties with the observed maximum as exceedances.
  lower_bound <- observed["statistic"] - 1e-12 * max(1, observed["statistic"])
  iid_stats <- numeric(B)
  block_stats <- matrix(NA_real_, nrow = B, ncol = length(block_lengths))
  for (b in seq_len(B)) {
    iid_stats[b] <- scan_univariate(sample(x, replace = FALSE))["statistic"]
    for (j in seq_along(block_lengths)) {
      block_stats[b, j] <- scan_univariate(
        block_permute(x, block_lengths[j])
      )["statistic"]
    }
  }
  answer <- c(
    statistic = unname(observed["statistic"]),
    location = unname(observed["location"]),
    p_iid = (1 + sum(iid_stats >= lower_bound)) / (B + 1)
  )
  for (j in seq_along(block_lengths)) {
    answer[paste0("p_block_", block_lengths[j])] <-
      (1 + sum(block_stats[, j] >= lower_bound)) / (B + 1)
  }
  answer
}

run_dependence_study <- function() {
  message("Running dependence calibration study")
  run_one_regime <- function(n_rep, change, seed_offset) {
    parallel_replicates(n_rep, seed_offset, function(i) {
      values <- resampling_pvalues(
        generate_ar1(change = change),
        B = cfg$dependence_B
      )
      as.data.table(as.list(values))
    })
  }

  null <- run_one_regime(cfg$dependence_null, 0, 7000000L)
  null[, scenario := "AR(1) null"]
  alternative <- run_one_regime(cfg$dependence_alt, 0.85, 8000000L)
  alternative[, scenario := "AR(1) mean change"]
  raw <- rbind(null, alternative)
  raw[, abs_error := abs(location - 60L)]
  raw[scenario == "AR(1) null", abs_error := NA_integer_]

  p_value_columns <- grep("^p_", names(raw), value = TRUE)
  long <- melt(
    raw,
    id.vars = c("replicate", "scenario", "statistic", "location", "abs_error"),
    measure.vars = p_value_columns,
    variable.name = "calibration",
    value.name = "p_value"
  )
  calibration_levels <- c("p_iid", paste0("p_block_", c(5, 10, 15, 20)))
  calibration_labels <- c("IID permutation", paste0("Block permutation (length ",
                                                       c(5, 10, 15, 20), ")"))
  long[, calibration := factor(calibration,
                               levels = calibration_levels,
                               labels = calibration_labels)]
  long[, reject := p_value <= alpha]

  summary <- long[, {
    successes <- sum(reject)
    ci <- wilson_interval(successes, .N)
    list(
      rejection = successes / .N,
      ci_low = ci[1], ci_high = ci[2], n_rep = .N,
      miss_rate = mean(!reject),
      mean_capped_error = mean(fifelse(reject, abs_error, 120)),
      median_error_detected = if (successes) as.numeric(median(abs_error[reject])) else NA_real_
    )
  }, by = .(scenario, calibration)]
  summary[scenario == "AR(1) null", `:=`(miss_rate = NA_real_, mean_capped_error = NA_real_)]

  list(raw = raw, long = long, summary = summary)
}

sym_log_features <- function(s, eigenvalues_only = FALSE) {
  eig <- eigen(s, symmetric = TRUE)
  log_values <- log(pmax(eig$values, 1e-10))
  if (eigenvalues_only) return(sort(log_values, decreasing = TRUE))
  log_matrix <- eig$vectors %*% diag(log_values) %*% t(eig$vectors)
  c(diag(log_matrix), sqrt(2) * log_matrix[lower.tri(log_matrix)])
}

simulate_spd_features <- function(n = 80L, internal_n, scenario, tau = 40L) {
  q <- 3L
  sigma0 <- diag(c(2, 1, 0.5))
  rotation <- diag(q)
  theta <- pi / 8
  rotation[1:2, 1:2] <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  sigma_rotation <- rotation %*% sigma0 %*% t(rotation)
  sigma_spectrum <- diag(c(2.5, 0.85, 0.40))

  full <- matrix(NA_real_, n, q * (q + 1L) / 2L)
  spectrum <- matrix(NA_real_, n, q)
  for (i in seq_len(n)) {
    sigma <- sigma0
    if (i > tau && scenario == "rotation") sigma <- sigma_rotation
    if (i > tau && scenario == "spectrum") sigma <- sigma_spectrum
    z <- matrix(rnorm(internal_n * q), nrow = internal_n, ncol = q)
    y <- z %*% chol(sigma)
    estimate <- crossprod(scale(y, center = TRUE, scale = FALSE)) / (internal_n - 1)
    estimate <- estimate + diag(1e-8, q)
    full[i, ] <- sym_log_features(estimate, FALSE)
    spectrum[i, ] <- sym_log_features(estimate, TRUE)
  }
  list(full = full, spectrum = spectrum)
}

compute_object_statistics <- function(features, cuts = seq.int(16L, 64L)) {
  one_representation <- function(x, label) {
    distance <- sqrt(pairwise_sqdist(x))
    result <- extract_scan(scan_matrix_cpp(distance, cuts, 1L), cuts, 40L)
    data.table(
      method = label,
      statistic = result["statistic"],
      location = as.integer(result["location"]),
      fixed = result["fixed"]
    )
  }
  rbind(
    one_representation(features$full, "Full log-Euclidean"),
    one_representation(features$spectrum, "Eigenvalue-only")
  )
}

run_object_study <- function() {
  message("Running SPD object-valued study")
  internal_sizes <- c(12L, 40L)
  alternatives <- c("spectrum", "rotation")
  all_raw <- list()
  threshold_rows <- list()
  seed_cursor <- 9000000L

  for (m in internal_sizes) {
    message("  internal sample size = ", m, ": calibration null")
    cal <- parallel_replicates(cfg$object_cal, seed_cursor, function(i) {
      compute_object_statistics(simulate_spd_features(internal_n = m, scenario = "null"))
    })
    seed_cursor <- seed_cursor + 100000L
    cal[, `:=`(phase = "calibration", scenario = "null", internal_n = m)]
    thresholds <- cal[, .(
      threshold = as.numeric(quantile(statistic, 1 - alpha, type = 8))
    ), by = method]
    thresholds[, internal_n := m]
    threshold_rows[[length(threshold_rows) + 1L]] <- thresholds
    all_raw[[length(all_raw) + 1L]] <- cal

    message("  internal sample size = ", m, ": evaluation null")
    evaluation <- parallel_replicates(cfg$object_eval, seed_cursor, function(i) {
      compute_object_statistics(simulate_spd_features(internal_n = m, scenario = "null"))
    })
    seed_cursor <- seed_cursor + 100000L
    evaluation[, `:=`(phase = "evaluation", scenario = "null", internal_n = m)]
    all_raw[[length(all_raw) + 1L]] <- evaluation

    for (scenario_name in alternatives) {
      message("  internal sample size = ", m, ": ", scenario_name)
      alternative <- parallel_replicates(cfg$object_alt, seed_cursor, function(i) {
        compute_object_statistics(
          simulate_spd_features(internal_n = m, scenario = scenario_name)
        )
      })
      seed_cursor <- seed_cursor + 100000L
      alternative[, `:=`(
        phase = "alternative", scenario = scenario_name, internal_n = m
      )]
      all_raw[[length(all_raw) + 1L]] <- alternative
    }
  }

  raw <- rbindlist(all_raw)
  thresholds <- rbindlist(threshold_rows)
  raw[thresholds, threshold := i.threshold, on = .(internal_n, method)]
  raw[, reject := statistic > threshold]
  raw[, abs_error := abs(location - 40L)]
  raw[, capped_error := fifelse(reject, abs_error, 80)]
  raw[scenario == "null", `:=`(abs_error = NA_integer_, capped_error = NA_real_)]

  summary <- raw[phase != "calibration", {
    successes <- sum(reject)
    ci <- wilson_interval(successes, .N)
    list(
      rejection = successes / .N,
      ci_low = ci[1], ci_high = ci[2], n_rep = .N,
      miss_rate = mean(!reject),
      mean_capped_error = mean(capped_error),
      median_error_detected = if (successes) as.numeric(median(abs_error[reject])) else NA_real_
    )
  }, by = .(scenario, internal_n, method)]
  summary[scenario == "null", miss_rate := NA_real_]

  list(raw = raw, thresholds = thresholds, summary = summary)
}

standalone_method <- function(x, method, cuts) {
  n <- nrow(x)
  p <- ncol(x)
  if (method %in% c("Max-CUSUM", "Dense projection")) {
    path <- cusum_paths(x, cuts)
    if (method == "Max-CUSUM") max(abs(path)) else max(abs(rowSums(path) / sqrt(p)))
  } else {
    d2 <- pairwise_sqdist(x)
    if (method == "Energy") {
      scan_matrix_cpp(sqrt(d2), cuts, 1L)$statistic
    } else if (method == "Gaussian MMD") {
      bandwidth2 <- median(d2[upper.tri(d2) & d2 > 0])
      scan_matrix_cpp(exp(-d2 / (2 * bandwidth2)), cuts, 2L)$statistic
    } else if (method == "5-NN graph") {
      scan_matrix_cpp(knn_graph(d2, 5L), cuts, 3L)$statistic
    } else stop("unknown method")
  }
}

run_runtime_study <- function() {
  message("Running implementation timing study")
  methods <- c("Max-CUSUM", "Dense projection", "Energy", "Gaussian MMD", "5-NN graph")
  dimensions <- c(5L, 50L, 200L)
  rows <- list()
  counter <- 0L
  for (p in dimensions) {
    cuts <- seq.int(24L, 96L)
    for (method in methods) {
      times <- numeric(cfg$runtime_reps)
      for (r in seq_len(cfg$runtime_reps)) {
        set.seed(master_seed + 12000000L + 1000L * p + 10L * r + match(method, methods))
        x <- generate_amoc(p = p, scenario = "null")
        inner_reps <- if (method %in% c("Max-CUSUM", "Dense projection")) 200L else 50L
        elapsed <- system.time({
          for (j in seq_len(inner_reps)) standalone_method(x, method, cuts)
        })["elapsed"]
        times[r] <- elapsed / inner_reps
      }
      counter <- counter + 1L
      rows[[counter]] <- data.table(
        p = p, method = method,
        median_milliseconds = 1000 * median(times),
        q25_milliseconds = 1000 * quantile(times, 0.25),
        q75_milliseconds = 1000 * quantile(times, 0.75),
        n_rep = length(times)
      )
    }
  }
  rbindlist(rows)
}

save_plot_pair <- function(plot_object, stem, width, height) {
  ggsave(file.path("figures", paste0(stem, suffix, ".pdf")), plot_object,
         width = width, height = height, units = "in", device = "pdf",
         useDingbats = FALSE)
  ggsave(file.path("figures", paste0(stem, suffix, ".png")), plot_object,
         width = width, height = height, units = "in", dpi = 220)
}

make_figures <- function(amoc, dependence, object) {
  alternative_order <- c("sparse", "dense", "orthogonal dense", "scale",
                         "correlation", "radial mixture")
  method_order <- c("Max-CUSUM", "Dense projection", "Energy", "Gaussian MMD", "5-NN graph")
  method_labels <- setNames(method_order, method_order)
  method_labels["5-NN graph"] <- "Cross-edge scan"
  power <- copy(amoc$summary[scenario %in% alternative_order])
  power[, scenario := factor(scenario, levels = alternative_order,
                             labels = c("Sparse mean", "Dense mean", "Dense mean, reflected", "Scale",
                                        "Correlation", "Radial mixture"))]
  power[, method := factor(method, levels = method_order)]
  setorder(power, p, scenario, method)
  power[, row_label := paste0(scenario, "; p=", p)]
  power[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  p_power <- ggplot(power, aes(method, row_label, fill = rejection)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.3f", rejection),
                  colour = ifelse(rejection < .50, "white", "black")), size = 2.6) +
    scale_colour_identity() +
    scale_fill_viridis_c(limits = c(0, 1), option = "C", name = "Power") +
    scale_x_discrete(labels = method_labels) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid = element_blank())
  save_plot_pair(p_power, "fig_representation_power", 7.2, 5.4)

  calibration <- copy(amoc$calibration)
  calibration[, method := factor(method, levels = method_order)]
  p_scan <- ggplot(calibration, aes(calibration, estimate, color = calibration)) +
    geom_hline(yintercept = alpha, linetype = 2, color = "grey35") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(0.4)) +
    facet_grid(p ~ method, labeller = labeller(method = method_labels)) +
    coord_cartesian(ylim = c(0, max(0.35, calibration$ci_high))) +
    labs(x = NULL, y = "Type I error rate", color = "Calibration",
         title = "(a) Fixed-split and full-scan thresholds") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", strip.text.x = element_text(size = 7),
          plot.title = element_text(face = "bold"))

  dep_null <- copy(dependence$summary[scenario == "AR(1) null"])
  p_dep <- ggplot(dep_null, aes(calibration, rejection, color = calibration)) +
    geom_hline(yintercept = alpha, linetype = 2, color = "grey35") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), size = 0.5) +
    coord_cartesian(ylim = c(0, max(0.35, dep_null$ci_high))) +
    scale_x_discrete(labels = function(labels) {
      ifelse(labels == "IID permutation", "Individual permutation", labels)
    }) +
    labs(x = NULL, y = "Type I error rate",
         title = "(b) Individual and block permutations under an AR(1) null") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none",
          plot.title = element_text(face = "bold"))

  pdf(NULL)
  combined <- gridExtra::arrangeGrob(p_scan, p_dep, ncol = 1, heights = c(1.5, 1))
  dev.off()
  pdf(file.path("figures", paste0("fig_calibration", suffix, ".pdf")),
      width = 7.2, height = 7.0, useDingbats = FALSE)
  grid::grid.draw(combined)
  dev.off()
  png(file.path("figures", paste0("fig_calibration", suffix, ".png")),
      width = 7.2, height = 7.0, units = "in", res = 220)
  grid::grid.draw(combined)
  dev.off()

  object_power <- copy(object$summary[scenario != "null"])
  object_power[, scenario := factor(scenario, levels = c("spectrum", "rotation"),
                                    labels = c("Eigenvalue change", "Rotation only"))]
  p_object <- ggplot(object_power,
                     aes(internal_n, rejection, color = method, shape = method)) +
    geom_hline(yintercept = alpha, linetype = 3, color = "grey55") +
    geom_line(linewidth = 0.55) +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), size = 0.45) +
    facet_wrap(~ scenario) +
    scale_x_continuous(breaks = c(12, 40)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Observations used to estimate each SPD matrix", y = "Estimated rejection probability",
         color = "Representation", shape = "Representation") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  save_plot_pair(p_object, "fig_object_power", 7.2, 3.8)
}

if (plot_only) {
  saved <- readRDS(file.path("analysis", "results", paste0("diagnostic_benchmark", suffix, ".rds")))
  make_figures(saved$amoc, saved$dependence, saved$object)
  message("Figures regenerated from saved ", mode, " results.")
  quit(save = "no", status = 0)
}

amoc <- run_amoc_study()
dependence <- run_dependence_study()
object <- run_object_study()
runtime <- run_runtime_study()

result_dir <- file.path("analysis", "results")
fwrite(amoc$summary, file.path(result_dir, paste0("amoc_summary", suffix, ".csv")))
fwrite(amoc$calibration, file.path(result_dir, paste0("scan_calibration", suffix, ".csv")))
fwrite(amoc$thresholds, file.path(result_dir, paste0("amoc_thresholds", suffix, ".csv")))
fwrite(dependence$summary, file.path(result_dir, paste0("dependence_summary", suffix, ".csv")))
fwrite(object$summary, file.path(result_dir, paste0("object_summary", suffix, ".csv")))
fwrite(object$thresholds, file.path(result_dir, paste0("object_thresholds", suffix, ".csv")))
fwrite(runtime, file.path(result_dir, paste0("runtime_summary", suffix, ".csv")))
saveRDS(
  list(config = cfg, amoc = amoc, dependence = dependence, object = object,
       runtime = runtime, mode = mode, master_seed = master_seed),
  file.path(result_dir, paste0("diagnostic_benchmark", suffix, ".rds")),
  compress = "xz"
)

make_figures(amoc, dependence, object)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
metadata <- c(
  paste("mode:", mode),
  paste("started:", format(start_time, tz = "Asia/Seoul")),
  paste("elapsed_seconds:", sprintf("%.2f", elapsed)),
  paste("master_seed:", master_seed),
  paste("workers:", cores),
  paste("R:", R.version.string),
  paste("platform:", R.version$platform),
  paste("system:", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " ")),
  paste("data.table:", as.character(packageVersion("data.table"))),
  paste("Rcpp:", as.character(packageVersion("Rcpp"))),
  paste("ggplot2:", as.character(packageVersion("ggplot2"))),
  paste("configuration:", paste(names(cfg), unlist(cfg), sep = "=", collapse = "; "))
)
writeLines(metadata, file.path(result_dir, paste0("run_metadata", suffix, ".txt")))

message(sprintf("Completed %s run in %.1f seconds", mode, elapsed))
