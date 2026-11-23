#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rcpp)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(root_dir)
sourceCpp(file.path("analysis", "benchmark_core.cpp"))

set.seed(104729)
x <- matrix(rnorm(24 * 4), 24, 4)
sq <- rowSums(x * x)
d2 <- pmax(outer(sq, sq, "+") - 2 * tcrossprod(x), 0)
d <- sqrt(d2)
bw2 <- median(d2[upper.tri(d2)])
k <- exp(-d2 / (2 * bw2))
cuts <- 5:19

energy_cpp <- scan_matrix_cpp(d, cuts, 1L)$values
mmd_cpp <- scan_matrix_cpp(k, cuts, 2L)$values

energy_r <- vapply(cuts, function(t) {
  left <- seq_len(t)
  right <- (t + 1):nrow(x)
  between <- mean(d[left, right, drop = FALSE])
  within_left <- (sum(d[left, left, drop = FALSE]) - sum(diag(d)[left])) /
    (t * (t - 1))
  nr <- length(right)
  within_right <- (sum(d[right, right, drop = FALSE]) - sum(diag(d)[right])) /
    (nr * (nr - 1))
  t * nr / nrow(x) * (2 * between - within_left - within_right)
}, numeric(1))

mmd_r <- vapply(cuts, function(t) {
  left <- seq_len(t)
  right <- (t + 1):nrow(x)
  nr <- length(right)
  within_left <- (sum(k[left, left, drop = FALSE]) - sum(diag(k)[left])) /
    (t * (t - 1))
  within_right <- (sum(k[right, right, drop = FALSE]) - sum(diag(k)[right])) /
    (nr * (nr - 1))
  between <- mean(k[left, right, drop = FALSE])
  t * nr / nrow(x) * (within_left + within_right - 2 * between)
}, numeric(1))

stopifnot(max(abs(energy_cpp - energy_r)) < 1e-10)
stopifnot(max(abs(mmd_cpp - mmd_r)) < 1e-10)

# Check the graph cross-edge mean and variance against complete enumeration.
n <- 8L
t <- 3L
adj <- matrix(0, n, n)
edge_list <- rbind(c(1, 2), c(1, 3), c(2, 4), c(3, 5), c(4, 5),
                   c(4, 6), c(5, 7), c(6, 8), c(7, 8))
for (i in seq_len(nrow(edge_list))) {
  u <- edge_list[i, 1]
  v <- edge_list[i, 2]
  adj[u, v] <- adj[v, u] <- 1
}
groups <- combn(n, t)
cross_counts <- apply(groups, 2, function(left) {
  right <- setdiff(seq_len(n), left)
  sum(adj[left, right, drop = FALSE])
})
edges <- sum(adj) / 2
degree <- rowSums(adj)
adjacent_pairs <- sum(degree * (degree - 1) / 2)
all_pairs <- edges * (edges - 1) / 2
disjoint_pairs <- all_pairs - adjacent_pairs
p_cross <- 2 * t * (n - t) / (n * (n - 1))
q_adjacent <- t * (n - t) / (n * (n - 1))
q_disjoint <- 4 * t * (t - 1) * (n - t) * (n - t - 1) /
  (n * (n - 1) * (n - 2) * (n - 3))
variance_formula <- edges * p_cross * (1 - p_cross) +
  2 * adjacent_pairs * (q_adjacent - p_cross^2) +
  2 * disjoint_pairs * (q_disjoint - p_cross^2)
stopifnot(abs(mean(cross_counts) - edges * p_cross) < 1e-12)
stopifnot(abs(var(cross_counts) * (length(cross_counts) - 1) / length(cross_counts) -
              variance_formula) < 1e-12)
graph_cpp <- apply(groups, 2, function(left) {
  order <- c(left, setdiff(seq_len(n), left))
  scan_matrix_cpp(adj[order, order], t, 3L)$statistic
})
stopifnot(max(abs(graph_cpp - (mean(cross_counts) - cross_counts) /
                   sqrt(mean((cross_counts - mean(cross_counts))^2)))) < 1e-10)

archive <- file.path("analysis", "results", "diagnostic_benchmark.rds")
if (file.exists(archive)) {
  library(data.table)
  result <- readRDS(archive)
  for (module in c("amoc", "object", "dependence")) {
    null <- result[[module]]$summary[grepl("null", scenario)]
    stopifnot(all(is.na(null$miss_rate)), all(is.na(null$mean_capped_error)),
              all(is.na(null$median_error_detected)))
  }
  paired <- dcast(result$amoc$raw[p == 200 & scenario %in% c("dense", "orthogonal dense")],
                  method + replicate ~ scenario, value.var = "statistic")
  stopifnot(nrow(paired) == 5L * result$config$amoc_alt,
            paired[method != "Dense projection", max(abs(dense - `orthogonal dense`))] < 1e-9)
}

message("All benchmark statistic checks passed.")
