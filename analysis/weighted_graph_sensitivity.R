#!/usr/bin/env Rscript

# One-cell sensitivity to the graph score (Chu and Chen, Ann. Statist. 2019,
# Sections 3.1--3.2 and Lemma 3.1, https://doi.org/10.1214/18-AOS1691; weighted two-sample count:
# Chen, Chen and Su, https://doi.org/10.1080/01621459.2017.1307757).
suppressPackageStartupMessages({library(data.table); library(Rcpp)})
script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", script[1])))
sourceCpp(file.path(here, "benchmark_core.cpp"))
out <- file.path(here, "results")
archive <- readRDS(file.path(out, "diagnostic_benchmark.rds"))
start_time <- Sys.time()

weighted_moments <- function(n, t, edges, degree) {
  r <- n-t
  # Exact permutation moments following (3.4) in Chu and Chen (2019).
  graph_factor <- edges - sum(degree^2)/(n-2) + 2*edges^2/((n-1)*(n-2))
  list(mean = edges*(t-1)*(r-1)/((n-1)*(n-2)),
       var = t*(t-1)*r*(r-1)/(n*(n-1)*(n-2)*(n-3))*graph_factor)
}

difference_moments <- function(n, t, edges, degree) {
  # R1 - R2 equals the sum of the t left-vertex degrees minus |G|.
  list(mean = (2*t/n-1)*edges, var = t*(n-t)/(n*(n-1))*sum((degree-2*edges/n)^2))
}

# Verify permutation moments, zero covariance of the two components, and
# equality to the full Mahalanobis statistic for every small-graph allocation.
test_edges <- rbind(c(1,2), c(1,3), c(2,4), c(3,5), c(4,5), c(4,6),
                    c(5,7), c(6,8), c(7,8))
for (t in 2:6) {
  within <- t(combn(8, t, function(left) {
    groups <- matrix(test_edges %in% left, ncol = 2)
    counts <- rowSums(groups)
    c(sum(counts == 2), sum(counts == 0))
  }))
  observed <- ((7-t)*within[,1] + (t-1)*within[,2])/6
  difference <- within[,1]-within[,2]
  expected <- weighted_moments(8, t, nrow(test_edges), tabulate(test_edges, 8))
  dm <- difference_moments(8, t, nrow(test_edges), tabulate(test_edges, 8))
  centered <- sweep(within, 2, colMeans(within))
  sigma <- crossprod(centered)/nrow(within)
  direct <- rowSums((centered %*% solve(sigma))*centered)
  components <- (observed-expected$mean)^2/expected$var + (difference-dm$mean)^2/dm$var
  stopifnot(abs(mean(observed) - expected$mean) < 1e-12,
            abs(mean((observed-mean(observed))^2) - expected$var) < 1e-12,
            abs(mean(difference) - dm$mean) < 1e-12,
            abs(mean((difference-dm$mean)^2) - dm$var) < 1e-12,
            abs(cov(observed, difference)) < 1e-12,
            max(abs(components-direct)) < 1e-10)
}

run_one <- function(seed, alternative) {
  set.seed(seed)
  n <- 120L
  cuts <- 24:96
  x <- matrix(rnorm(n*200), n)
  if (alternative) {
    s <- ifelse(runif(60) < .8, .65, sqrt((1-.8*.65^2)/.2))
    x[61:120, ] <- x[61:120, ] * s
  }
  sq <- rowSums(x*x)
  d2 <- pmax(outer(sq, sq, "+")-2*tcrossprod(x), 0)
  diag(d2) <- Inf
  nearest <- apply(d2, 1, function(row) order(row)[1:5])
  graph <- matrix(0, n, n)
  graph[cbind(rep(seq_len(n), each = 5), as.vector(nearest))] <- 1
  graph <- 1*((graph+t(graph)) > 0)
  uv <- which(upper.tri(graph) & graph > 0, arr.ind = TRUE)
  m <- nrow(uv)
  degree <- tabulate(uv, n)
  left <- cumsum(tabulate(uv[,2], n))[cuts]
  right <- m-cumsum(tabulate(uv[,1], n))[cuts]
  moments <- weighted_moments(n, cuts, m, degree)
  dm <- difference_moments(n, cuts, m, degree)
  stopifnot(all(moments$var > 0), all(dm$var > 0))
  values <- (((n-cuts-1)*left + (cuts-1)*right)/(n-2) -
             moments$mean)/sqrt(moments$var)
  z_difference <- (left-right-dm$mean)/sqrt(dm$var)
  generalized <- values^2 + z_difference^2
  original <- scan_matrix_cpp(graph, cuts, 3L)
  c(statistic = max(values), location = cuts[which.max(values)],
    fixed = values[cuts == 60], z_difference_fixed = z_difference[cuts == 60],
    generalized_stat = max(generalized), generalized_cut = cuts[which.max(generalized)],
    generalized_fixed = generalized[cuts == 60], original = original$statistic,
    original_cut = original$location)
}

seeds <- list(calibration = 1400000L + seq_len(10000),
              evaluation = 1500000L + seq_len(2000),
              alternative = 2000000L + seq_len(500))
workers <- as.integer(Sys.getenv("CPD_CORES", "4"))
weighted <- rbindlist(lapply(names(seeds), function(phase_name) {
  values <- parallel::mclapply(archive$master_seed + seeds[[phase_name]],
    function(seed) run_one(seed, phase_name == "alternative"),
    mc.cores = workers, mc.set.seed = FALSE)
  result <- as.data.table(do.call(rbind, values))
  result[, `:=`(phase = phase_name, replicate = seq_len(.N))]
  result
}))
original <- archive$amoc$raw[p == 200 & method == "5-NN graph" &
                               scenario %in% c("null", "radial mixture")]
paired <- merge(weighted, original[, .(phase, replicate, original_stat = statistic,
                                        old_cut = location, original_fixed = fixed)],
                by = c("phase", "replicate"))
stopifnot(nrow(paired) == 12500L, max(abs(paired$original-paired$original_stat)) < 1e-8,
          all(paired$original_cut == paired$old_cut),
          max(abs(paired$fixed-paired$original_fixed)) < 1e-8)
values <- rbind(weighted[, .(phase, replicate, method = "Weighted within-edge",
                             statistic, location, fixed)],
                weighted[, .(phase, replicate, method = "Generalized edge-count",
                             statistic = generalized_stat, location = generalized_cut,
                             fixed = generalized_fixed)],
                original[, .(phase, replicate, method = "Cross-edge deficit",
                             statistic, location, fixed)])
values[, threshold := as.numeric(quantile(statistic[phase == "calibration"],
                                          .95, type = 8)), by = method]
values[, reject := statistic > threshold]
summary <- values[phase != "calibration", {
  rate <- mean(reject)
  z <- qnorm(.975)
  center <- (rate+z^2/(2*.N))/(1+z^2/.N)
  half <- z*sqrt(rate*(1-rate)/.N+z^2/(4*.N^2))/(1+z^2/.N)
  locations <- location[phase == "alternative" & reject]
  .(n_rep = .N, detections = sum(reject), rejection = rate,
    ci_low = center-half, ci_high = center+half,
    threshold = threshold[1], median_detected_location = as.numeric(median(locations)),
    median_detected_error = as.numeric(median(abs(locations-60))),
    median_score_at_60 = median(fixed))
}, by = .(method, phase)]
fwrite(summary, file.path(out, "weighted_graph_sensitivity.csv"))
saveRDS(list(scores = values, difference_at_60 = weighted[, .(phase, replicate,
              z_difference_fixed)]),
        file.path(out, "weighted_graph_sensitivity.rds"), compress = "xz")
print(summary)
message(sprintf("Elapsed %.1f seconds on %i workers",
                as.numeric(difftime(Sys.time(), start_time, units = "secs")), workers))
