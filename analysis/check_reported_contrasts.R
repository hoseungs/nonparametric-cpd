#!/usr/bin/env Rscript

# Reconstruct the manuscript's paired and signed summaries from raw scans.
suppressPackageStartupMessages(library(data.table))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
amoc <- readRDS(file.path(script_dir, "results", "diagnostic_benchmark.rds"))$amoc$raw

fixed <- amoc[phase == "evaluation",
              .(rate = mean(fixed > fixed_threshold)), by = .(p, method)]
cat(sprintf("Held-out fixed-cut rejection across %d cells: %.4f to %.4f\n",
            nrow(fixed), min(fixed$rate), max(fixed$rate)))

radial <- amoc[p == 200 & scenario == "radial mixture"]
g <- radial[method == "5-NN graph"]
ci <- quantile(g[reject == TRUE, location], c(.25, .75), names = FALSE)
cat(sprintf(paste0("Graph: %d/%d detected; %d/%d scores below zero at the true cut",
                   " (median %.3f)\n"), sum(g$reject), nrow(g), sum(g$fixed < 0),
            nrow(g), median(g$fixed)))
cat(sprintf(paste0("Detected graph peak: median %.0f, IQR %.0f to %.0f; ",
                   "energy peak: median %.0f\n"),
            g[reject == TRUE, median(location)], ci[1], ci[2],
            radial[method == "Energy" & reject == TRUE, median(location)]))

aligned <- amoc[p == 200 & scenario == "dense", .(method, replicate,
                aligned_stat = statistic, aligned_reject = reject)]
reflected <- amoc[scenario == "orthogonal dense", .(method, replicate,
                  reflected_stat = statistic, reflected_reject = reject)]
pair_keys <- c("method", "replicate")
stopifnot(!anyDuplicated(aligned, by = pair_keys), !anyDuplicated(reflected, by = pair_keys))
pairs <- merge(aligned, reflected, by = pair_keys, all = TRUE)
stopifnot(!anyNA(pairs))
invariant <- pairs[method != "Dense projection",
                   .(max_error = max(abs(aligned_stat - reflected_stat)),
                     discordant = sum(aligned_reject != reflected_reject)),
                   by = method]
stopifnot(all(invariant$max_error < 1e-9), all(invariant$discordant == 0))
cat(sprintf("Four invariant scans: largest paired statistic difference %.3g; %d discordant decisions\n",
            max(invariant$max_error), sum(invariant$discordant)))
dense <- pairs[method == "Dense projection"]
delta <- as.integer(dense$aligned_reject) - as.integer(dense$reflected_reject)
cat(sprintf(paste0("Fixed projection: %d/%d aligned, %d/%d reflected; ",
                   "paired difference %.3f, paired MCSE %.4f\n"),
            sum(dense$aligned_reject), nrow(dense), sum(dense$reflected_reject),
            nrow(dense), mean(delta), sd(delta) / sqrt(length(delta))))
