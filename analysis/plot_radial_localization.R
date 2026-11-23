#!/usr/bin/env Rscript

# Subdistribution of the reported split: all 500 sequences are in the denominator.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", script[1])))
results <- file.path(here, "results")
graphs <- readRDS(file.path(results, "weighted_graph_sensitivity.rds"))$scores
main <- readRDS(file.path(results, "diagnostic_benchmark.rds"))$amoc$raw
scores <- rbind(graphs[phase == "alternative", .(method, replicate, location, reject)],
                main[p == 200 & scenario == "radial mixture" & method == "Energy",
                     .(method, replicate, location, reject)])
display_names <- c("Cross-edge deficit" = "Cross-edge scan",
                   "Weighted within-edge" = "Weighted edge count",
                   "Generalized edge-count" = "Generalized edge count",
                   "Energy" = "Energy")
stopifnot(nrow(scores) == 2000L, !anyDuplicated(scores, by = c("method", "replicate")),
          all(scores[, .N, by = method]$N == 500L))
scores[, method := factor(method, levels = names(display_names), labels = display_names)]
curve <- scores[, {
  grid <- 24:96
  .(split = grid, fraction = cumsum(tabulate(match(location[reject], grid),
                                             nbins = length(grid)))/.N)
}, by = method]
plot <- ggplot(curve, aes(split, fraction, colour = method, linetype = method)) +
  geom_vline(xintercept = 60, colour = "grey45", linetype = "dashed", linewidth = .45) +
  geom_step(linewidth = .75) +
  scale_colour_manual(values = c("#0072B2", "#882255", "#00734A", "#BB4B00")) +
  scale_linetype_manual(values = c("solid", "22", "longdash", "dotdash")) +
  scale_y_continuous(breaks = seq(0, 1, .25), limits = c(-.025, 1.025),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(24, 96, 12), labels = c("24", "36", "48", "60 (true)", "72", "84", "96"),
                     expand = expansion(mult = .015)) +
  labs(x = expression("Candidate split" ~ italic(t)),
       y = expression(hat(P)(plain(reject) ~ plain(and) ~ hat(tau) <= italic(t))),
       colour = NULL, linetype = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(.55, "in"),
        panel.grid.minor = element_blank())
stem <- file.path(dirname(here), "figures", "fig_radial_localization")
ggsave(paste0(stem, ".pdf"), plot, width = 7.2, height = 3.5, useDingbats = FALSE)
ggsave(paste0(stem, ".png"), plot, width = 7.2, height = 3.5, dpi = 220)
