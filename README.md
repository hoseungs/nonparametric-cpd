# Nonparametric change-point detection for modern data

This repository provides the R codes for the paper:

> **Nonparametric Change-Point Detection for Modern Data: Representations and Statistical Guarantees.**

> Hoseung Song [arxiv](https://arxiv.org/pdf/2611.32158)

Section 5 compares how scans are calibrated, which changes they detect, and where they place the change. This repository includes the saved simulation results and the code for Figures 1–3 and the covariance-matrix experiment. 

## Find a result

| In the paper | Result | Files |
| --- | --- | --- |
| Section 5.2, Figure 1 | Fixed-split versus full-scan calibration; individual versus block permutations for an AR(1) sequence | [Figure 1](figures/fig_calibration.pdf), [scan rejection rates](analysis/results/scan_calibration.csv), [permutation results](analysis/results/dependence_summary.csv) |
| Section 5.3, Figure 2 | Power of five scans under sparse and dense mean shifts, correlation and scale changes, and a radial mixture | [Figure 2](figures/fig_representation_power.pdf), [power and location errors](analysis/results/amoc_summary.csv) |
| Section 5.3, Figure 3 | Detection and estimated locations on the same 500 radial-mixture sequences using three scores on a common graph and an energy scan | [Figure 3](figures/fig_radial_localization.pdf), [graph-score comparison](analysis/results/weighted_graph_sensitivity.csv) |
| Section 5.4 | Spectral and orientation changes in sequences of estimated covariance matrices | [rejection rates and location errors](analysis/results/object_summary.csv), [critical values](analysis/results/object_thresholds.csv) |

## Quick example

Download or clone the repository and open R in its top-level directory. The summary tables can be read with base R. For instance, a rotation changes the orientation of a covariance matrix while preserving its eigenvalues:

```r
objects <- read.csv("analysis/results/object_summary.csv")
subset(objects, scenario == "rotation", select = c(internal_n, method, rejection))
```

With 12 measurements per matrix, the full log-Euclidean scan detected the rotation in 72.4% of sequences, compared with 4.6% for the eigenvalue-only scan. With 40 measurements, the corresponding rates were 100% and 3.2%.

## Reproduce the figures from the saved results

Install the R dependencies once:

```r
install.packages(c("data.table", "Rcpp", "ggplot2", "gridExtra"))
```

From the repository root, run:

```sh
CPD_MODE=full CPD_PLOT_ONLY=1 Rscript analysis/run_diagnostic_benchmark.R
Rscript analysis/plot_radial_localization.R
Rscript analysis/check_reported_contrasts.R
```

The first command draws Figures 1 and 2 and an additional covariance-matrix plot. The second draws Figure 3 using both saved simulation archives. Plots are written to `figures/` as PDFs and PNGs. The third command prints the fixed-split rejection range, the radial-mixture estimates, and the paired contrast for the dense projection (443 detections out of 500 before reflection and 34 after reflection).

## Rerun the simulations

The simulation and statistic-check scripts compile [the scan implementation](analysis/benchmark_core.cpp) with Rcpp, so they require a C++ toolchain for R. Run a numerical check and a small pilot first:

```sh
Rscript analysis/test_benchmark.R
CPD_MODE=pilot CPD_CORES=4 Rscript analysis/run_diagnostic_benchmark.R
```

Pilot files are marked `_pilot`. To reproduce the paper's configuration from scratch, run the following commands in order. They write new full-run results to `analysis/results/`:

```sh
CPD_MODE=full CPD_CORES=4 Rscript analysis/run_diagnostic_benchmark.R
CPD_CORES=4 Rscript analysis/weighted_graph_sensitivity.R
Rscript analysis/plot_radial_localization.R
```

The second command reconstructs the 12,500 graphs in the radial-mixture comparison from the original seeds: 10,000 null sequences to set thresholds, 2,000 independent null sequences to evaluate them, and 500 alternatives. It recalibrates the weighted and generalized edge-count scans separately and checks the original scan against the main archive.

The shell examples use POSIX syntax. On Windows, set `CPD_CORES=1` for the simulation scripts. In PowerShell, set variables before calling `Rscript`, for example `$env:CPD_CORES = "1"` and `$env:CPD_MODE = "full"`. Set `CPD_PLOT_ONLY` to `0` for simulations if it was previously set to `1` for plotting.

## Simulation settings

The main run used master seed `20261118`, R 4.4.2, and four workers. The recorded elapsed time for the main benchmark was 101 seconds on an ARM64 Mac; the separate graph-score comparison takes additional time. See [the run metadata](analysis/results/run_metadata.txt) for the recorded configuration and package versions.

| Experiment | Design |
| --- | --- |
| Euclidean scans | Length 120, change after observation 60, candidate splits 24–96, dimensions 5, 50 and 200. Each dimension uses 10,000 training nulls, 2,000 evaluation nulls, and 500 sequences per alternative. |
| Serial dependence | Gaussian AR(1) with coefficient 0.5, 1,000 null sequences, 500 alternatives, and 199 rearrangements per sequence and permutation scheme. Block lengths are 5, 10, 15 and 20. |
| Covariance matrices | Length 80, change after observation 40, 12 or 40 measurements per matrix. Each measurement count uses 2,000 training nulls, 2,000 evaluation nulls, and 500 sequences per alternative. |

The [main archive](analysis/results/diagnostic_benchmark.rds) contains replicate-level results, thresholds, configuration and summaries. The [graph archive](analysis/results/weighted_graph_sensitivity.rds) contains the paired graph scores and estimated splits. Both can be opened with `readRDS()`. The scans share simulated sequences, which is why comparisons such as the reflection experiment can be made within replicate.
