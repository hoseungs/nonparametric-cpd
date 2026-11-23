#include <Rcpp.h>
using namespace Rcpp;

// Rectangle sums for a matrix whose rows and columns stay in temporal order.
// The prefix matrix makes a complete scan O(n^2 + |T|), rather than O(n^3).

// [[Rcpp::export]]
List scan_matrix_cpp(NumericMatrix a, IntegerVector cuts, int statistic_type) {
  const int n = a.nrow();
  if (a.ncol() != n) stop("a must be square");

  NumericMatrix prefix(n + 1, n + 1);
  for (int i = 1; i <= n; ++i) {
    double row_sum = 0.0;
    for (int j = 1; j <= n; ++j) {
      row_sum += a(i - 1, j - 1);
      prefix(i, j) = prefix(i - 1, j) + row_sum;
    }
  }

  const double total = prefix(n, n);
  double edges = 0.0, adjacent_pairs = 0.0, disjoint_pairs = 0.0;
  if (statistic_type == 3) {
    if (n < 4) stop("graph standardization requires at least four observations");
    edges = total / 2.0;
    for (int i = 0; i < n; ++i) {
      const double degree = prefix(i+1, n) - prefix(i, n);
      adjacent_pairs += degree * (degree - 1.0) / 2.0;
    }
    disjoint_pairs = edges * (edges - 1.0) / 2.0 - adjacent_pairs;
  }
  NumericVector diagonal_prefix(n + 1);
  for (int i = 1; i <= n; ++i) {
    diagonal_prefix[i] = diagonal_prefix[i - 1] + a(i - 1, i - 1);
  }
  NumericVector values(cuts.size());
  double best = R_NegInf;
  int best_cut = NA_INTEGER;

  for (int h = 0; h < cuts.size(); ++h) {
    const int t = cuts[h];
    if (t <= 0 || t >= n) stop("cuts must lie between 1 and n - 1");
    const int r = n - t;
    if (statistic_type != 3 && (t < 2 || r < 2))
      stop("U-statistics require at least two observations per group");
    const double ll = prefix(t, t);
    const double lr = prefix(t, n) - prefix(t, t);
    const double rr = total - prefix(t, n) - prefix(n, t) + prefix(t, t);
    double value = NA_REAL;

    if (statistic_type == 1) {
      // Unbiased energy distance, multiplied by the usual split-balance factor.
      const double energy = 2.0 * lr / (static_cast<double>(t) * r) -
        (ll - diagonal_prefix[t]) / (static_cast<double>(t) * (t - 1.0)) -
        (rr - (diagonal_prefix[n] - diagonal_prefix[t])) /
          (static_cast<double>(r) * (r - 1.0));
      value = static_cast<double>(t) * r / n * energy;
    } else if (statistic_type == 2) {
      // Unbiased squared MMD, with the same split-balance factor.
      const double mmd2 =
        (ll - diagonal_prefix[t]) / (static_cast<double>(t) * (t - 1.0)) +
        (rr - (diagonal_prefix[n] - diagonal_prefix[t])) /
          (static_cast<double>(r) * (r - 1.0)) -
        2.0 * lr / (static_cast<double>(t) * r);
      value = static_cast<double>(t) * r / n * mmd2;
    } else if (statistic_type == 3) {
      // Permutation-standardized deficit of cross-group edges in a fixed
      // undirected similarity graph.  The variance accounts for pairs of
      // edges that share a vertex.
      const double cross_probability =
        2.0 * static_cast<double>(t) * r /
        (static_cast<double>(n) * (n - 1.0));
      const double joint_adjacent =
        static_cast<double>(t) * r /
        (static_cast<double>(n) * (n - 1.0));
      const double joint_disjoint =
        4.0 * static_cast<double>(t) * (t - 1.0) * r * (r - 1.0) /
        (static_cast<double>(n) * (n - 1.0) * (n - 2.0) * (n - 3.0));
      const double variance =
        edges * cross_probability * (1.0 - cross_probability) +
        2.0 * adjacent_pairs *
          (joint_adjacent - cross_probability * cross_probability) +
        2.0 * disjoint_pairs *
          (joint_disjoint - cross_probability * cross_probability);
      const double expected_cross = edges * cross_probability;
      value = (expected_cross - lr) / std::sqrt(std::max(variance, 1e-12));
    } else {
      stop("unknown statistic_type");
    }

    values[h] = value;
    if (R_finite(value) && value > best) {
      best = value;
      best_cut = t;
    }
  }

  return List::create(
    _["statistic"] = best,
    _["location"] = best_cut,
    _["values"] = values
  );
}
