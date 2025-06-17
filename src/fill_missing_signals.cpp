#include <Rcpp.h>

// [[Rcpp::export(name = .detect_missing_cpp)]]
Rcpp::IntegerVector detect_missing_cpp(Rcpp::NumericVector intensity,
                                       int gap = 1,
                                       double noise = 100) {
  std::vector<int> result;
  int n = intensity.size();

  for (int i = 0; i < n; ++i) {
    if (intensity[i] == 0) {
      // 检查是否在有效范围内
      if (i >= gap && i < n - gap) {
        bool left_ok = false;
        bool right_ok = false;

        // 检查左侧
        for (int j = 1; j <= gap; ++j) {
          if (intensity[i - j] > noise) {
            left_ok = true;
            break;
          }
        }

        // 检查右侧
        for (int j = 1; j <= gap; ++j) {
          if (intensity[i + j] > noise) {
            right_ok = true;
            break;
          }
        }

        if (left_ok && right_ok) {
          result.push_back(i + 1); // 只在最后转换为1-based
        }
      }
    }
  }

  return Rcpp::wrap(result);
}
