// spectral_denoising

#include <Rcpp.h>
#include <algorithm>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericMatrix electronic_denoising_cpp(Rcpp::NumericMatrix msms) {
  // 检查输入是否为NA或无效
  if (msms.nrow() == 0 || msms.ncol() != 2) {
    Rcpp::NumericMatrix empty(1, 2);
    empty(0, 0) = NA_REAL;
    empty(0, 1) = NA_REAL;
    return empty;
  }

  // 提取m/z和强度
  Rcpp::NumericVector mass = msms(Rcpp::_, 0);
  Rcpp::NumericVector intensity = msms(Rcpp::_, 1);

  // 创建索引向量并排序
  std::vector<size_t> idx(mass.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

  // 按强度升序排序
  std::sort(idx.begin(), idx.end(),
            [&intensity](size_t i1, size_t i2) {return intensity[i1] < intensity[i2];});

  // 重新排序mass和intensity
  Rcpp::NumericVector mass_sorted(mass.size());
  Rcpp::NumericVector intensity_sorted(intensity.size());
  for (size_t i = 0; i < idx.size(); ++i) {
    mass_sorted[i] = mass[idx[i]];
    intensity_sorted[i] = intensity[idx[i]];
  }

  // 确认的峰
  std::vector<double> mass_confirmed;
  std::vector<double> intensity_confirmed;

  while (intensity_sorted.size() > 0) {
    double seed_intensity = Rcpp::max(intensity_sorted);
    double threshold = seed_intensity * 0.999;

    // 找到第一个大于等于threshold的索引
    size_t idx_left = std::lower_bound(intensity_sorted.begin(),
                                       intensity_sorted.end(),
                                       threshold) - intensity_sorted.begin();

    // 确保idx_left不会导致无效范围
    size_t current_intensity_sorted_size = static_cast<size_t>(intensity_sorted.size());
    if (idx_left >= current_intensity_sorted_size) {
      break;
    }

    // 提取临时峰
    std::vector<double> mass_temp;
    std::vector<double> intensity_temp;

    size_t current_mass_sorted_size = static_cast<size_t>(mass_sorted.size());
    for (size_t i = idx_left; i < current_mass_sorted_size; ++i) {
      mass_temp.push_back(mass_sorted[i]);
      intensity_temp.push_back(intensity_sorted[i]);
    }

    // 如果峰数量<=3则保留
    if (mass_temp.size() <= 3) {
      mass_confirmed.insert(mass_confirmed.end(), mass_temp.begin(), mass_temp.end());
      intensity_confirmed.insert(intensity_confirmed.end(), intensity_temp.begin(), intensity_temp.end());
    }

    // 更新数组为剩余部分 - 使用更安全的方式
    if (idx_left == 0) {
      break; // 没有剩余元素了
    }

    Rcpp::NumericVector new_mass(idx_left);
    Rcpp::NumericVector new_intensity(idx_left);

    for (size_t i = 0; i < idx_left; ++i) {
      new_mass[i] = mass_sorted[i];
      new_intensity[i] = intensity_sorted[i];
    }

    mass_sorted = new_mass;
    intensity_sorted = new_intensity;
  }

  // 如果没有确认的峰，返回NA矩阵
  if (mass_confirmed.empty()) {
    Rcpp::NumericMatrix empty(1, 2);
    empty(0, 0) = NA_REAL;
    empty(0, 1) = NA_REAL;
    return empty;
  }

  // 创建结果矩阵
  Rcpp::NumericMatrix result(mass_confirmed.size(), 2);

  // 按m/z排序
  std::vector<size_t> order(mass_confirmed.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;

  std::sort(order.begin(), order.end(),
            [&mass_confirmed](size_t i1, size_t i2) {return mass_confirmed[i1] < mass_confirmed[i2];});

  // 填充结果矩阵
  for (size_t i = 0; i < order.size(); ++i) {
    result(i, 0) = mass_confirmed[order[i]];
    result(i, 1) = intensity_confirmed[order[i]];
  }

  // 设置列名
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("mz", "intensity");

  return result;
}
