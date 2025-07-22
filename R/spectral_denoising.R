# Spectral Denoising
# Barry Song
# 250722

.simulate_msms <- function() {
  # 真实信号参数
  n_peaks <- 10
  mass_real <- sort(runif(n_peaks, 100, 1000))
  intensity_real <- abs(rnorm(n_peaks, mean = 1000, sd = 300))

  # 电子噪音参数
  n_noise <- 50
  mass_noise <- runif(n_noise, 100, 1000)
  intensity_noise <- abs(rnorm(n_noise, mean = 50, sd = 0.1))

  # 合并信号和噪音
  mass <- c(mass_real, mass_noise)
  intensity <- c(intensity_real, intensity_noise)

  # 随机排序
  set.seed(123)
  idx <- sample(seq_along(mass))
  msms <- cbind(mass[idx], intensity[idx])

  return(msms)
}

# msms: `matrix()` with mz and intensity
#' @rdname spectral_denoising
#' @title Spectral Denoising
#' @description
#' Noise reduction for MS/MS, including electronic noise and chemical noise.
#' Electronic denoising method refers from reference (Denoising Search doubles the number of metabolite and exposome annotations in human plasma using an Orbitrap Astral mass spectrometer)
#' Chemical denoising method is similar with MS2Purifier
#' @param msms `matrix()` with mz and intensity
#'
#' @returns `matrix()` with mz and intensity
#' @export
#'
#' @examples
#' set.seed(110)
#' msms_data <- .simulate_msms()
#' result_matrix <- electronic_denoising(msms_data)
#'
#' plot(msms_data[,1], msms_data[,2], type = "h",
#'      main = "原始质谱数据 (含电子噪音)",
#'      xlab = "m/z", ylab = "Intensity")
#' plot(result_matrix[, "mz"], result_matrix[, "intensity"], type = "h",
#'      main = "去噪后的质谱数据 (矩阵格式)",
#'      xlab = "m/z", ylab = "Intensity")
electronic_denoising <- function(msms){
  electronic_denoising_cpp(msms = msms)
}

