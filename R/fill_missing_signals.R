# intensity: intensity vector
# gap: number of missing signals allowed in a chromatographic peak
# noise: noise value of chromatogram
.detect_missing <- function(intensity, gap = 1, noise = 100){
  .detect_missing_cpp(intensity = intensity, gap = gap, noise = noise)
}

#' @title Fill missing signals for chromatographic peak
#' @description
#' Use approx function to fill missing signals for chromatographic peak
#'
#' @param intensity `numeric()`, intensity vector
#' @param gap `integer(1)`, number of missing signals allowed in a chromatographic peak
#' @param noise `numeric(1)`, noise value of chromatogram
#'
#' @returns `numeric()`, filled intensity vector
#' @export
fill_missing_signals <- function(intensity, gap = 1, noise = 100){
  missing_idx <- .detect_missing(intensity = intensity, gap = gap, noise = noise)
  if(length(missing_idx) == 0) return(intensity)
  x <- 1:length(intensity)
  intensity[missing_idx] <- approx(x[-missing_idx], intensity[-missing_idx], xout = missing_idx, rule = 2)$y
  return(intensity)
}
