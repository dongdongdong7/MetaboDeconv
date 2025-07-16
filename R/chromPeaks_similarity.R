# Calculation of chromPeaks similarity
# Barry Song
# 20250620

ccf_peak2peak <- function(intensity1, rtime1, intensity2, rtime2){
  rtime_interval <- round(mean(diff(rtime1)), 1)
  rtime_c <- c(rtime1, rtime2)
  rtime_range <- round(c(min(rtime_c), max(rtime_c)), 1)
  rtime_approx <- seq(rtime_range[1], rtime_range[2], rtime_interval)
  intensity1_ <- approx(rtime1, intensity1, xout = rtime_approx, rule = 2)$y
  intensity2_ <- approx(rtime2, intensity2, xout = rtime_approx, rule = 2)$y
  ccfRes <- ccf(intensity1_, intensity2_, plot = FALSE, type = "correlation")
  imax <- which.max(ccfRes$acf)
  c(ccfRes$acf[imax], ccfRes$lag[imax])
}

# ms1Peak is a matrix with rtime and intensity
# ms2PeaksList is a list contains ms2Peak which is a matrix with rtime and intensity
ccf_ms1Peak_2_ms2Peaks <- function(ms1Peak, ms2PeaksList){
  # calculate rtime interval
  rtime_interval <- round(mean(diff(ms1Peak[, "rtime"])), 1)
  rtime_range <- round(c(min(ms1Peak[, "rtime"]), max(ms1Peak[, "rtime"])), 1)
  rtime_approx <- seq(rtime_range[1], rtime_range[2], rtime_interval)
  ms1Intensity <- approx(x = ms1Peak[, "rtime"], y = ms1Peak[, "intensity"], xout = rtime_approx, rule = 2)$y
  ms2IntensityList <- lapply(ms2PeaksList, function(ms2Peak) {
    approx(x = ms2Peak[, "rtime"], y = ms2Peak[, "intensity"], xout = rtime_approx, rule = 2)$y
  })
  res_ <- mapply(function(ms1Intensity_, ms2Intensity_) {
    ccfRes <- ccf(ms1Intensity_, ms2Intensity_, plot = FALSE)
    imax <- which.max(ccfRes$acf)
    if(length(imax) == 0) return(c(0, 0))
    c(ccfRes$acf[imax], ccfRes$lag[imax])
  }, ms2Intensity_ = ms2IntensityList,
  MoreArgs = list(ms1Intensity_ = ms1Intensity),
  SIMPLIFY = TRUE)
  res_ <- t(res_)
  colnames(res_) <- c("score", "lag")
  return(res_)
}


#' @title ccf_ms1PeakAms2Peaks
#' @description
#' Calculate ccf score matrix for ms1Peak and its corresponding ms2Peaks
#'
#' @param ms1PeakDT `data.table()`, a one row data.table with ms1 peak information
#'  (mz, mzmin, mzmax, rt, rtmin, rtmax, maxo, intensity, rtime)
#' @param ms2PeaksDT `data.table()`, a data.table with ms2 peaks information corresponding ms1 peak
#'
#' @returns `matrix()` with score, lag and maxo
#' @export
#'
#' @examples
#' ccfMatrix <- ccf_ms1PeakAms2Peaks(ms1PeakDT = ms1PeakDT, ms2PeaksDT = ms2PeaksDT)
ccf_ms1PeakAms2Peaks <- function(ms1PeakDT, ms2PeaksDT){
  # calculate rtime interval
  rtime_interval <- round(mean(diff(ms1PeakDT[1, rtime][[1]])), 1)
  rtime_range <- round(c(min(ms1PeakDT[1, rtime][[1]]), max(ms1PeakDT[1, rtime][[1]])), 1)
  # Rentention Time Vector
  rtime_approx <- seq(rtime_range[1], rtime_range[2], rtime_interval)

  # Intensity vector
  ms1Intensity <- approx(x = ms1PeakDT[1, rtime][[1]], y = ms1PeakDT[1, intensity][[1]],
                         xout = rtime_approx, rule = 2)$y
  ms2IntensityList <- mapply(function(intensity, rtime, xout){
    approx(x = rtime, y = intensity, xout = xout, rule = 2)$y
  }, intensity = ms2PeaksDT[, intensity], rtime = ms2PeaksDT[, rtime], MoreArgs = list(xout = rtime_approx), SIMPLIFY = FALSE)

  # ccf matrix
  res_ <- mapply(function(ms1Intensity_, ms2Intensity_, maxo_) {
    ccfRes <- ccf(ms1Intensity_, ms2Intensity_, plot = FALSE)
    imax <- which.max(ccfRes$acf)
    if(length(imax) == 0) return(c(0, 0, 0))
    c(score = round(ccfRes$acf[imax], 2), lag = ccfRes$lag[imax], maxo = maxo_)
  }, ms2Intensity_ = ms2IntensityList, maxo_ = ms2PeaksDT[, maxo], MoreArgs = list(ms1Intensity_ = ms1Intensity), SIMPLIFY = TRUE)
  t(res_)
}
