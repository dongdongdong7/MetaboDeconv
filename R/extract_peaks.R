# Extract ms1 peaks from target mz or extract ms2 peaks from ms1PeaksDT
# Barry Song
# 250620

#' @rdname extract_peaks
#' @title Extract ms1 or ms2 peaks
#' @description
#' Extract ms1 peaks from target mz or extract ms2 peaks from ms1PeaksDT
#'
#' @param targetMz target m/z
#' @param dps_ms1 `data.table()` with mz, intensity, scanIndex
#' @param scanIndex_ms1 scanIndex of ms1 spectra
#' @param rtime_ms1 retention time of ms1 spectra
#' @param ppm ppm
#' @param peakwidth peakwidth
#' @param snthresh snthresh
#' @param noise noise
#' @param r2thresh r2thresh
#' @param peakWidthExtend peakWidthExtend
#' @param fill_missing whether to fill in missing values in chromatogram
#'
#' @returns ms1PeaksDT `data.table()` with mz, mzmin, mzmax, rt, rtmin, rtmax, into, intb, maxo, r2, intensity(list), rtime(list)
#' @export
#' @examples
#' ms1PeaksDT <- extract_ms1Peaks()
extract_ms1Peaks <- function(targetMz,
                             dps_ms1, scanIndex_ms1, rtime_ms1,
                             ppm = 5,
                             peakwidth = c(5, 20), snthresh = 3, noise = 100, r2thresh = 0.8,
                             peakWidthExtend = 1,
                             fill_missing = FALSE){
  targetMzRange <- lapply(targetMz, function(x) {
    tol_mz <- MsCoreUtils::ppm(x = x, ppm = ppm)
    c(x - tol_mz, x + tol_mz)
  })
  pb <- progress::progress_bar$new(
    format = "[:bar] :elapsed | progress: :percent",
    total = length(targetMzRange),
    clear = FALSE,
    width = 60
  )
  ms1PeaksDTList <- lapply(targetMzRange, function(mzrange_) {
    pb$tick()
    dps_ <- dps_ms1[data.table::between(x = mz, lower = mzrange_[1], upper = mzrange_[2]), ]
    if(nrow(dps_) < 3) return(NULL)
    dps_ <- dps_[, .(
      mz = mean(mz, na.rm = TRUE),
      intensity = sum(intensity, na.rm = TRUE)
    ), by = scanIndex]
    intensity_ <- numeric(length(scanIndex_ms1))
    consensus_mz <- mean(c(dps_[which.max(dps_[, intensity]), mz], median(dps_[, mz])), na.rm = TRUE)
    intensity_[match(dps_[, scanIndex], scanIndex_ms1)] <- dps_[, intensity]
    if(fill_missing) intensity_ <- fill_missing_signals(intensity = intensity_, gap = 1, noise = noise)
    ps_ <- findChromPeaks::findChromPeaks_CWT(int = intensity_, rt = rtime_ms1,
                                              peakwidth = peakwidth, snthresh = snthresh, noise = noise,
                                              r2thresh = r2thresh)
    if(nrow(ps_) == 0) return(NULL)
    data.table::rbindlist(lapply(1:nrow(ps_), function(i) {
      peakWidth_p <- ps_[i, "rtmax"] - ps_[i, "rtmin"]
      rtmin_p <- ps_[i, "rtmin"]
      rtmax_p <- ps_[i, "rtmax"]
      rtmin_p_extend <- max(rtime_ms1[1], rtmin_p - peakWidth_p * peakWidthExtend)
      rtmax_p_extend <- min(tail(rtime_ms1, 1), rtmax_p + peakWidth_p * peakWidthExtend)
      l_p <- rtime_ms1 >= rtmin_p & rtime_ms1 <= rtmax_p
      l_p_extend <- rtime_ms1 >= rtmin_p_extend & rtime_ms1 <= rtmax_p_extend

      intensity_ini <- numeric(length(scanIndex_ms1))
      intensity_ini[l_p] <- intensity_[l_p]
      intensity_p_extend <- intensity_ini[l_p_extend]
      rtime_p_extend <- rtime_ms1[l_p_extend]
      data.table::data.table(
        mz = consensus_mz,
        mzmin = mzrange_[1],
        mzmax = mzrange_[2],
        rt = ps_[i, "rt"],
        rtmin = ps_[i, "rtmin"],
        rtmax = ps_[i, "rtmax"],
        into = ps_[i, "into"],
        intb = ps_[i, "intb"],
        maxo = ps_[i, "maxo"],
        r2 = ps_[i, "r2"],
        intensity = list(intensity_p_extend),
        rtime = list(rtime_p_extend)
      )
    }))
  })
  ms1PeakDT <- data.table::rbindlist(ms1PeaksDTList)
  ms1PeakDT
}

#' @rdname extract_peaks
#' @param ms1PeakDT `data.table()`, one row of ms1PeaksDT
#' @param fragmentMz `numeric()`, fragment m/z of ms1 peak
#' @param dps_ms2 dps_ms2 `data.table()` with mz, intensity and scanIndex. Data points of ms2 spectra.
#' @param scanIndex_ms2 scanIndex_ms2 scanIndex of ms2 spectra
#' @param rtime_ms2 rtime_ms2 retention time of ms2 spectra
#' @return ms2PeaksDT `data.table()` with mz, mzmin, mzmax, rt, rtmin, rtmax, maxo, intensity(list), rtime(list)
#' @export
extract_ms2Peaks <- function(ms1PeakDT, fragmentMz,
                             dps_ms2, scanIndex_ms2, rtime_ms2,
                             ppm = 5,
                             peakwidth = c(5, 20), snthresh = 3, noise = 100, r2thresh = 0.8,
                             fill_missing = FALSE){
  rtime_p_ini <- ms1PeakDT[1, rtime][[1]]
  l_p <- rtime_ms2 >= rtime_p_ini[1] & rtime_ms2 <= tail(rtime_p_ini, 1)
  scanIndex_p <- scanIndex_ms2[l_p]
  rtime_p <- rtime_ms2[l_p]
  dps_ms2 <- dps_ms2[scanIndex %in% scanIndex_p, ]

  precursorMzRange <- c(ms1PeakDT[1, mzmin], ms1PeakDT[1, mzmax])
  fragmentMzRangeList <- lapply(fragmentMz, function(mz_) {
    tol_mz <- MsCoreUtils::ppm(mz_, ppm = ppm)
    c(mz_ - tol_mz, mz_ + tol_mz)
  })
  ms2PeaksDTList <- lapply(fragmentMzRangeList, function(mzrange_) {
    dps_ <- dps_ms2[data.table::between(mz, mzrange_[1], mzrange_[2])]
    if(nrow(dps_) < 3) return(NULL)
    dps_ <- dps_[, .(
      mz = mean(mz, na.rm = TRUE),
      intensity = sum(intensity, na.rm = TRUE)
    ), by = scanIndex]
    intensity_ <- numeric(length(scanIndex_p))
    max_intensity_mz <- dps_[which.max(intensity), mz]
    median_mz <- median(dps_$mz, na.rm = TRUE)
    consensus_mz <- mean(c(max_intensity_mz, median_mz), na.rm = TRUE)
    intensity_[match(dps_$scanIndex, scanIndex_p)] <- dps_$intensity
    if(fill_missing) intensity_ <- fill_missing_signals(intensity = intensity_, gap = 1, noise = noise)
    ps_ <- findChromPeaks::findChromPeaks_CWT(int = intensity_, rt = rtime_p,
                                              snthresh = snthresh, peakwidth = peakwidth, noise = noise,
                                              r2thresh = r2thresh)
    if(nrow(ps_) == 0) return(NULL)
    data.table::rbindlist(lapply(1:nrow(ps_), function(i) {
      l_p <- rtime_p >= ps_[i, "rtmin"] & rtime_p <= ps_[i, "rtmax"]

      intensity_new <- numeric(length(scanIndex_p))
      intensity_new[l_p] <- intensity_[l_p]
      data.table::data.table(
        mz = consensus_mz,
        mzmin = mzrange_[1],
        mzmax = mzrange_[2],
        rt = ps_[i, "rt"],
        rtmin = ps_[i, "rtmin"],
        rtmax = ps_[i, "rtmax"],
        maxo = ps_[i, "maxo"],
        r2 = ps_[i, "r2"],
        intensity = list(intensity_new),
        rtime = list(rtime_p)
      )
    }))
  })
  data.table::rbindlist(ms2PeaksDTList)
}
