# Extract ms2Peaks from a ms1Peak
# Barry Song
# 250617

# ms1PeakDT: mz, mzmin, mzmax, rt, rtmin, rtmax, intensity, rtime
extract_ms2Peaks <- function(ms1PeakDT, fragmentMz,
                             dps_ms2, scanIndex_ms2, rtime_ms2,
                             peakwidth = c(5, 20), ppm = 5, snthresh = 3){
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
    ps_ <- findChromPeaks::findChromPeaks_CWT(int = intensity_, rt = rtime_p,
                                              snthresh = snthresh, peakwidth = peakwidth)
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
        intensity = list(intensity_new),
        rtime = list(rtime_p)
      )
    }))
  })
  data.table::rbindlist(ms2PeaksDTList)
}
