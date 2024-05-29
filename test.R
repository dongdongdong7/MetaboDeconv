devtools::document()

# Load XcmsExperiment
load("./test_data/swath_data.RData")
load("./test_data/swath_spectra.RData")
# One sample
swath_data # ndata 第n个样本
swath_spectra

featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
                                       xcms::chromPeakData(swath_data)),
                                 rownames = "feature_id")
featureTable_ms1 <- featureTable %>%
  dplyr::filter(ms_level == 1)
featureTable_ms2 <- featureTable %>%
  dplyr::filter(ms_level == 2)

fenamiphos_mz <- 304.113077
feature_fenamiphos <- featureTable_ms1 %>%
  dplyr::filter(dplyr::near(mz, fenamiphos_mz, tol = 0.001))
fenamiphos_rt <- feature_fenamiphos$rt
peakWidth <- feature_fenamiphos$rtmax - feature_fenamiphos$rtmin
rtRange <- c(fenamiphos_rt - 2.5 * peakWidth, fenamiphos_rt + 2.5 * peakWidth)
mfeatureTable_ms2 <- featureTable_ms2 %>%
  dplyr::filter(rt >= rtRange[1] & rt <= rtRange[2]) %>%
  dplyr::filter(fenamiphos_mz >= isolationWindowLowerMz & fenamiphos_mz <= isolationWindowUpperMz)
spectra_ms2 <- xcms::spectra(swath_data) %>%
  Spectra::filterMsLevel(2L)

i <- 4
mfeatureTable_ms2_i <- mfeatureTable_ms2[i, ]
mfeatureSpectra_ms2_i <- spectra_ms2 %>%
  Spectra::filterRt(c(mfeatureTable_ms2_i$rtmin, mfeatureTable_ms2_i$rtmax)) %>%
  Spectra::filterMzRange(mz = c(mfeatureTable_ms2_i$mzmin, mfeatureTable_ms2_i$mzmax), keep = TRUE)
rtSpectra_i <- Spectra::rtime(mfeatureSpectra_ms2_i)
peaksData_i <- Spectra::peaksData(mfeatureSpectra_ms2_i)
chrDf_i <- purrr::list_rbind(lapply(1:length(rtSpectra_i), function(j) {
  rt_tmp <- rtSpectra_i[j]
  peakMat <- peaksData_i[[j]]
  if(nrow(peakMat) == 0){
    return(NULL)
  }else{
    peakDf <- as.data.frame(peakMat)
    peakDf$rt <- rt_tmp
    return(peakDf)
  }
}))
# Multi sample
load("./test_data/data_MSE.RData")
