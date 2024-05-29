#' @title get_chrDfList
#' @description
#' Get all related ms2 features chromatogram information of the mth ms1 feature.
#'
#' @param ndata nth sample data.
#' @param m mth feature
#'
#' @return A list contains chrDf which belongs to mth ms1 feature.
#' @export
#'
#' @examples
#' load("./test_data/swath_data.RData")
#' mchrDfList <- get_chrDfList(ndata = swath_data, m = 34)
get_chrDfList <- function(ndata, m){
  spectra_ms1 <- xcms::spectra(ndata) %>%
    Spectra::filterMsLevel(1L)
  spectra_ms2 <- xcms::spectra(ndata) %>%
    Spectra::filterMsLevel(2L)
  featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                         xcms::chromPeakData(ndata)),
                                   rownames = "feature_id")
  featureTable_ms1 <- featureTable %>%
    dplyr::filter(ms_level == 1)
  featureTable_ms2 <- featureTable %>%
    dplyr::filter(ms_level == 2)

  mfeatureTable_ms1 <- featureTable_ms1[m, ]
  peakWidth <- mfeatureTable_ms1$rtmax - mfeatureTable_ms1$rtmin
  rtRange <- c(mfeatureTable_ms1$rt - 2.5 * peakWidth,
               mfeatureTable_ms1$rt + 2.5 * peakWidth)
  mfeatureTable_ms2 <- featureTable_ms2 %>%
    dplyr::filter(rt >= rtRange[1] & rt <= rtRange[2]) %>%
    dplyr::filter(mfeatureTable_ms1$mz >= isolationWindowLowerMz & mfeatureTable_ms1$mz <= isolationWindowUpperMz)
  # remove ms2 feature mz >= precursorMz
  mfeatureTable_ms2 <- mfeatureTable_ms2 %>%
    dplyr::filter(mzmax < mfeatureTable_ms1$mzmin)
  # remove noise
  mfeatureTable_ms2 <- mfeatureTable_ms2 %>%
    dplyr::filter(maxo > max(mfeatureTable_ms2$maxo) * 0.05)
  if(nrow(mfeatureTable_ms2) == 0){
    chrDfList_ms2 <- NULL
  }else{
    loop <- function(i){
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
      return(chrDf_i)
    }
    chrDfList_ms2 <- lapply(1:nrow(mfeatureTable_ms2), function(i) {
      loop(i)
    })
    names(chrDfList_ms2) <- mfeatureTable_ms2$feature_id
  }
  mfeatureSpectra_ms1 <- spectra_ms1 %>%
    Spectra::filterRt(c(mfeatureTable_ms1$rtmin, mfeatureTable_ms1$rtmax)) %>%
    Spectra::filterMzRange(mz = c(mfeatureTable_ms1$mzmin, mfeatureTable_ms1$mzmax), keep = TRUE)
  rtSpectra <- Spectra::rtime(mfeatureSpectra_ms1)
  peaksData <- Spectra::peaksData(mfeatureSpectra_ms1)
  chrDf_ms1 <- purrr::list_rbind(lapply(1:length(rtSpectra), function(j) {
    rt_tmp <- rtSpectra[j]
    peakMat <- peaksData[[j]]
    if(nrow(peakMat) == 0){
      return(NULL)
    }else{
      peakDf <- as.data.frame(peakMat)
      peakDf$rt <- rt_tmp
      return(peakDf)
    }
  }))
  chrDfList_ms1 <- list(chrDf_ms1)
  names(chrDfList_ms1) <- mfeatureTable_ms1$feature_id
  chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
  return(chrDfList)
}
#' @title plot_chrDfList
#'
#' @param chrDfList
#'
#' @return
#' @export
#'
#' @examples
plot_chrDfList <- function(chrDfList){

}


