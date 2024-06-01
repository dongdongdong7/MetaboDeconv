#' @title generate_chrDf
#' @description
#' Generate a chrDf.
#'
#' @param ndata nth sample data.
#' @param cpid chromatographic peak id.
#' @param noise noise value.
#' @param smooth whether or not smooth.
#' @param size smoothMean size.
#'
#' @return chrDf.
#' @export
#'
#' @examples
#' load("./test_data/swath_data.RData")
#' chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
#'                                          xcms::chromPeakData(swath_data)),
#'                                    rownames = "cpid")
#' chromPeakTable_ms1 <- chromPeakTable %>%
#'   dplyr::filter(ms_level == 1)
#' chromPeakTable_ms2 <- chromPeakTable %>%
#'   dplyr::filter(ms_level == 2)
#' chrDf_i <- generate_chrDf(ndata = swath_data, cpid = "CP34", smooth = TRUE, size = 3)
generate_chrDf <- function(ndata, cpid, noise = 10, smooth = TRUE, size = 3){
  chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                           xcms::chromPeakData(ndata)),
                                     rownames = "cpid")
  chromPeak <- chromPeakTable[chromPeakTable$cpid == cpid, ]
  ms_level <- chromPeak$ms_level
  spectra <- xcms::spectra(ndata) %>%
    Spectra::filterMsLevel(ms_level)
  cpSpectra <- spectra %>%
    Spectra::filterRt(c(chromPeak$rtmin, chromPeak$rtmax)) %>%
    Spectra::filterMzRange(mz = c(chromPeak$mzmin, chromPeak$mzmax))
  rtSpectra <- Spectra::rtime(cpSpectra)
  peaksData <- Spectra::peaksData(cpSpectra)
  chrDf <- purrr::list_rbind(lapply(1:length(cpSpectra), function(j) {
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
  chrDf <- chrDf %>%
    dplyr::filter(intensity > noise)
  if(smooth) chrDf$intensity <- smoothMean(chrDf$intensity, size = size)
  return(chrDf)
}
#' @title generate_chrDfList
#' @description
#' Generate a chrDfList, contains ms1 chrDf and associated ms2 chrDf. And chromatographic peak will be
#' detect from chromPeakTable_ms1 and chromPeakTable_ms2.
#' @param ndata nth data.
#' @param cpid chromatographic peak id.
#' @param chromPeakTable_ms1 1 level chromPeakTable which have been filtered.(must be generate by ndata)
#' @param chromPeakTable_ms2 2 level chromPeakTable which have been filtered.(must be generate by ndata)
#' @param noise noise value.
#' @param smooth whether or not smooth.
#' @param size smoothMean size.
#' @param factor factor * peakwidth -> rtRange.
#'
#' @return A chrDf.
#' @export
#'
#' @examples
#' load("./test_data/swath_data.RData")
#' chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
#'                                         xcms::chromPeakData(swath_data)),
#'                                   rownames = "cpid")
#' chromPeakTable_ms1 <- chromPeakTable %>%
#'  dplyr::filter(ms_level == 1)
#' chromPeakTable_ms2 <- chromPeakTable %>%
#'  dplyr::filter(ms_level == 2)
#' # filter chromPeakTable
#' chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
#'   dplyr::filter(maxo >= 1000)
#' chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
#'   dplyr::filter(maxo > 100)
#' chrDfList <- generate_chrDfList(ndata = swath_data, cpid = "CP34",
#'                                 chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2)
#' plot_chrDfList(chrDfList)
generate_chrDfList <- function(ndata, cpid, chromPeakTable_ms1,chromPeakTable_ms2,
                               noise = 10, smooth = TRUE, size = 3, factor = 0.5){
  if(!cpid %in% chromPeakTable_ms1$cpid) stop(paste0(cpid, " is not in chromPeakTable_ms1!"))
  mchromPeakTable_ms1 <- chromPeakTable_ms1[chromPeakTable_ms1$cpid == cpid, ]
  peakWidth <- mchromPeakTable_ms1$rtmax - mchromPeakTable_ms1$rtmin
  rtRange <- c(mchromPeakTable_ms1$rt - factor * peakWidth,
               mchromPeakTable_ms1$rt + factor * peakWidth)
  mchromPeakTable_ms2 <- chromPeakTable_ms2 %>%
    dplyr::filter(rt >= rtRange[1] & rt <= rtRange[2]) %>%
    dplyr::filter(mchromPeakTable_ms1$mz >= isolationWindowLowerMz & mchromPeakTable_ms1$mz <= isolationWindowUpperMz)
  # remove ms2 feature mz >= precursorMz
  mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
    dplyr::filter(mzmax < mchromPeakTable_ms1$mzmin)
  # remove noise
  mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
    dplyr::filter(maxo > max(mchromPeakTable_ms2$maxo) * 0.05)
  if(nrow(mchromPeakTable_ms2) == 0) return(NULL)
  else{
    chrDfList_ms2 <- lapply(mchromPeakTable_ms2$cpid, function(x) {
      chrDf_x <- generate_chrDf(ndata = ndata, cpid = x, noise = noise, smooth = smooth, size = size)
    })
    names(chrDfList_ms2) <- mchromPeakTable_ms2$cpid
  }
  chrDfList_ms1 <- list(generate_chrDf(ndata = ndata, cpid = mchromPeakTable_ms1$cpid,
                                  noise = noise, smooth = smooth, size = size))
  names(chrDfList_ms1) <- mchromPeakTable_ms1$cpid
  chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
  return(chrDfList)
}
#' @title filterChrDf
#' @description
#' Filter less relevant ms2 features
#'
#' @param chrDfList A chrDfList contains ms1 and ms2 chrDf
#'
#' @return A new chrDfList.
#' @export
#'
#' @examples
#' load("./test_data/swath_data.RData")
#' chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
#'                                         xcms::chromPeakData(swath_data)),
#'                                   rownames = "cpid")
#' chromPeakTable_ms1 <- chromPeakTable %>%
#'  dplyr::filter(ms_level == 1)
#' chromPeakTable_ms2 <- chromPeakTable %>%
#'  dplyr::filter(ms_level == 2)
#' # filter chromPeakTable
#' chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
#'   dplyr::filter(maxo >= 1000)
#' chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
#'   dplyr::filter(maxo > 100)
#' mchrDfList <- generate_chrDfList(ndata = swath_data, cpid = "CP34",
#'                                 chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2)
#' plot_chrDfList(mchrDfList)
#' mchrDfList_new <- filterChrDf(chrDfList = mchrDfList, weight_rt = 0.5, weight_shape = 0.5,st = 0.70)
#' plot_chrDfList(chrDfList = mchrDfList_new)
filterChrDf <- function(chrDfList, weight_rt = 0.5, weight_shape = 0.5, st = 0.8){
  chrDfList_ms1 <- chrDfList$ms1
  chrDf_ms1 <- chrDfList_ms1[[1]]
  chrDfList_ms2 <- chrDfList$ms2
  score <- sapply(1:length(chrDfList_ms2), function(i) {
    cor_rt <- calCor_rt(chrDf1 = chrDf_ms1, chrDf2 = chrDfList_ms2[[i]])
    cor_shape <- calCor_shape(chrDf1 = chrDf_ms1, chrDf2 = chrDfList_ms2[[i]])
    plot_chrDf(chrDfList_ms2[[i]])
    cor_both <- weight_rt * cor_rt + weight_shape * cor_shape
  })
  chrDfList_ms2 <- chrDfList_ms2[which(score > st)]
  return(list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2))
}
