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
      chrDf_i$smooth_int <- smoothMean(chrDf_i$intensity, size = 3)
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
  chrDf_ms1$smooth_int <- smoothMean(chrDf_ms1$intensity, size = 3)
  chrDfList_ms1 <- list(chrDf_ms1)
  names(chrDfList_ms1) <- mfeatureTable_ms1$feature_id
  chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
  return(chrDfList)
}
#' @title smoothMean
#' @description
#' Smooth function using mean.
#'
#' @param int intensity.
#' @param size mean size.
#'
#' @return intensity vector
#' @export
#'
#' @examples
#' intensity <- c(1701.196, 2121.159, 2401.849, 2043.653, 1542.964,  723.803)
#' smoothMean(int = intensity)
smoothMean <- function(int, size = 3){
  if(size %% 2 == 0){
    stop("size should be a singular!")
  }
  half_size <- floor(size / 2)
  N <- length(int)
  special_idx <- c(seq_len(half_size), rev((N:(N - half_size + 1))))
  smooth_int <- sapply(1:N, function(i) {
    if(i %in% special_idx){
      return(int[i])
    }
    y <- sapply((i - half_size):(i + half_size), function(j){
      int[j]
    })
    smooth_int <- mean(y)
  })
  return(smooth_int)
}
#' @title plot_chrDfList
#' @description
#' plot_chrDfList
#'
#' @param chrDfList chrDfList
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_chrDfList(chrDfList = mchrDfList)
plot_chrDfList <- function(chrDfList){
  chrDfList_ms1 <- chrDfList$ms1
  feature_ms1_name <- names(chrDfList_ms1)
  chrDf_ms1 <- chrDfList_ms1[[1]]
  chrDfList_ms2 <- chrDfList$ms2
  dfList <- lapply(1:length(chrDfList_ms2), function(i) {
    df_tmp <- chrDfList_ms2[[i]]
    df_tmp$level <- "2"
    df_tmp$feature <- i
    return(df_tmp)
  })
  df_tmp <- chrDf_ms1
  df_tmp$level <- "1"
  df_tmp$feature <- length(dfList) + 1
  dfList <- append(dfList, list(df_tmp))
  df <- purrr::list_rbind(dfList)
  group <- unique(df$level)
  group_number <- length(group)
  group_color <- RColorBrewer::brewer.pal(6, "Set1")[1:group_number]
  names(group_color) <- group
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = rt, group = feature, color = level)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = smooth_int), linewidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention Time", y = "Intensity") +
    ggplot2::scale_color_manual(values = group_color)
  return(p)
}


