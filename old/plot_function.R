#' @title plot_chrDf
#' @description
#' Plot a chromatogram using chromatographic data frame.
#'
#' @param chrDf chrDf.
#' @param linewidth linewidth.
#'
#' @return ggplot object.
#' @export
#'
#' @examples
#' plot_chrDf(chrDf = chrDf_i, linewidth = 2)
plot_chrDf <- function(chrDf, linewidth = 1){
  title <- paste0(round(min(chrDf$mz), 4), " - ", round(max(chrDf$mz), 4))
  p <- ggplot2::ggplot(data = chrDf, mapping = ggplot2::aes(x = rt)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = intensity), linewidth = linewidth) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention Time", y = "Intensity", title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 15),
                  axis.title = ggplot2::element_text(size = 15),
                  axis.text = ggplot2::element_text(size = 10))
  return(p)
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
plot_chrDfList <- function(chrDfList){
  chrDfList_ms1 <- chrDfList$ms1
  feature_ms1_name <- names(chrDfList_ms1)
  chrDf_ms1 <- chrDfList_ms1[[1]]
  chrDfList_ms2 <- chrDfList$ms2
  if(is.null(chrDfList_ms2)) dfList <- list()
  else{
    dfList <- lapply(1:length(chrDfList_ms2), function(i) {
      df_tmp <- chrDfList_ms2[[i]]
      df_tmp$level <- "2"
      df_tmp$feature <- i
      return(df_tmp)
    })
  }
  df_tmp <- chrDf_ms1
  df_tmp$level <- "1"
  df_tmp$feature <- length(dfList) + 1
  dfList <- append(dfList, list(df_tmp))
  df <- purrr::list_rbind(dfList)
  df <- df %>%
    dplyr::arrange(level)
  group <- unique(df$level)
  group_number <- length(group)
  group_color <- RColorBrewer::brewer.pal(6, "Set1")[1:group_number]
  names(group_color) <- group
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = rt, group = feature, color = level)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = intensity), linewidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention Time", y = "Intensity") +
    ggplot2::scale_color_manual(values = group_color)
  return(p)
}
