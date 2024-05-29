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
mchrDfList <- get_chrDfList(ndata = swath_data, m = 34)
plot_chrDfList(chrDfList = mchrDfList)



# Multi sample
load("./test_data/data_MSE.RData")
