devtools::document()

# Load XcmsExperiment
load("./test_data/swath_data.RData")
load("./test_data/swath_spectra.RData")
fenamiphos <- Spectra::Spectra(
  system.file("mgf", "metlin-72445.mgf", package = "xcms"),
  source = MsBackendMgf::MsBackendMgf())
fenamiphos_spMat <- sp2spMat(fenamiphos[2])
MetaboSpectra::plotSpectra(fenamiphos_spMat)
# One sample
swath_data # ndata 第n个样本
swath_spectra
chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
                                        xcms::chromPeakData(swath_data)),
                                  rownames = "cpid")
chromPeakTable_ms1 <- chromPeakTable %>%
 dplyr::filter(ms_level == 1)
chromPeakTable_ms2 <- chromPeakTable %>%
 dplyr::filter(ms_level == 2)
# filter chromPeakTable
chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
  dplyr::filter(maxo >= 1000)
chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
  dplyr::filter(maxo > 100)
chromPeakTable_ms1 <- Deconv4ndata(ndata = swath_data, chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2,
                                   thread = 4, st = 0.7)
Spectra::plotSpectra(chromPeakTable_ms1[9, ]$spectra[[1]])

featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
                                       xcms::chromPeakData(swath_data)),
                                 rownames = "feature_id")
featureTable_ms1 <- featureTable %>%
  dplyr::filter(ms_level == 1)
featureTable_ms2 <- featureTable %>%
  dplyr::filter(ms_level == 2)

fenamiphos_mz <- 304.113077
fenamiphos_ms1_peak <- xcms::chromPeaks(swath_data, mz = fenamiphos_mz, ppm = 2)
mchrDfList <- get_chrDfList(ndata = swath_data, m = 34, smooth = TRUE, noise = 10)
plot_chrDfList(chrDfList = mchrDfList)
mchrDfList_new <- filterChrDf(chrDfList = mchrDfList, weight_rt = 0.5, weight_shape = 0.5,st = 0.70)
plot_chrDfList(chrDfList = mchrDfList_new)
sp_ms2 <- chrDfList2spectra(chrDfList = mchrDfList_new, ndata = swath_data)
sp_ms2_spMat <- sp2spMat(sp_ms2)
sp_ms2_spMat <- MetaboSpectra::clean_spMat(sp_ms2_spMat)
MetaboSpectra::plotSpectra(sp_ms2_spMat)
fenamiphos_swath_spectrum <- swath_spectra[
  swath_spectra$peak_id == rownames(fenamiphos_ms1_peak)]
fenamiphos_swath_spectrum_spMat <- sp2spMat(fenamiphos_swath_spectrum)
fenamiphos_swath_spectrum_spMat <- MetaboSpectra::clean_spMat(fenamiphos_swath_spectrum_spMat)
MetaboSpectra::plotComparableSpectra(sp_ms2_spMat, fenamiphos_spMat)
MetaboSpectra::compare_spMat_entropy(sp_ms2_spMat,fenamiphos_spMat)
MetaboSpectra::compare_spMat_ndotproduct(sp_ms2_spMat,fenamiphos_spMat)
MetaboSpectra::plotComparableSpectra(fenamiphos_swath_spectrum_spMat, fenamiphos_spMat)
MetaboSpectra::compare_spMat_entropy(fenamiphos_swath_spectrum_spMat,fenamiphos_spMat)
MetaboSpectra::compare_spMat_ndotproduct(fenamiphos_swath_spectrum_spMat,fenamiphos_spMat)
# Multi sample
load("./test_data/QC/data_MSE.RData")
nLength <- length(data_MSE)
featureTable_ms1List <- lapply(1:nLength, function(n) {
  ndata <- data_MSE[n]
  chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                           xcms::chromPeakData(ndata)),
                                     rownames = "cpid")
  chromPeakTable_ms1 <- chromPeakTable %>%
    dplyr::filter(ms_level == 1)
  chromPeakTable_ms2 <- chromPeakTable %>%
    dplyr::filter(ms_level == 2)
  # filter chromPeakTable
  chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
    dplyr::filter(maxo >= 10000)
  chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
    dplyr::filter(maxo > 1000)
  plot_chrDf(generate_chrDf(ndata, "CP00347", noise = 10, smooth = TRUE, size = 3))
  mchrDfList <- generate_chrDfList(ndata = ndata, cpid = "CP00347",
                                  chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2)
  plot_chrDfList(mchrDfList)
  mchrDfList_new <- filterChrDf(chrDfList = mchrDfList, weight_rt = 0.5, weight_shape = 0.5,st = 0.70)
  plot_chrDfList(chrDfList = mchrDfList_new)
  chromPeakTable_ms1 <- Deconv4ndata(ndata = ndata, chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2,
                                     thread = 1, st = 0.7)
  Spectra::plotSpectra(chromPeakTable_ms1[9, ]$spectra[[1]])
})
