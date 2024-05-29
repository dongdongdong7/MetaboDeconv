library(xcms)
library(MsExperiment)
swath_file <- system.file("TripleTOF-SWATH",
                          "PestMix1_SWATH.mzML",
                          package = "msdata")

swath_data <- readMsExperiment(swath_file)
swath_data <- filterRt(swath_data, rt = c(230, 610))
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaks(swath_data, param = cwp)
cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp)
swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.9)
save(swath_data, file = "./test_data/swath_data.RData")
save(swath_spectra, file = "./test_data/swath_spectra.RData")
