# Chromatographic peaks deconvolution example using MetaboDeconv
# For the algorithm principle, please refer to DecoMetDIA
# Barry Song
# 250716

library(magrittr)
library(data.table)

ppm <- 30
snthresh <- 3
peakwidth <- c(5, 20)
r2thresh <- 0.8
noise <- 100
peakWidthExtend <- 1
fill_missing <- TRUE
scoreTh <- 0.95
lagTh <- 1
rt_diff_tol <- 10

# ms1PeaksDT: mz mzmin mzmax intensity rtime
# ms1PeakDT (one row of ms1PeaksDT): mz mzmin mzmax intensity rtime
# ms2PeaksDT (ms1PeakDT - ms2PeaksDT): mz mzmin mzmax intensity rtime

demo_path <- "../demo_data/pos/plate1_series05.mzML"

sps_dia <- Spectra::Spectra(object = demo_path,
                            source = Spectra::MsBackendMzR())
sps_dia <- sps_dia %>%
  Spectra::filterEmptySpectra() %>%
  Spectra::filterIntensity(noise)

sps_dia_ms1 <- sps_dia %>%
  Spectra::filterMsLevel(1L)
scanIndex_dia_ms1 <- Spectra::scanIndex(sps_dia_ms1)
rtime_dia_ms1 <- Spectra::rtime(sps_dia_ms1)
sps_dia_ms2 <- sps_dia %>%
  Spectra::filterMsLevel(2L)

# data points from ms1
dps_dia_ms1 <- MetaboSpectra::sps2dps(sps = sps_dia_ms1, ppm = ppm, noise = noise)

# Load standard_car
data("standard_car", package = "MetaboRI")
# L-Carnitine (HMDB0000062)
mz_adduct <- standard_car[1, ]$precursorMz

# Extract ms1 peak
ms1PeaksDT <- extract_ms1Peaks(targetMz = mz_adduct,
                               dps_ms1 = dps_dia_ms1, scanIndex_ms1 = scanIndex_dia_ms1, rtime_ms1 = rtime_dia_ms1,
                               ppm = ppm, peakwidth = peakwidth, snthresh = snthresh, noise = noise, r2thresh = r2thresh,
                               peakWidthExtend = peakWidthExtend,
                               fill_missing = fill_missing)
plot(ms1PeaksDT[1, rtime][[1]], ms1PeaksDT[1, intensity][[1]], type = "o")
# Extract ms2 peak
i <- 1 # ms1PeakDT <- ms1PeaksDT[i, ]
sps_i <- sps_dia_ms2 %>%
  Spectra::filterIsolationWindow(ms1PeaksDT[i, mz]) %>%
  Spectra::filterRt(c(ms1PeaksDT[i, rtime][[1]][1], tail(ms1PeaksDT[i, rtime][[1]], 1)))
dps_i <- MetaboSpectra::sps2dps(sps = sps_i, ppm = ppm, noise = noise)
scanIndex_i <- Spectra::scanIndex(sps_i)
rtime_i <- Spectra::rtime(sps_i)
scanIndex_i_apex <- scanIndex_i[which.min(abs(rtime_i - ms1PeaksDT[i, rt]))]
fragmentMz_i <- dps_i[scanIndex == scanIndex_i_apex, mz]
fragmentMz_i <- fragmentMz_i[fragmentMz_i < mz_adduct + MsCoreUtils::ppm(mz_adduct, ppm = ppm)] # Remove fragmentMz > precursorMz
ms2PeaksDT_i <- MetaboDeconv::extract_ms2Peaks(ms1PeakDT = ms1PeaksDT[i, ], fragmentMz = fragmentMz_i,
                                               dps_ms2 = dps_i, scanIndex_ms2 = scanIndex_i, rtime_ms2 = rtime_i,
                                               ppm = ppm, peakwidth = peakwidth, noise = noise,
                                               snthresh = snthresh, r2thresh = r2thresh, fill_missing = fill_missing)
ccfMatrix_i <- MetaboDeconv::ccf_ms1PeakAms2Peaks(ms1PeakDT = ms1PeaksDT[i, ], ms2PeaksDT = ms2PeaksDT_i)
ms2DT_i <- ms2PeaksDT_i[which(ccfMatrix_i[, "score"] >= scoreTh & abs(ccfMatrix_i[, "lag"]) <= lagTh), ]
# Only retain ms2 peaks shifted to the right:
# ms2DT_i <- ms2PeaksDT_i[which(ccfMatrix_i[, "score"] >= scoreTh & abs(ccfMatrix_i[, "lag"]) <= lagTh & ccfMatrix_i[, "lag"] >= 0), ]
spMat_ms2 <- matrix(c(ms2DT_i[, mz], ms2DT_i[, maxo]), ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
spMat_ms2 <- MetaboSpectra::clean_spMat(spMat_ms2)
spMat_ms2 <- MetaboSpectra::round_spMat(spMat_ms2, digits = 0)
standard_spMat <- MetaboSpectra::get_spMat(standard_car[1, ])
standard_spMat <- MetaboSpectra::clean_spMat(standard_spMat)
MetaboSpectra::compare_spMat_ndotproduct(x = spMat_ms2, y = standard_spMat, joinpeak = "outer")
MetaboSpectra::compare_spMat_entropy(x = spMat_ms2, y = standard_spMat)
MetaboSpectra::plotSpectra(spMat_ms2)
MetaboSpectra::plotSpectra(standard_spMat)

plot(x = ms1PeaksDT[i, rtime][[1]], y = ms1PeaksDT[i, intensity][[1]], type = "l",
     xlab = "Retention Time", ylab = "Intensity", col = "red", xlim = c(32, 38))
for(j in 1:nrow(ms2DT_i)){
  lines(ms2DT_i[j, rtime][[1]], ms2DT_i[j, intensity][[1]], col = "blue")
}

plot(x = ms1PeaksDT[i, rtime][[1]], y = ms1PeaksDT[i, intensity][[1]] / max(ms1PeaksDT[i, intensity][[1]]), type = "l",
     xlab = "Retention Time", ylab = "Intensity", col = "red", xlim = c(32, 38))
for(j in 1:nrow(ms2DT_i)){
  lines(ms2DT_i[j, rtime][[1]], ms2DT_i[j, intensity][[1]] / max(ms2DT_i[j, intensity][[1]]), col = "blue")
}
