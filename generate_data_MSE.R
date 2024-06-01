library(xcms)
library(MsExperiment)

# 1. Load data
dir_path <- "D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/MetaboDeconv/test_data/QC/"
file_path <- list.files(dir_path, pattern = ".mzML")
file_path <- paste0(dir_path, file_path)
pd <- data.frame(
  sample_name = sub(basename(file_path), pattern = ".mzML",
                    replacement = "", fixed = TRUE),
  sample_group = rep("UNKNOWN", 3),
  sample_type = c(rep("QC", 3)),
  inject_order = c(3, 1, 2),
  sample_path = file_path,
  stringsAsFactors = FALSE
)
pd <- pd %>%
  dplyr::arrange(inject_order)
data_MSE <- readMsExperiment(pd$sample_path, sampleData = pd)

load("./test_data/QC/best_params.RData")
cwp <- CentWaveParam(snthresh = best_params$snthresh, noise = best_params$noise, ppm = best_params$ppm,
                     peakwidth = c(best_params$min_peakwidth, best_params$max_peakwidth),
                     prefilter = c(best_params$prefilter, best_params$value_of_prefilter),
                     mzdiff = best_params$mzdiff)
data_MSE <- findChromPeaks(data_MSE, param = cwp,
                           chunkSize = length(data_MSE),
                           BPPARAM = BiocParallel::SnowParam(workers = length(data_MSE), type = "SOCK"))
pdp_subs <- PeakDensityParam(sampleGroups = sampleData(data_MSE)$sample_type,
                             minFraction = best_params$minFraction)
data_MSE <- groupChromPeaks(data_MSE, param = pdp_subs)
pgp_subs <- PeakGroupsParam(
  minFraction = best_params$minFraction,
  span = best_params$span,
  smooth = best_params$smooth,
  family = best_params$family,
  subset = c(1,3),
  subsetAdjust = "average")
data_MSE <- adjustRtime(data_MSE, param = pgp_subs,
                          chunkSize = length(data_MSE),
                        BPPARAM = BiocParallel::SnowParam(workers = length(data_MSE), type = "SOCK"))

pdp <- PeakDensityParam(sampleGroups = sampleData(data_MSE)$sample_type,
                        minFraction = best_params$minFraction,
                        bw = best_params$bw)
data_MSE <- groupChromPeaks(data_MSE, param = pdp)
data_MSE <- fillChromPeaks(data_MSE, param = ChromPeakAreaParam(),
                             chunkSize = length(data_MSE),
                             BPPARAM = BiocParallel::SnowParam(workers = length(data_MSE), type = "SOCK"))
data_MSE <- findChromPeaksIsolationWindow(data_MSE,
                                          param = cwp,
                                          chunkSize = length(data_MSE),
                                          BPPARAM = BiocParallel::SnowParam(workers = length(data_MSE), type = "SOCK"))
save(data_MSE, file = "./test_data/QC/data_MSE.RData")
