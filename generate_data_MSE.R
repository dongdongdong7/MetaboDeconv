library(magrittr)
# 创建一个用于测试的MSE数据
# 如何减少feature的数量, CPC
# 1. Load data
dir_path <- "D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/QC/"
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
data <- MsExperiment::readMsExperiment(pd$sample_path, sampleData = pd)

# Peak picking
optPara_ms1_xcms <- readRDS("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/optPara_ms1_xcms.rds")
optPara_ms2_xcms <- readRDS("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/optPara_ms2_xcms.rds")
cwp_opt_ms1 <- xcms::CentWaveParam(snthresh = optPara_ms1_xcms["signal/noise threshold"],
                                   noise = optPara_ms1_xcms["noise filter"],
                                   ppm = optPara_ms1_xcms["ppm"],
                                   peakwidth = c(optPara_ms1_xcms["minimum peakwidth"], optPara_ms1_xcms["maximum peakwidth"]),
                                   prefilter = c(optPara_ms1_xcms["prefilter peaks"], optPara_ms1_xcms["noise filter"]))
cwp_opt_ms2 <- xcms::CentWaveParam(snthresh = optPara_ms2_xcms["signal/noise threshold"],
                                   noise = optPara_ms2_xcms["noise filter"],
                                   ppm = optPara_ms2_xcms["ppm"],
                                   peakwidth = c(optPara_ms2_xcms["minimum peakwidth"], optPara_ms2_xcms["maximum peakwidth"]),
                                   prefilter = c(optPara_ms2_xcms["prefilter peaks"], optPara_ms2_xcms["noise filter"]))
data_opt_ms1 <- xcms::findChromPeaks(data, param = cwp_opt_ms1, chunkSize = 3L, msLevel = 1L,
                                     BPPARAM = BiocParallel::SnowParam(workers = 3L, type = "SOCK"))

# Peak shape filter
neatms_res <- readr::read_csv("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/neatms_output_ms1.csv")
colnames(neatms_res) <- c("row", "row_id", "sample", "mz", "rt", "height", "area", "label")
neatms_res$row_id <- neatms_res$row_id + 1
neatms_res <- neatms_res %>% dplyr::filter(label == "High_quality")
chromPeakTable_ms1 <- MetaboProcess::getChromPeakTable(data = data_opt_ms1, style = "xcms")
chromPeakTable_ms1 <- chromPeakTable_ms1[neatms_res$row_id, ]
data_filter <- data_opt_ms1
chromPeaks_new <- as.data.frame(chromPeakTable_ms1[, 2:12])
rownames(chromPeaks_new) <- chromPeakTable_ms1$cpid
xcms::chromPeaks(data_filter) <- as.matrix(chromPeaks_new)
chromPeakData_new <- as.data.frame(chromPeakTable_ms1[, 13:14])
rownames(chromPeakData_new) <- chromPeakTable_ms1$cpid
xcms::chromPeakData(data_filter) <- chromPeakData_new

# Peak annotation
resTable <- readRDS("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/resTable.rds")
resTable_filter <- resTable %>%
  dplyr::filter(isotope == "M0") %>%
  dplyr::filter(an1 == "[M+H]+" | is.na(an1))
chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
  dplyr::filter(cpid %in% resTable_filter$cpid)
chromPeaks_new <- as.data.frame(chromPeakTable_ms1[, 2:12])
rownames(chromPeaks_new) <- chromPeakTable_ms1$cpid
xcms::chromPeaks(data_filter) <- as.matrix(chromPeaks_new)
chromPeakData_new <- as.data.frame(chromPeakTable_ms1[, 13:14])
rownames(chromPeakData_new) <- chromPeakTable_ms1$cpid
xcms::chromPeakData(data_filter) <- chromPeakData_new

# MS2 assign
data_filter <- xcms::findChromPeaksIsolationWindow(data_filter,
                                          param = cwp_opt_ms2,
                                          chunkSize = length(data_filter),
                                          BPPARAM = BiocParallel::SnowParam(workers = length(data_filter), type = "SOCK"))
saveRDS(data_filter, file = "D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/data_filter.rds")
chromPeakTable <- MetaboProcess::getChromPeakTable(data = data_filter, style = "xcms")
chromPeakTable_ms1 <- chromPeakTable %>%
  dplyr::filter(ms_level == 1)
chromPeakTable_ms2 <- chromPeakTable %>%
  dplyr::filter(ms_level == 2)

## Generate a chrDf
test <- generate_chrDf(ndata = data_filter[1], cpid = "CP105468")
MetaboProcess::plot_chrDf(test)

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

chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(data_MSE),
                                         xcms::chromPeakData(data_MSE)),
                                   rownames = "cpid")
chromPeakTable_ms1 <- chromPeakTable %>%
  dplyr::filter(ms_level == 1)

featureInfo <- xcms::featureDefinitions(data_MSE) %>%
  dplyr::as_tibble(rownames = "ftid")# 13023
featureValue <- dplyr::as_tibble(xcms::featureValues(data_MSE, value = "into"))
featureTable <- dplyr::as_tibble(cbind(featureInfo, featureValue))
cpidList <- lapply(1:nrow(featureTable), function(i) {
  peakidx <- featureTable[i, ]$peakidx[[1]]
  return(chromPeakTable_ms1[peakidx, ]$cpid)
})
featureTable$cpid <- cpidList

featureTable[14, ]
featureTable[14, ]$cpid
chromPeakTable[chromPeakTable$cpid %in% featureTable[14, ]$cpid[[1]], ]
featureTable[14, c(12,13,14)]
plot_chrDf(generate_chrDf(ndata = data_MSE, cpid = "CP00332", smooth = TRUE, size = 5))

data_xcmsSet <- xcms:::.XCMSnExp2xcmsSet(data_MSE)
xsa <- CAMERA::xsAnnotate(data_xcmsSet)
xsaF <- CAMERA::groupFWHM(xsa, perfwhm=0.6)
xsaC <- CAMERA::groupCorr(xsaF)
xsaFI <- CAMERA::findIsotopes(xsaC)
#With provided rule table
file <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
rules <- read.csv(file)

xsaFA <- CAMERA::findAdducts(xsaFI, polarity="positive", rules = rules)
peakList <- CAMERA::getPeaklist(xsaFA) # 13023
rownames(peakList) <- rownames(data_peakList)
