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
data_opt_ms1 <- xcms::findChromPeaks(data, param = cwp_opt_ms1, chunkSize = 1L, msLevel = 1L,
                                     BPPARAM = BiocParallel::SnowParam(workers = 1L, type = "SOCK"))
chromPeakTable_ms1 <- MetaboProcess::getChromPeakTable(data = data_opt_ms1, style = "xcms")
i <- 600
tmp <- xcms::chromPeakChromatograms(data_opt_ms1, peaks =  chromPeakTable_ms1$cpid[i], aggregationFun = "sum")
max(tmp[1,1]@intensity, na.rm = TRUE)
xcms::plot(tmp)
chromPeakTable_ms1[i, ]$maxo

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
data_filter <- readRDS("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/240722/data_filter.rds")
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
test <- generate_chrDf(ndata = data_filter[1], cpid = "CP34")
MetaboProcess::plot_chrDf(test)

## cluster peak
cluseter_peak <- function(ndata = ndata, cpid = NA, factor = 0.5, smooth = TRUE, size = 3, cosTh = 0.8, corTh = 0.8){
  if(is.na(cpid)) stop("Please set cpid!")
  chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                           xcms::chromPeakData(ndata)),
                                     rownames = "cpid")
  chromPeakTable_ms1 <- chromPeakTable %>% dplyr::filter(ms_level == 1)
  chromPeakTable_ms2 <- chromPeakTable %>% dplyr::filter(ms_level == 2)
  mchromPeakTable_ms1 <- chromPeakTable[chromPeakTable$cpid == cpid, ]
  if(mchromPeakTable_ms1$ms_level != 1) stop("msLevel should be 1!")
  peakWidth <- mchromPeakTable_ms1$rtmax - mchromPeakTable_ms1$rtmin
  rtRange <- c(mchromPeakTable_ms1$rt - factor * peakWidth,
               mchromPeakTable_ms1$rt + factor * peakWidth)
  mchromPeakTable_ms2 <- chromPeakTable_ms2 %>%
    dplyr::filter(rtmin >= rtRange[1] & rtmax <= rtRange[2]) %>%
    dplyr::filter(mchromPeakTable_ms1$mz >= isolationWindowLowerMz & mchromPeakTable_ms1$mz <= isolationWindowUpperMz)
  # remove ms2 feature mz >= precursorMz
  mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
    dplyr::filter(mzmax < mchromPeakTable_ms1$mzmin)
  # remove noise
  mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
    dplyr::filter(maxo > max(mchromPeakTable_ms2$maxo) * 0.05)
  #browser()
  if(nrow(mchromPeakTable_ms2) == 0) return(NULL)
  else{
    chrDfList_ms2 <- lapply(mchromPeakTable_ms2$cpid, function(x) {
      chrDf_x <- generate_chrDf(ndata = ndata, cpid = x, smooth = smooth, size = size)
    })
    names(chrDfList_ms2) <- mchromPeakTable_ms2$cpid
  }
  chrDf_ms1 <- generate_chrDf(ndata = ndata, cpid = mchromPeakTable_ms1$cpid,
                              smooth = smooth, size = size)
  chrDfList_ms1 <- list(chrDf_ms1)
  names(chrDfList_ms1) <- mchromPeakTable_ms1$cpid
  #browser()
  if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
  score <- sapply(1:length(chrDfList_ms2), function(i) {
    vec1_rt <- chrDf_ms1$rt;vec2_rt <- chrDfList_ms2[[i]]$rt
    vec1_int <- chrDf_ms1$intensity;vec2_int <- chrDfList_ms2[[i]]$intensity
    tmp <- align_vectors(vec1 = vec1_rt, vec2 = vec2_rt)
    vec1_rt <- tmp$aligned_vec1;vec2_rt <- tmp$aligned_vec2
    vec1_rt[!is.na(vec1_rt)] <- vec1_int;vec2_rt[!is.na(vec2_rt)] <- vec2_int
    vec1_int <- vec1_rt;vec2_int <- vec2_rt
    cosine_similarity(vec1_int, vec2_int)
  })
  chrDfList_ms2 <- chrDfList_ms2[which(score > cosTh)]
  if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
  score <- sapply(1:length(chrDfList_ms2), function(i) {
    calCor_shape(chrDf_ms1, chrDfList_ms2[[i]])
  })
  chrDfList_ms2 <- chrDfList_ms2[which(score > corTh)]
  chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
  return(chrDfList)
}
align_vectors <- function(vec1, vec2) {
  # 创建一个矩阵来存储两个向量之间的差值
  diff_matrix <- outer(vec1, vec2, FUN = function(x, y) abs(x - y))

  # 找到最小差值的位置
  min_indices <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
  best_i <- min_indices[1, 1]
  best_j <- min_indices[1, 2]

  # 计算对齐后的长度
  max_length <- max(length(vec1) + best_j - 1, length(vec2) + best_i - 1)

  # 对齐向量
  aligned_vec1 <- c(rep(NA, best_j - 1), vec1, rep(NA, max_length - (length(vec1) + best_j - 1)))
  aligned_vec2 <- c(rep(NA, best_i - 1), vec2, rep(NA, max_length - (length(vec2) + best_i - 1)))

  list(aligned_vec1 = aligned_vec1, aligned_vec2 = aligned_vec2)
}
cosine_similarity <- function(vec1, vec2) {
  # 计算点积
  dot_product <- sum(vec1 * vec2, na.rm = TRUE)

  # 计算范数
  norm_vec1 <- sqrt(sum(vec1^2, na.rm = TRUE))
  norm_vec2 <- sqrt(sum(vec2^2, na.rm = TRUE))

  # 计算余弦相似度
  cosine_sim <- dot_product / (norm_vec1 * norm_vec2)

  return(cosine_sim)
}
fill_na_with_mean <- function(vec) {
  # 获取向量的长度
  n <- length(vec)

  # 遍历向量中的每个元素
  for (i in 1:n) {
    # 如果当前元素是NA
    if (is.na(vec[i])) {
      # 找到前一个非NA值
      forward <- NA
      if (i > 1) {
        for (j in (i-1):1) {
          if (!is.na(vec[j])) {
            forward <- vec[j]
            break
          }
        }
      }

      # 找到后一个非NA值
      backward <- NA
      if (i < n) {
        for (j in (i+1):n) {
          if (!is.na(vec[j])) {
            backward <- vec[j]
            break
          }
        }
      }

      # 计算前后两个非NA值的平均值，并填充当前的NA
      if (!is.na(forward) && !is.na(backward)) {
        vec[i] <- (forward + backward) / 2
      } else if (!is.na(forward)) {
        vec[i] <- forward
      } else if (!is.na(backward)) {
        vec[i] <- backward
      }
    }
  }

  return(vec)
}
chromPeakTable_ms1[chromPeakTable_ms1$maxo > 10000, ][201, ]
cpid <- chromPeakTable_ms1[chromPeakTable_ms1$maxo > 10000, ][201, ]$cpid
test <- generate_chrDf(data_filter[1], cpid = cpid)
plot_chrDf(test)
system.time(clusterPeaks <- cluseter_peak(ndata = data_filter[1], cpid = cpid, cosTh = 0.4, corTh = 0.8))
clusterPeaks <- cluseter_peak(ndata = data_filter[1], cpid = cpid, cosTh = 0.8, corTh = 0.8)
plot_chrDfList(clusterPeaks)

## clusterPeak2spectra
peak2spectra <- function(clusterPeaks){
  #browser()
  peak_ms1 <- clusterPeaks$ms1[[1]]
  peak_ms1$mz <- fill_na_with_mean(peak_ms1$mz)
  idx_ms1 <- which.max(peak_ms1$intensity)
  rt_ms1 <- peak_ms1$rt[idx_ms1]
  mz_ms1 <- peak_ms1$mz[idx_ms1]
  int_ms1 <- peak_ms1$intensity[idx_ms1]
  peaks_ms2 <- clusterPeaks$ms2
  df_ms1 <- data.frame(mz = mz_ms1, intensity = int_ms1)
  tmp <- lapply(peaks_ms2, function(x) {
    x$mz <- fill_na_with_mean(x$mz)
    idx_ms2 <- which.max(x$intensity)
    mz_ms2 <- x$mz[idx_ms2]
    int_ms2 <- x$intensity[idx_ms2]
    return(data.frame(mz = mz_ms2, intensity = int_ms2))
  })
  spDf <- rbind(purrr::list_rbind(tmp), df_ms1) %>%
    dplyr::arrange(mz)
  spd <- dplyr::tibble(msLevel = 2L, rtime = rt_ms1)
  spd$mz <- list(spDf$mz)
  spd$intensity <- list(spDf$intensity)
  sp <- Spectra::Spectra(spd)
  return(sp)
}
sp <- peak2spectra(clusterPeaks)
Spectra::plotSpectra(sp)

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
