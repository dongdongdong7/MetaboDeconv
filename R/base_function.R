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
calCor_shape <- function(chrDf1, chrDf2){
  align_vectors <- function(A, B) {
    # find apex
    peak_A <- which.max(A)
    peak_B <- which.max(B)

    # length of two vector
    left_len <- max(peak_A, peak_B) - 1
    right_len <- max(length(A) - peak_A, length(B) - peak_B)

    # new length
    aligned_length <- left_len + right_len + 1
    new_A <- rep(NA, aligned_length)
    new_B <- rep(NA, aligned_length)

    # alignment A
    start_A <- left_len - (peak_A - 1) + 1
    end_A <- start_A + length(A) - 1
    new_A[start_A:end_A] <- A

    # alignment B
    start_B <- left_len - (peak_B - 1) + 1
    end_B <- start_B + length(B) - 1
    new_B[start_B:end_B] <- B

    return(list(aligned_A = new_A, aligned_B = new_B))
  }
  int1 <- chrDf1$intensity
  int2 <- chrDf2$intensity
  if(all(diff(int2) > 0) | all(diff(int2) < 0)){
    warning("This feature is detected as incremental or decremental!")
    return(0)
  }
  intTmp <- align_vectors(A = int1, B = int2)
  int1 <- intTmp[[1]];int2 <- intTmp[[2]]
  cor_shape <- round(cor(int1, int2, method = "pearson", use = "pairwise.complete.obs"), 4)
  return(cor_shape)
}
sp2spMat <- function(sp){
  spMat <- Spectra::peaksData(sp)[[1]]
  if(nrow(spMat) != 1) spMat <- spMat[order(spMat[, 1]), ]
  return(spMat)
}
spMat2sp <- function(spMat, msLevel = 1L, rtime = 0){
  spd <- S4Vectors::DataFrame(msLevel = msLevel, rtime = rtime)
  spd$mz <- list(spMat[, "mz"])
  spd$intensity <- list(spMat[,"intensity"])
  return(Spectra::Spectra(spd))
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

generate_chrDf <- function(ndata, cpid = NA, noise = 0,smooth = TRUE, size = 3){
  if(is.na(cpid)) stop("Please set cpid!")
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
  #browser()
  chrDf <- purrr::list_rbind(lapply(1:length(cpSpectra), function(j) {
    peakMat <- peaksData[[j]]
    if(nrow(peakMat) == 0){
      mz <- NA;intensity <- 0
      return(data.frame(mz = mz, intensity = intensity, row.names = NULL))
    }else{
      idx <- which.max(peakMat[, "intensity"])
      #mz <- mean(peakMat[, "mz"]);intensity <- mean(peakMat[, "intensity"])
      mz <- peakMat[, "mz"][idx];intensity <- peakMat[, "intensity"][idx]
      peakDf <- data.frame(mz = mz, intensity = intensity, row.names = NULL)
      return(peakDf)
    }
  }))
  chrDf$rt <- rtSpectra
  chrDf <- chrDf %>%
    dplyr::filter(intensity > noise)
  if(smooth) chrDf$intensity <- smoothMean(chrDf$intensity, size = size)
  return(chrDf)
}
cluseter_peak <- function(ndata = ndata, cpid = NA, factor = 0.5, noise1 = 0, noise2 = 0, smooth = TRUE, size = 3, cosTh = 0.8, corTh = 0.8){
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
  if(nrow(mchromPeakTable_ms2) == 0) return(NULL)
  else{
    chrDfList_ms2 <- lapply(mchromPeakTable_ms2$cpid, function(x) {
      chrDf_x <- generate_chrDf(ndata = ndata, cpid = x, smooth = smooth, size = size, noise = noise2)
    })
    names(chrDfList_ms2) <- mchromPeakTable_ms2$cpid
  }
  chrDf_ms1 <- generate_chrDf(ndata = ndata, cpid = mchromPeakTable_ms1$cpid,
                              smooth = smooth, size = size, noise = noise1)
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

#' @title Deconv4ndata
#' @description
#' Deconv function for nth sample.
#'
#' @param ndata nth sample, a XcmsExperiment object.
#' @param smooth whether or not smooth.
#' @param size size smoothMean size.
#' @param factor factor * peakwidth -> rtRange.
#' @param cosTh cos threshold.
#' @param corTh cor threshold.
#' @param thread thread.
#' @param noise1 noise for ms1 peak.
#' @param noise2 noise for ms2 peak.
#'
#' @return A tibble with spectra.
#' @export
#'
#' @examples
#' load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/swath_data.RData")
#' chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
#'                                          xcms::chromPeakData(swath_data)),
#'                                    rownames = "cpid")
#' chromPeakTable_ms1 <- chromPeakTable %>%
#'   dplyr::filter(ms_level == 1)
#' chromPeakTable_ms2 <- chromPeakTable %>%
#'   dplyr::filter(ms_level == 2)
#' # filter chromPeakTable
#' chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
#'   dplyr::filter(maxo >= 1000)
#' chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
#'   dplyr::filter(maxo > 100)
#' chromPeakTable <- rbind(chromPeakTable_ms1, chromPeakTable_ms2)
#' chromPeaks_new <- as.data.frame(chromPeakTable[, 2:12])
#' rownames(chromPeaks_new) <- chromPeakTable$cpid
#' xcms::chromPeaks(swath_data) <- as.matrix(chromPeaks_new)
#' chromPeakData_new <- as.data.frame(chromPeakTable[, 13:ncol(chromPeakTable)])
#' rownames(chromPeakData_new) <- chromPeakTable$cpid
#' xcms::chromPeakData(swath_data) <- chromPeakData_new
#' chromPeakTable_ms1 <- Deconv4ndata(ndata = swath_data, thread = 3, cosTh = 0.6, corTh = 0.6,noise1 = 100, noise2 = 10)
#' Spectra::plotSpectra(chromPeakTable_ms1[9, ]$spectra[[1]])
Deconv4ndata <- function(ndata, smooth = TRUE, size = 3, factor = 0.5, noise1 = 0, noise2 = 0, cosTh = 0.8, corTh = 0.8, thread = 1){
  chromPeakTable <- chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                                             xcms::chromPeakData(ndata)),
                                                       rownames = "cpid")
  chromPeakTable_ms1 <- chromPeakTable %>%
    dplyr::filter(ms_level == 1)
  chromPeakTable_ms2 <- chromPeakTable %>%
    dplyr::filter(ms_level == 2)
  #browser()
  loop <- function(x){
    clusterPeaks <- cluseter_peak(ndata = ndata, cpid = x, noise1 = noise1, noise2 = noise2, cosTh = cosTh, corTh = corTh)
    if(is.null(clusterPeaks) | is.null(clusterPeaks$ms2) | length(clusterPeaks$ms2) == 0) return(NULL)
    sp <- peak2spectra(clusterPeaks)
    return(sp)
  }

  pb <- utils::txtProgressBar(max = nrow(chromPeakTable_ms1), style = 3)
  #browser()
  if(thread == 1){
    spectraList <- lapply(1:nrow(chromPeakTable_ms1), function(m) {
      utils::setTxtProgressBar(pb, m)
      x <- chromPeakTable_ms1[m, ]$cpid
      #print(m)
      loop(x)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    envir <- environment(cluseter_peak)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(smoothMean)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(peak2spectra)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(generate_chrDf)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    spectraList <- foreach::`%dopar%`(foreach::foreach(m = 1:nrow(chromPeakTable_ms1),
                                                       .packages = c("dplyr"),
                                                       .options.snow = opts),
                                      {
                                        x <- chromPeakTable_ms1[m, ]$cpid
                                        loop(x)
                                      })
    snow::stopCluster(cl)
    gc()
  }else stop("thread is wrong!")
  chromPeakTable_ms1$spectra <- spectraList
  return(chromPeakTable_ms1)
}
