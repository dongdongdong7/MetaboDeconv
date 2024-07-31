# devtools::document()
# load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/swath_data.RData")
# load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/swath_spectra.RData")
# chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
#                                          xcms::chromPeakData(swath_data)),
#                                    rownames = "cpid")
# chromPeakTable_ms1 <- chromPeakTable %>%
#   dplyr::filter(ms_level == 1)
# chromPeakTable_ms2 <- chromPeakTable %>%
#   dplyr::filter(ms_level == 2)
# # filter chromPeakTable
# chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
#   dplyr::filter(maxo >= 1000)
# chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
#   dplyr::filter(maxo > 100)
# chromPeakTable <- rbind(chromPeakTable_ms1, chromPeakTable_ms2)
# chromPeaks_new <- as.data.frame(chromPeakTable[, 2:12])
# rownames(chromPeaks_new) <- chromPeakTable$cpid
# xcms::chromPeaks(swath_data) <- as.matrix(chromPeaks_new)
# chromPeakData_new <- as.data.frame(chromPeakTable[, 13:ncol(chromPeakTable)])
# rownames(chromPeakData_new) <- chromPeakTable$cpid
# xcms::chromPeakData(swath_data) <- chromPeakData_new
# chrDfList <- getChromPeaksDf(ndata = swath_data, cpid = chromPeakTable$cpid, cpN = 100, noise1 = 100, noise2 = 10, thread = 3)
# clusterPeaks0 <- cluster_peak(ndata = swath_data,chrDfList = chrDfList, cpid = "CP34", factor = 2, method = "direct",noise_threshold = 0.01,
#                              cosTh = -1, corTh = -1)
# plot_chrDfList(clusterPeaks0)
# clusterPeaks <- cluster_peak(ndata = swath_data,chrDfList = chrDfList, cpid = "CP34", factor = 1, method = "direct",noise_threshold = 0.01,
#                              cosTh = 0.9, corTh = 0.9)
# plot_chrDfList(clusterPeaks)
# sp <- peak2spectra(clusterPeaks)
# Spectra::plotSpectra(sp)
# chromPeakTable_ms1 <- Deconv4ndata(ndata = swath_data, factor = 1,cosTh = 0.8, corTh = 0.8,noise1 = 100, noise2 = 10, noise_threshold = 0.01, method = "direct", thread = 3, cpN = 100)
# DIA_spMat <- sp2spMat(chromPeakTable_ms1[9, ]$spectra[[1]])
# DIA_spMat1 <- MetaboSpectra::clean_spMat(DIA_spMat)
# DIA_spMat2 <- MetaboSpectra::clean_spMat(DIA_spMat, normalize_intensity = TRUE)
# MetaboSpectra::plotSpectra(DIA_spMat1)
# fenamiphos <- Spectra::Spectra(
#       system.file("mgf", "metlin-72445.mgf", package = "xcms"),
#       source = MsBackendMgf::MsBackendMgf())
# fenamiphos_spMat <- sp2spMat(fenamiphos[2])
# fenamiphos_spMat1 <- MetaboSpectra::clean_spMat(fenamiphos_spMat, noise_threshold = 0.01)
# fenamiphos_spMat2 <- MetaboSpectra::clean_spMat(fenamiphos_spMat, noise_threshold = 0.01, normalize_intensity = TRUE)
# MetaboSpectra::plotSpectra(fenamiphos_spMat1)
# MetaboSpectra::plotComparableSpectra(DIA_spMat1, fenamiphos_spMat1, num = 30, tol_da2 = 0.05)
# MetaboSpectra::compare_spMat_entropy(DIA_spMat2,fenamiphos_spMat2)
# MetaboSpectra::compare_spMat_ndotproduct(DIA_spMat1,fenamiphos_spMat1, joinpeak = "outer")
# fenamiphos_mz <- 304.113077
# fenamiphos_ms1_peak <- xcms::chromPeaks(swath_data, mz = fenamiphos_mz, ppm = 2)
# fenamiphos_ms1_peak
# fenamiphos_swath_spectrum <- swath_spectra[
#   swath_spectra$peak_id == rownames(fenamiphos_ms1_peak)]
# fenamiphos_swath_spMat <- sp2spMat(fenamiphos_swath_spectrum)
# fenamiphos_swath_spMat1 <- MetaboSpectra::clean_spMat(fenamiphos_swath_spMat, noise_threshold = 0.01)
# fenamiphos_swath_spMat2 <- MetaboSpectra::clean_spMat(fenamiphos_swath_spMat, noise_threshold = 0.01, normalize_intensity = TRUE)
# MetaboSpectra::plotSpectra(fenamiphos_swath_spMat1)
# MetaboSpectra::plotComparableSpectra(fenamiphos_swath_spMat1, fenamiphos_spMat1, num = 30, tol_da2 = 0.05)
# MetaboSpectra::compare_spMat_entropy(fenamiphos_swath_spMat2,fenamiphos_spMat2)
# MetaboSpectra::compare_spMat_ndotproduct(fenamiphos_swath_spMat1,fenamiphos_spMat1, joinpeak = "outer")

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
calCor_shape <- function(chrDf1, chrDf2, method = "direct"){
  if(method == "apex"){
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
    intTmp <- align_vectors(A = chrDf1$intensity, B = chrDf2$intensity)
    int1 <- intTmp[[1]];int2 <- intTmp[[2]]
  }else if(method == "direct"){
    int1 <- chrDf1$intensity
    int2 <- chrDf2$intensity
  }else stop("method wrong!")
  cor_shape <- round(cor(int1, int2, method = "pearson", use = "pairwise.complete.obs"), 4)
  if(is.na(cor_shape)) cor_shape <- 0
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
cosine_similarity <- function(vec1, vec2) {
  vec1[is.na(vec1)] <- 0
  vec2[is.na(vec2)] <- 0

  # 计算点积
  dot_product <- sum(vec1 * vec2)

  # 计算范数
  norm_vec1 <- sqrt(sum(vec1^2))
  norm_vec2 <- sqrt(sum(vec2^2))

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
approx_chrDf <- function(chrDf_x, chrDf_y){
  tmp <- approx(x = chrDf_x$rt, y = chrDf_x$intensity, xout = chrDf_y$rt)
  extra_chrDf <- dplyr::tibble(intensity = tmp$y, rt = tmp$x)
  extra_chrDf$mz <- mean(chrDf_x$mz)
  chrDf_x_new <- rbind(chrDf_x, extra_chrDf) %>% dplyr::arrange(rt)

  tmp <- approx(x = chrDf_y$rt, y = chrDf_y$intensity, xout = chrDf_x$rt)
  extra_chrDf <- dplyr::tibble(intensity = tmp$y, rt = tmp$x)
  extra_chrDf$mz <- mean(chrDf_y$mz)
  chrDf_y_new <- rbind(chrDf_y, extra_chrDf) %>% dplyr::arrange(rt)

  return(list(chrDf_x = chrDf_x_new, chrDf_y = chrDf_y_new))
}
getChromPeaksDf <- function(ndata, cpid = NA, cpN = NA, noise1 = 0, noise2 = 0, smooth = FALSE, size = 3, expandRt = 0, expandMz = 0, thread = 1){
  #browser()
  if(any(is.na(cpid))) stop("cpid has NA!")
  if(is.na(cpN)){
    tmp <- xcms::chromPeakChromatograms(ndata, peaks = cpid, expandRt = expandRt, expandMz = expandMz)
    tmp <- lapply(1:nrow(tmp), function(j) {tmp[j, 1]})
    names(tmp) <- cpid
    return(tmp)
  }
  else{
    blocks <- lapply(seq(1,length(cpid),cpN), function(start) {
      end <- min(start + cpN - 1, length(cpid))
      c(start,end)
    })
    loop <- function(i){
      start <- blocks[[i]][1];end <- blocks[[i]][2]
      tmp <- xcms::chromPeakChromatograms(ndata, peaks = cpid[start:end], expandRt = expandRt, expandMz = expandMz, progressbar = FALSE)
      tmpList <- lapply(1:nrow(tmp), function(j) {tmp[j, 1]})
      names(tmpList) <- cpid[start:end]
      return(tmpList)
    }
    pb <- utils::txtProgressBar(max = length(blocks), style = 3)
    if(thread == 1){
      tmp <- lapply(1:length(blocks), function(i) {
        utils::setTxtProgressBar(pb, i)
        loop(i)
      })
      tmp <- purrr::list_flatten(tmp)
    }else if(thread > 1){
      cl <- snow::makeCluster(thread)
      doSNOW::registerDoSNOW(cl)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                   n))
      tmp <- foreach::`%dopar%`(foreach::foreach(i = 1:length(blocks),
                                                 .packages = c("xcms"),
                                                 .options.snow = opts),
                                {
                                  loop(i)
                                })
      snow::stopCluster(cl)
      gc()
      tmp <- purrr::list_flatten(tmp)
    }
  }
  chrDfList <- lapply(1:length(tmp), function(i) {
    XChr <- tmp[[i]]
    msLevel <- XChr@msLevel
    if(msLevel == 1) noise <- noise1
    else if(msLevel == 2) noise <- noise2
    else stop("msLevel wrong!")
    intensity <- XChr@intensity
    intensity[which(intensity <= noise)] <- NA
    rt <- XChr@rtime
    chrDf <- dplyr::tibble(intensity = intensity, rt = rt) %>%
      dplyr::filter(!is.na(intensity))
    chrDf$mz <- mean(XChr@mz)
    if(smooth) chrDf$intensity <- smoothMean(chrDf$intensity, size = size)
    return(chrDf)
  })
  names(chrDfList) <- names(tmp)
  return(chrDfList)
}
cluster_peak <- function(ndata = ndata, chrDfList, cpid = NA, factor = 1, cosTh = 0.8, corTh = 0.8, method = "direct", noise_threshold = 0.05){
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
    dplyr::filter(maxo > max(mchromPeakTable_ms2$maxo) * noise_threshold)
  chrDfList_ms1 <- chrDfList[mchromPeakTable_ms1$cpid]
  if(nrow(mchromPeakTable_ms2) == 0) return(NULL)
  else{
    chrDfList_ms2 <- chrDfList[mchromPeakTable_ms2$cpid]
    chrDfLength <- sapply(chrDfList_ms2, function(x) {
      if(is.null(x)) return(0)
      else{
        return(nrow(x))
      }
    })
    chrDfList_ms2 <- chrDfList_ms2[chrDfLength!=0]
  }
  if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
  #browser()
  scoreList <- lapply(1:length(chrDfList_ms2), function(i) {
    #plot_chrDfList(list(ms1 = chrDfList_ms1[1], ms2 = chrDfList_ms2[i]))
    res <- approx_chrDf(chrDf_x = chrDfList_ms1[[1]], chrDf_y = chrDfList_ms2[[i]])
    #plot_chrDfList(list(ms1 = list(res$chrDf_x), ms2 = list(res$chrDf_y)))
    cosScore <- cosine_similarity(res$chrDf_x$intensity, res$chrDf_y$intensity)
    corScore <- calCor_shape(chrDf1 = res$chrDf_x, chrDf2 = res$chrDf_y, method = method)
    return(c(cosScore = cosScore, corScore = corScore))
  })
  cosScore <- sapply(scoreList, function(x) {x["cosScore"]})
  corScore <- sapply(scoreList, function(x) {x["corScore"]})
  chrDfList_ms2 <- chrDfList_ms2[which(cosScore > cosTh & corScore > corTh)]
  if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
  chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
  #plot_chrDfList(chrDfList)
  return(chrDfList)
}
peak2spectra <- function(clusterPeaks){
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
#' @param noise_threshold noise threshold.
#' @param method calCor_shape parameter, for align method.
#' @param cpN cpN is getChromPeaksDf's para for parallel.
#'
#' @return A tibble with spectra.
#' @export
#'
#' @examples
#' load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/swath_data.RData")
#' load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/generate_data/test_data/swath_spectra.RData")
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
#' chrDfList <- getChromPeaksDf(ndata = swath_data, cpid = chromPeakTable$cpid, cpN = 100, noise1 = 100, noise2 = 10, thread = 3)
#' clusterPeaks <- cluster_peak(ndata = swath_data,chrDfList = chrDfList, cpid = "CP34", factor = 1, method = "direct",noise_threshold = 0.01,
#'                              cosTh = 0.9, corTh = 0.9)
#' sp <- peak2spectra(clusterPeaks)
#' Spectra::plotSpectra(sp)
#' chromPeakTable_ms1 <- Deconv4ndata(ndata = swath_data, factor = 1,cosTh = 0.8, corTh = 0.8,noise1 = 100, noise2 = 10, noise_threshold = 0.01, method = "direct", thread = 1)
#' DIA_spMat <- sp2spMat(chromPeakTable_ms1[9, ]$spectra[[1]])
#' DIA_spMat1 <- MetaboSpectra::clean_spMat(DIA_spMat)
#' DIA_spMat2 <- MetaboSpectra::clean_spMat(DIA_spMat, normalize_intensity = TRUE)
#' MetaboSpectra::plotSpectra(DIA_spMat1)
#' fenamiphos <- Spectra::Spectra(
#'       system.file("mgf", "metlin-72445.mgf", package = "xcms"),
#'       source = MsBackendMgf::MsBackendMgf())
#' fenamiphos_spMat <- sp2spMat(fenamiphos[2])
#' fenamiphos_spMat1 <- MetaboSpectra::clean_spMat(fenamiphos_spMat, noise_threshold = 0.01)
#' fenamiphos_spMat2 <- MetaboSpectra::clean_spMat(fenamiphos_spMat, noise_threshold = 0.01, normalize_intensity = TRUE)
#' MetaboSpectra::plotSpectra(fenamiphos_spMat1)
#' MetaboSpectra::plotComparableSpectra(DIA_spMat1, fenamiphos_spMat1, num = 30, tol_da2 = 0.05)
#' MetaboSpectra::compare_spMat_entropy(DIA_spMat2,fenamiphos_spMat2)
#' MetaboSpectra::compare_spMat_ndotproduct(DIA_spMat1,fenamiphos_spMat1, joinpeak = "inner")
Deconv4ndata <- function(ndata, smooth = FALSE, size = 3, factor = 0.5, noise1 = 0, noise2 = 0, cosTh = 0.8, corTh = 0.8, noise_threshold = 0.05, method = "direct",thread = 1, cpN = NA){
  chromPeakTable <- chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
                                                             xcms::chromPeakData(ndata)),
                                                       rownames = "cpid")
  chromPeakTable_ms1 <- chromPeakTable %>%
    dplyr::filter(ms_level == 1)
  #browser()
  chrDfList <- getChromPeaksDf(ndata = ndata, cpid = chromPeakTable$cpid,cpN = cpN, noise1 = noise1, noise2 = noise2, smooth = smooth, size = size, thread = thread)
  loop <- function(x){
    clusterPeaks <- cluster_peak(ndata = ndata, chrDfList = chrDfList, cpid = x, factor = factor, cosTh = cosTh, corTh = corTh, method = method, noise_threshold = noise_threshold)
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
      loop(x)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    envir <- environment(cluster_peak)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(smoothMean)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(peak2spectra)
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
