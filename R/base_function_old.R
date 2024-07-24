# #' @title smoothMean
# #' @description
# #' Smooth function using mean.
# #'
# #' @param int intensity.
# #' @param size mean size.
# #'
# #' @return intensity vector
# #' @export
# #'
# #' @examples
# #' intensity <- c(1701.196, 2121.159, 2401.849, 2043.653, 1542.964,  723.803)
# #' smoothMean(int = intensity)
# smoothMean <- function(int, size = 3){
#   if(size %% 2 == 0){
#     stop("size should be a singular!")
#   }
#   half_size <- floor(size / 2)
#   N <- length(int)
#   special_idx <- c(seq_len(half_size), rev((N:(N - half_size + 1))))
#   smooth_int <- sapply(1:N, function(i) {
#     if(i %in% special_idx){
#       return(int[i])
#     }
#     y <- sapply((i - half_size):(i + half_size), function(j){
#       int[j]
#     })
#     smooth_int <- mean(y)
#   })
#   return(smooth_int)
# }
# #' @title get_chrDfList
# #' @description
# #' Get all related ms2 features chromatogram information of the mth ms1 feature.
# #'
# #' @param ndata nth sample data.
# #' @param m mth feature.
# #' @param smooth Smoothing or not.
# #' @param noise noise.
# #' @param factor factor * peakwidth -> rtRange
# #' @param sn sn threshold.
# #'
# #' @return A list contains chrDf which belongs to mth ms1 feature.
# #' @export
# #'
# #' @examples
# #' load("./test_data/swath_data.RData")
# #' mchrDfList <- get_chrDfList(ndata = swath_data, m = 34)
# get_chrDfList <- function(ndata, m, smooth = TRUE, noise = 10, factor = 1.5, sn = 10, maxo = 1000){
#   spectra_ms1 <- xcms::spectra(ndata) %>%
#     Spectra::filterMsLevel(1L)
#   spectra_ms2 <- xcms::spectra(ndata) %>%
#     Spectra::filterMsLevel(2L)
#   featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
#                                          xcms::chromPeakData(ndata)),
#                                    rownames = "feature_id")
#   featureTable_ms1 <- featureTable %>%
#     dplyr::filter(ms_level == 1)
#   featureTable_ms2 <- featureTable %>%
#     dplyr::filter(ms_level == 2)
#
#   mfeatureTable_ms1 <- featureTable_ms1[m, ]
#   if(mfeatureTable_ms1$sn < sn | mfeatureTable_ms1$maxo < maxo) return(NULL)
#   peakWidth <- mfeatureTable_ms1$rtmax - mfeatureTable_ms1$rtmin
#   rtRange <- c(mfeatureTable_ms1$rt - factor * peakWidth,
#                mfeatureTable_ms1$rt + factor * peakWidth)
#   mfeatureTable_ms2 <- featureTable_ms2 %>%
#     dplyr::filter(rt >= rtRange[1] & rt <= rtRange[2]) %>%
#     dplyr::filter(mfeatureTable_ms1$mz >= isolationWindowLowerMz & mfeatureTable_ms1$mz <= isolationWindowUpperMz)
#   # remove ms2 feature mz >= precursorMz
#   mfeatureTable_ms2 <- mfeatureTable_ms2 %>%
#     dplyr::filter(mzmax < mfeatureTable_ms1$mzmin)
#   # remove noise
#   mfeatureTable_ms2 <- mfeatureTable_ms2 %>%
#     dplyr::filter(maxo > max(mfeatureTable_ms2$maxo) * 0.05)
#   if(nrow(mfeatureTable_ms2) == 0){
#     chrDfList_ms2 <- NULL
#   }else{
#     chrDfList_ms2 <- lapply(1:nrow(mfeatureTable_ms2), function(i) {
#       mfeatureTable_ms2_i <- mfeatureTable_ms2[i, ]
#       mfeatureSpectra_ms2_i <- spectra_ms2 %>%
#         Spectra::filterRt(c(mfeatureTable_ms2_i$rtmin, mfeatureTable_ms2_i$rtmax)) %>%
#         Spectra::filterMzRange(mz = c(mfeatureTable_ms2_i$mzmin, mfeatureTable_ms2_i$mzmax), keep = TRUE)
#       rtSpectra_i <- Spectra::rtime(mfeatureSpectra_ms2_i)
#       peaksData_i <- Spectra::peaksData(mfeatureSpectra_ms2_i)
#       chrDf_i <- purrr::list_rbind(lapply(1:length(rtSpectra_i), function(j) {
#         rt_tmp <- rtSpectra_i[j]
#         peakMat <- peaksData_i[[j]]
#         if(nrow(peakMat) == 0){
#           return(NULL)
#         }else{
#           peakDf <- as.data.frame(peakMat)
#           peakDf$rt <- rt_tmp
#           return(peakDf)
#         }
#       }))
#       chrDf_i <- chrDf_i %>%
#         dplyr::filter(intensity > noise)
#       if(smooth) chrDf_i$smooth_int <- smoothMean(chrDf_i$intensity, size = 3)
#       else chrDf_i$smooth_int <- chrDf_i$intensity
#       return(chrDf_i)
#     })
#     names(chrDfList_ms2) <- mfeatureTable_ms2$feature_id
#   }
#   mfeatureSpectra_ms1 <- spectra_ms1 %>%
#     Spectra::filterRt(c(mfeatureTable_ms1$rtmin, mfeatureTable_ms1$rtmax)) %>%
#     Spectra::filterMzRange(mz = c(mfeatureTable_ms1$mzmin, mfeatureTable_ms1$mzmax), keep = TRUE)
#   rtSpectra <- Spectra::rtime(mfeatureSpectra_ms1)
#   peaksData <- Spectra::peaksData(mfeatureSpectra_ms1)
#   chrDf_ms1 <- purrr::list_rbind(lapply(1:length(rtSpectra), function(j) {
#     rt_tmp <- rtSpectra[j]
#     peakMat <- peaksData[[j]]
#     if(nrow(peakMat) == 0){
#       return(NULL)
#     }else{
#       peakDf <- as.data.frame(peakMat)
#       peakDf$rt <- rt_tmp
#       return(peakDf)
#     }
#   }))
#   chrDf_ms1 <- chrDf_ms1 %>%
#     dplyr::filter(intensity > noise)
#   if(smooth) chrDf_ms1$smooth_int <- smoothMean(chrDf_ms1$intensity, size = 3)
#   else chrDf_ms1$smooth_int <- chrDf_ms1$intensity
#   chrDfList_ms1 <- list(chrDf_ms1)
#   names(chrDfList_ms1) <- mfeatureTable_ms1$feature_id
#   chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
#   return(chrDfList)
# }
# #' @title calCor_rt
# #' @description
# #' Calculate the retention time correlation of two chrDf
# #'
# #' @param chrDf1 chrDf1
# #' @param chrDf2 chrDf2
# #'
# #' @return A double .
# #' @export
# #'
# #' @examples
# #' load("./test_data/swath_data.RData")
# #' mchrDfList <- get_chrDfList(ndata = swath_data, m = 34)
# #' calCor_rt(chrDf1 = mchrDfList$ms1[[1]], chrDf2 = mchrDfList$ms2[[5]])
# calCor_rt <- function(chrDf1, chrDf2){
#   deltaRt <- mean(c(diff(chrDf1$rt), diff(chrDf2$rt)))
#   rt1 <- mean(chrDf1$rt)
#   rt2 <- mean(chrDf2$rt)
#   sa <- abs(rt1 - rt2) / deltaRt
#   # exp(1) ^ -abs(rt1 - rt2)
#   # exp(1) ^ -sa
#   cor_rt <- round(exp(1) ^ -abs(rt1 - rt2), 4)
#   return(cor_rt)
# }
# #' @title calCor_shape
# #' @description
# #' Calculate the peak correlation of two chrDf
# #'
# #' @param chrDf1 chrDf1
# #' @param chrDf2 chrDf2
# #'
# #' @return A double.
# #' @export
# #'
# #' @examples
# #' calCor_shape(chrDf1 = mchrDfList$ms1[[1]], chrDf2 = mchrDfList$ms2[[5]])
# calCor_shape <- function(chrDf1, chrDf2){
#   align_vectors <- function(A, B) {
#     # find apex
#     peak_A <- which.max(A)
#     peak_B <- which.max(B)
#
#     # length of two vector
#     left_len <- max(peak_A, peak_B) - 1
#     right_len <- max(length(A) - peak_A, length(B) - peak_B)
#
#     # new length
#     aligned_length <- left_len + right_len + 1
#     new_A <- rep(NA, aligned_length)
#     new_B <- rep(NA, aligned_length)
#
#     # alignment A
#     start_A <- left_len - (peak_A - 1) + 1
#     end_A <- start_A + length(A) - 1
#     new_A[start_A:end_A] <- A
#
#     # alignment B
#     start_B <- left_len - (peak_B - 1) + 1
#     end_B <- start_B + length(B) - 1
#     new_B[start_B:end_B] <- B
#
#     return(list(aligned_A = new_A, aligned_B = new_B))
#   }
#   int1 <- chrDf1$intensity
#   int2 <- chrDf2$intensity
#   if(all(diff(int2) > 0) | all(diff(int2) < 0)){
#     warning("This feature is detected as incremental or decremental!")
#     return(0)
#   }
#   intTmp <- align_vectors(A = int1, B = int2)
#   int1 <- intTmp[[1]];int2 <- intTmp[[2]]
#   cor_shape <- round(cor(int1, int2, method = "pearson", use = "pairwise.complete.obs"), 4)
#   return(cor_shape)
# }
# #' @title chrDfList2spectra
# #' @description
# #' Transfor chrDfList to a Spectra object
# #'
# #' @param chrDfList chrDfList
# #' @param ndata nth sample data.
# #'
# #' @return A Spectra object.
# #' @export
# #'
# #' @examples
# #' chrDfList2spectra(chrDfList = mchrDfList_new, ndata = swath_data)
# chrDfList2spectra <- function(chrDfList, ndata){
#   featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
#                                          xcms::chromPeakData(ndata)),
#                                    rownames = "feature_id")
#   featureTable_ms1 <- featureTable %>%
#     dplyr::filter(ms_level == 1)
#   featureTable_ms2 <- featureTable %>%
#     dplyr::filter(ms_level == 2)
#   if(length(chrDfList$ms2) == 0) return(NULL)
#   featureName_ms1 <- names(chrDfList$ms1)
#   featureName_ms2 <- names(chrDfList$ms2)
#   mfeatureTable_ms1 <- featureTable_ms1[featureTable_ms1$feature_id == featureName_ms1, ]
#   newfeatureTable_ms2 <- featureTable_ms2[featureTable_ms2$feature_id %in% featureName_ms2, ]
#   spd <- dplyr::tibble(msLevel = 2L, rtime = mfeatureTable_ms1$rt)
#   spd$mz <- list(c(mfeatureTable_ms1$mz, newfeatureTable_ms2$mz))
#   spd$intensity <- list(c(mfeatureTable_ms1$maxo, newfeatureTable_ms2$maxo))
#   sp <- Spectra::Spectra(spd)
#   return(sp)
# }
# #' @title sp2spMat
# #' @description
# #' Transformer a Spectra object to spMat.
# #'
# #' @param sp A Spectra object
# #'
# #' @return A spMat
# #' @export
# #'
# #' @examples
# #' sp2spMat(sp = sp_ms2)
# sp2spMat <- function(sp){
#   spMat <- Spectra::peaksData(sp)[[1]]
#   if(nrow(spMat) != 1) spMat <- spMat[order(spMat[, 1]), ]
#   return(spMat)
# }
# #' @title spMat2sp
# #' @description
# #' Transfor spMat 2 sp.
# #'
# #' @param spMat A spMat.
# #'
# #' @return A Spectra object.
# #' @export
# #'
# #' @examples
# #' spMat2sp(sp_ms2_spMat, msLevel = 2L, rtime = 20)
# spMat2sp <- function(spMat, msLevel = 1L, rtime = 0){
#   spd <- S4Vectors::DataFrame(msLevel = msLevel, rtime = rtime)
#   spd$mz <- list(spMat[, "mz"])
#   spd$intensity <- list(spMat[,"intensity"])
#   return(Spectra::Spectra(spd))
# }
# #' @title Deconv4ndata_old
# #' @description
# #' Deconv function for nth sample.
# #'
# #'
# #' @param ndata nth sample data.
# #' @param weight_rt weight_rt.
# #' @param weight_shape weight_shape.
# #' @param st socre threshold.
# #' @param thread prarllel thread.
# #' @param chromPeakTable_ms1 1 level chromPeakTable which have been filtered.(must be generate by ndata)
# #' @param chromPeakTable_ms2 2 level chromPeakTable which have been filtered.(must be generate by ndata)
# #' @param noise noise value.
# #' @param smooth whether or not smooth.
# #' @param size smoothMean size.
# #' @param factor factor * peakwidth -> rtRange.
# #'
# #' @return A new featureTable_ms1 with ms2 spectra.
# #' @export
# #'
# #' @examples
# #' load("./test_data/swath_data.RData")
# #' chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
# #'                                         xcms::chromPeakData(swath_data)),
# #'                                   rownames = "cpid")
# #' chromPeakTable_ms1 <- chromPeakTable %>%
# #'  dplyr::filter(ms_level == 1)
# #' chromPeakTable_ms2 <- chromPeakTable %>%
# #'  dplyr::filter(ms_level == 2)
# #' # filter chromPeakTable
# #' chromPeakTable_ms1 <- chromPeakTable_ms1 %>%
# #'   dplyr::filter(maxo >= 1000)
# #' chromPeakTable_ms2 <- chromPeakTable_ms2 %>%
# #'   dplyr::filter(maxo > 100)
# #' chromPeakTable_ms1 <- Deconv4ndata(ndata = swath_data, chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2)
# Deconv4ndata_old <- function(ndata, chromPeakTable_ms1,chromPeakTable_ms2,
#                          noise = 10, smooth = TRUE, size = 3, factor = 0.5,
#                          weight_rt = 0.5, weight_shape = 0.5, st = 0.8,
#                          thread = 1){
#   loop <- function(x){
#     mchrDfList <- generate_chrDfList(ndata = ndata, cpid = x,
#                                      chromPeakTable_ms1 = chromPeakTable_ms1, chromPeakTable_ms2 = chromPeakTable_ms2,
#                                      noise = noise, smooth = smooth, size = size, factor = factor)
#     if(is.null(mchrDfList)) return(NULL)
#     mchrDfList_new <- filterChrDf(chrDfList = mchrDfList, weight_rt = weight_rt, weight_shape = weight_shape,st = st)
#     if(is.null(mchrDfList_new$ms2) | length(mchrDfList_new$ms2) == 0) return(NULL)
#     sp_ms2 <- chrDfList2spectra(chrDfList = mchrDfList_new, ndata = ndata)
#     return(sp_ms2)
#   }
#   pb <- utils::txtProgressBar(max = nrow(chromPeakTable_ms1), style = 3)
#   if(thread == 1){
#     spectraList <- lapply(1:nrow(chromPeakTable_ms1), function(m) {
#       utils::setTxtProgressBar(pb, m)
#       x <- chromPeakTable_ms1[m, ]$cpid
#       print(x)
#       loop(x)
#     })
#   }else if(thread > 1){
#     cl <- snow::makeCluster(thread)
#     doSNOW::registerDoSNOW(cl)
#     opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
#                                                                  n))
#     envir <- environment(generate_chrDfList)
#     parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
#     spectraList <- foreach::`%dopar%`(foreach::foreach(m = 1:nrow(chromPeakTable_ms1),
#                                                        .packages = c("dplyr"),
#                                                        .options.snow = opts),
#                                       {
#                                         x <- chromPeakTable_ms1[m, ]$cpid
#                                         loop(x)
#                                       })
#     snow::stopCluster(cl)
#     gc()
#   }else stop("thread is wrong!")
#   chromPeakTable_ms1$spectra <- spectraList
#   return(chromPeakTable_ms1)
# }
#
# cluseter_peak <- function(ndata = ndata, cpid = NA, factor = 0.5, noise = 0,smooth = TRUE, size = 3, cosTh = 0.8, corTh = 0.8){
#   if(is.na(cpid)) stop("Please set cpid!")
#   chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
#                                            xcms::chromPeakData(ndata)),
#                                      rownames = "cpid")
#   chromPeakTable_ms1 <- chromPeakTable %>% dplyr::filter(ms_level == 1)
#   chromPeakTable_ms2 <- chromPeakTable %>% dplyr::filter(ms_level == 2)
#   mchromPeakTable_ms1 <- chromPeakTable[chromPeakTable$cpid == cpid, ]
#   if(mchromPeakTable_ms1$ms_level != 1) stop("msLevel should be 1!")
#   peakWidth <- mchromPeakTable_ms1$rtmax - mchromPeakTable_ms1$rtmin
#   rtRange <- c(mchromPeakTable_ms1$rt - factor * peakWidth,
#                mchromPeakTable_ms1$rt + factor * peakWidth)
#   mchromPeakTable_ms2 <- chromPeakTable_ms2 %>%
#     dplyr::filter(rtmin >= rtRange[1] & rtmax <= rtRange[2]) %>%
#     dplyr::filter(mchromPeakTable_ms1$mz >= isolationWindowLowerMz & mchromPeakTable_ms1$mz <= isolationWindowUpperMz)
#   # remove ms2 feature mz >= precursorMz
#   mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
#     dplyr::filter(mzmax < mchromPeakTable_ms1$mzmin)
#   # remove noise
#   mchromPeakTable_ms2 <- mchromPeakTable_ms2 %>%
#     dplyr::filter(maxo > max(mchromPeakTable_ms2$maxo) * 0.05)
#   if(nrow(mchromPeakTable_ms2) == 0) return(NULL)
#   else{
#     chrDfList_ms2 <- lapply(mchromPeakTable_ms2$cpid, function(x) {
#       chrDf_x <- generate_chrDf(ndata = ndata, cpid = x, smooth = smooth, size = size, noise = noise)
#     })
#     names(chrDfList_ms2) <- mchromPeakTable_ms2$cpid
#   }
#   chrDf_ms1 <- generate_chrDf(ndata = ndata, cpid = mchromPeakTable_ms1$cpid,
#                               smooth = smooth, size = size)
#   chrDfList_ms1 <- list(chrDf_ms1)
#   names(chrDfList_ms1) <- mchromPeakTable_ms1$cpid
#   #browser()
#   if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
#   score <- sapply(1:length(chrDfList_ms2), function(i) {
#     vec1_rt <- chrDf_ms1$rt;vec2_rt <- chrDfList_ms2[[i]]$rt
#     vec1_int <- chrDf_ms1$intensity;vec2_int <- chrDfList_ms2[[i]]$intensity
#     tmp <- align_vectors(vec1 = vec1_rt, vec2 = vec2_rt)
#     vec1_rt <- tmp$aligned_vec1;vec2_rt <- tmp$aligned_vec2
#     vec1_rt[!is.na(vec1_rt)] <- vec1_int;vec2_rt[!is.na(vec2_rt)] <- vec2_int
#     vec1_int <- vec1_rt;vec2_int <- vec2_rt
#     cosine_similarity(vec1_int, vec2_int)
#   })
#   chrDfList_ms2 <- chrDfList_ms2[which(score > cosTh)]
#   if(length(chrDfList_ms2) == 0) return(list(ms1 = chrDfList_ms1, ms2 = NULL))
#   score <- sapply(1:length(chrDfList_ms2), function(i) {
#     calCor_shape(chrDf_ms1, chrDfList_ms2[[i]])
#   })
#   chrDfList_ms2 <- chrDfList_ms2[which(score > corTh)]
#   chrDfList <- list(ms1 = chrDfList_ms1, ms2 = chrDfList_ms2)
#   return(chrDfList)
# }
# align_vectors <- function(vec1, vec2) {
#   # 创建一个矩阵来存储两个向量之间的差值
#   diff_matrix <- outer(vec1, vec2, FUN = function(x, y) abs(x - y))
#
#   # 找到最小差值的位置
#   min_indices <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
#   best_i <- min_indices[1, 1]
#   best_j <- min_indices[1, 2]
#
#   # 计算对齐后的长度
#   max_length <- max(length(vec1) + best_j - 1, length(vec2) + best_i - 1)
#
#   # 对齐向量
#   aligned_vec1 <- c(rep(NA, best_j - 1), vec1, rep(NA, max_length - (length(vec1) + best_j - 1)))
#   aligned_vec2 <- c(rep(NA, best_i - 1), vec2, rep(NA, max_length - (length(vec2) + best_i - 1)))
#
#   list(aligned_vec1 = aligned_vec1, aligned_vec2 = aligned_vec2)
# }
# cosine_similarity <- function(vec1, vec2) {
#   # 计算点积
#   dot_product <- sum(vec1 * vec2, na.rm = TRUE)
#
#   # 计算范数
#   norm_vec1 <- sqrt(sum(vec1^2, na.rm = TRUE))
#   norm_vec2 <- sqrt(sum(vec2^2, na.rm = TRUE))
#
#   # 计算余弦相似度
#   cosine_sim <- dot_product / (norm_vec1 * norm_vec2)
#
#   return(cosine_sim)
# }
# fill_na_with_mean <- function(vec) {
#   # 获取向量的长度
#   n <- length(vec)
#
#   # 遍历向量中的每个元素
#   for (i in 1:n) {
#     # 如果当前元素是NA
#     if (is.na(vec[i])) {
#       # 找到前一个非NA值
#       forward <- NA
#       if (i > 1) {
#         for (j in (i-1):1) {
#           if (!is.na(vec[j])) {
#             forward <- vec[j]
#             break
#           }
#         }
#       }
#
#       # 找到后一个非NA值
#       backward <- NA
#       if (i < n) {
#         for (j in (i+1):n) {
#           if (!is.na(vec[j])) {
#             backward <- vec[j]
#             break
#           }
#         }
#       }
#
#       # 计算前后两个非NA值的平均值，并填充当前的NA
#       if (!is.na(forward) && !is.na(backward)) {
#         vec[i] <- (forward + backward) / 2
#       } else if (!is.na(forward)) {
#         vec[i] <- forward
#       } else if (!is.na(backward)) {
#         vec[i] <- backward
#       }
#     }
#   }
#
#   return(vec)
# }
# peak2spectra <- function(clusterPeaks){
#   #browser()
#   peak_ms1 <- clusterPeaks$ms1[[1]]
#   peak_ms1$mz <- fill_na_with_mean(peak_ms1$mz)
#   idx_ms1 <- which.max(peak_ms1$intensity)
#   rt_ms1 <- peak_ms1$rt[idx_ms1]
#   mz_ms1 <- peak_ms1$mz[idx_ms1]
#   int_ms1 <- peak_ms1$intensity[idx_ms1]
#   peaks_ms2 <- clusterPeaks$ms2
#   df_ms1 <- data.frame(mz = mz_ms1, intensity = int_ms1)
#   tmp <- lapply(peaks_ms2, function(x) {
#     x$mz <- fill_na_with_mean(x$mz)
#     idx_ms2 <- which.max(x$intensity)
#     mz_ms2 <- x$mz[idx_ms2]
#     int_ms2 <- x$intensity[idx_ms2]
#     return(data.frame(mz = mz_ms2, intensity = int_ms2))
#   })
#   spDf <- rbind(purrr::list_rbind(tmp), df_ms1) %>%
#     dplyr::arrange(mz)
#   spd <- dplyr::tibble(msLevel = 2L, rtime = rt_ms1)
#   spd$mz <- list(spDf$mz)
#   spd$intensity <- list(spDf$intensity)
#   sp <- Spectra::Spectra(spd)
#   return(sp)
# }
#
# Deconv4ndata <- function(ndata, smooth = TRUE, size = 3, factor = 0.5, noise = 0, cosTh = 0.8, corTh = 0.8, thread = 1){
#   chromPeakTable <- chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(ndata),
#                                                              xcms::chromPeakData(ndata)),
#                                                        rownames = "cpid")
#   chromPeakTable_ms1 <- chromPeakTable %>%
#     dplyr::filter(ms_level == 1)
#   chromPeakTable_ms2 <- chromPeakTable %>%
#     dplyr::filter(ms_level == 2)
#   #browser()
#   loop <- function(x){
#     clusterPeaks <- cluseter_peak(ndata = ndata, cpid = x, noise = noise,cosTh = cosTh, corTh = corTh)
#     if(is.null(clusterPeaks) | is.null(clusterPeaks$ms2) | length(clusterPeaks$ms2) == 0) return(NULL)
#     sp <- peak2spectra(clusterPeaks)
#     return(sp)
#   }
#
#   pb <- utils::txtProgressBar(max = nrow(chromPeakTable_ms1), style = 3)
#   #browser()
#   if(thread == 1){
#     spectraList <- lapply(1:nrow(chromPeakTable_ms1), function(m) {
#       utils::setTxtProgressBar(pb, m)
#       x <- chromPeakTable_ms1[m, ]$cpid
#       #print(m)
#       loop(x)
#     })
#   }else if(thread > 1){
#     cl <- snow::makeCluster(thread)
#     doSNOW::registerDoSNOW(cl)
#     opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
#                                                                  n))
#     envir <- environment(cluseter_peak)
#     parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
#     envir <- environment(smoothMean)
#     parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
#     envir <- environment(peak2spectra)
#     parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
#     envir <- environment(generate_chrDf)
#     parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
#     spectraList <- foreach::`%dopar%`(foreach::foreach(m = 1:nrow(chromPeakTable_ms1),
#                                                        .packages = c("dplyr"),
#                                                        .options.snow = opts),
#                                       {
#                                         x <- chromPeakTable_ms1[m, ]$cpid
#                                         loop(x)
#                                       })
#     snow::stopCluster(cl)
#     gc()
#   }else stop("thread is wrong!")
#   chromPeakTable_ms1$spectra <- spectraList
#   return(chromPeakTable_ms1)
# }
# #Deconv4ndata_new(ndata = data_filter[1], thread = 3)
