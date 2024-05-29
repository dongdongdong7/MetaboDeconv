# MetaboDeconv

## Important variable

### featureTable

一个包含有全部feature的table，包含这个feature所处的SWATH窗口。

```R
> featureTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(swath_data),
+                                        xcms::chromPeakData(swath_data)),
+                                  rownames = "feature_id")
> featureTable
# A tibble: 370 × 18
   feature_id    mz mzmin mzmax    rt rtmin rtmax   into   intb  maxo    sn sample ms_level is_filled
   <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>    <int> <lgl>    
 1 CP01        220.  220.  220.  239.  237.  242.  1281.  1277.  312.    89      1        1 FALSE    
 2 CP02        219.  219.  219.  240.  236.  244. 11080. 11010. 2683.    55      1        1 FALSE    
 3 CP03        153.  153.  153.  330.  326.  335.  2056.  2047.  572.   323      1        1 FALSE    
 4 CP04        235.  235.  235.  330.  326.  335.  2706.  2697.  755.   160      1        1 FALSE    
 5 CP05        299.  299.  299.  346.  341.  352.  3988.  3978.  714.   698      1        1 FALSE    
 6 CP06        298.  298.  298.  346.  341.  351. 22961. 22888. 5342.   109      1        1 FALSE    
 7 CP07        306.  306.  306.  357.  353.  362.  1307.  1300.  263.    93      1        1 FALSE    
 8 CP08        228.  228.  228.  357.  353.  363.  3995.  3986.  903.   902      1        1 FALSE    
 9 CP09        589.  589.  589.  357.  353.  363.  2942.  2932.  583.   448      1        1 FALSE    
10 CP10        301.  301.  301.  365.  363.  370.   920.   913.  301.   115      1        1 FALSE    
# ℹ 360 more rows
# ℹ 4 more variables: isolationWindow <dbl>, isolationWindowTargetMZ <dbl>,
#   isolationWindowLowerMz <dbl>, isolationWindowUpperMz <dbl>
# ℹ Use `print(n = ...)` to see more rows
```

一级和二级```featureTable```

```R
featureTable_ms1 <- featureTable %>%
  dplyr::filter(ms_level == 1)
featureTable_ms2 <- featureTable %>%
  dplyr::filter(ms_level == 2)
```

```mfeatureTable_ms2```是第m个一级feature相关的二级特征，同理应该还有一个变量```mfeatureTable_ms1```。

```R
fenamiphos_mz <- 304.113077
feature_fenamiphos <- featureTable_ms1 %>%
  dplyr::filter(dplyr::near(mz, fenamiphos_mz, tol = 0.001))
fenamiphos_rt <- feature_fenamiphos$rt
peakWidth <- feature_fenamiphos$rtmax - feature_fenamiphos$rtmin
rtRange <- c(fenamiphos_rt - 2.5 * peakWidth, fenamiphos_rt + 2.5 * peakWidth)
mfeatureTable_ms2 <- featureTable_ms2 %>%
  dplyr::filter(rt >= rtRange[1] & rt <= rtRange[2]) %>%
  dplyr::filter(fenamiphos_mz >= isolationWindowLowerMz & fenamiphos_mz <= isolationWindowUpperMz)
```

