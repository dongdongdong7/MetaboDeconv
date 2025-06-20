% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromPeaks_similarity.R
\name{ccf_ms1PeakAms2Peaks}
\alias{ccf_ms1PeakAms2Peaks}
\title{ccf_ms1PeakAms2Peaks}
\usage{
ccf_ms1PeakAms2Peaks(ms1PeakDT, ms2PeaksDT)
}
\arguments{
\item{ms1PeakDT}{`data.table()`, a one row data.table with ms1 peak information
(mz, mzmin, mzmax, rt, rtmin, rtmax, maxo, intensity, rtime)}

\item{ms2PeaksDT}{`data.table()`, a data.table with ms2 peaks information corresponding ms1 peak}
}
\value{
`matrix()` with score, lag and maxo
}
\description{
Calculate ccf score matrix for ms1Peak and its corresponding ms2Peaks
}
\examples{
TODO: write examples
}
