\name{bw_detect_format}
\alias{bw_detect_format}
\title{Detect a supported genomic signal format}
\usage{bw_detect_format(file)}
\arguments{\item{file}{A local signal file path.}}
\value{One of \code{"bigwig"}, \code{"wig"}, or \code{"bedgraph"}.}
\description{Uses file suffixes, BigWig magic bytes, and text headers to detect supported formats.}
