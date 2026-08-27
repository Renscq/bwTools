\name{read_bwg}
\alias{read_bwg}
\title{Read BigWig, WIG, or bedGraph signal tracks}
\usage{
read_bwg(files, format = c("auto", "bigwig", "wig", "bedgraph"),
  sample_names = NULL, strand = "*", mode = c("memory", "lazy"))
}
\arguments{
\item{files}{One or more local signal files.}
\item{format}{Input format or \code{"auto"}.}
\item{sample_names}{Optional sample identifiers.}
\item{strand}{Strand labels associated with the input files.}
\item{mode}{\code{"memory"} or \code{"lazy"}; lazy access is currently available for BigWig only.}
}
\value{A BwgTrack-compatible object using 1-based closed in-memory coordinates.}
\description{Reads supported genomic signal formats. BigWig parsing is implemented with native R random-access binary I/O and zlib decompression.}
