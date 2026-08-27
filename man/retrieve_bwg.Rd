\name{retrieve_bwg}
\alias{retrieve_bwg}
\title{Retrieve a genomic signal interval}
\usage{
retrieve_bwg(object, chrom, start, end, samples = NULL, strand = "ignore")
}
\arguments{
\item{object}{A BwgTrack-compatible object.}
\item{chrom}{Chromosome name.}
\item{start}{Start coordinate in 1-based closed coordinates.}
\item{end}{End coordinate in 1-based closed coordinates.}
\item{samples}{Optional sample identifiers.}
\item{strand}{Optional strand selector.}
}
\value{A \code{data.table} containing matching intervals.}
\description{Retrieves signal from an in-memory object or performs indexed random access for lazy BigWig objects.}
