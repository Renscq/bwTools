\name{bw_track}
\alias{bw_track}
\alias{is_bwg_track}
\title{Construct and inspect BwgTrack-compatible signal objects}
\usage{
bw_track(samples, data = NULL, seqinfo = NULL, meta = list())

is_bwg_track(x)
}
\arguments{
\item{samples}{Sample metadata containing at least \code{sample_id}.}
\item{data}{Optional in-memory signal table.}
\item{seqinfo}{Optional chromosome metadata.}
\item{meta}{Optional metadata list.}
\item{x}{An R object.}
}
\value{\code{bw_track()} returns a BwgTrack-compatible S3 object. \code{is_bwg_track()} returns a logical scalar.}
\description{Creates the canonical genomic signal object used by bwTools and compatible with GeneTrackR signal consumers.}
