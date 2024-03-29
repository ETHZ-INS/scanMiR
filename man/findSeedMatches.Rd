% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanning.R
\name{findSeedMatches}
\alias{findSeedMatches}
\title{Predicting and characterizing miRNA binding sites}
\usage{
findSeedMatches(
  seqs,
  seeds,
  shadow = 0L,
  onlyCanonical = FALSE,
  maxLogKd = c(-1, -1.5),
  keepMatchSeq = FALSE,
  minDist = 7L,
  p3.extra = FALSE,
  p3.params = list(maxMirLoop = 7L, maxTargetLoop = 9L, maxLoopDiff = 4L, mismatch =
    TRUE, GUwob = TRUE),
  agg.params = .defaultAggParams(),
  ret = c("GRanges", "data.frame", "aggregated"),
  BP = NULL,
  verbose = NULL,
  n_seeds = NULL,
  useTmpFiles = FALSE,
  keepTmpFiles = FALSE
)
}
\arguments{
\item{seqs}{A character vector or `DNAStringSet` of DNA sequences in which to
look.}

\item{seeds}{A character vector of 7-nt seeds to look for. If RNA, will be
reversed and complemented before matching. If DNA, they are assumed to be
the target sequence to look for. Alternatively, a list of objects of class
`KdModel` or an object of class `KdModelList` can be given.}

\item{shadow}{Integer giving the shadow, i.e. the number of nucleotides
hidden at the beginning of the sequence (default 0).}

\item{onlyCanonical}{Logical; whether to restrict the search only to
canonical binding sites.}

\item{maxLogKd}{Maximum log_kd value to keep. This has a major impact on the
number of sites returned, and hence on the memory requirements. Set to Inf
to disable (_not_ recommended when running large scans!).}

\item{keepMatchSeq}{Logical; whether to keep the sequence (including flanking
dinucleotides) for each seed match (default FALSE).}

\item{minDist}{Integer specifying the minimum distance between matches of the
same miRNA (default 7). Closer matches will be reduced to the
highest-affinity. To disable the removal of overlapping features, use
`minDist=-Inf`.}

\item{p3.extra}{Logical; whether to keep extra information about 3'
alignment. Disable (default) this when running large scans, otherwise you
might hit your system's memory limits.}

\item{p3.params}{Named list of parameters for 3' alignment with slots
`maxMirLoop` (integer, default = 7), `maxTargetLoop` (integer, default = 9),
`maxLoopDiff` (integer, default = 4), `mismatch`
(logical, default = TRUE) and `GUwob` (logical, default = TRUE).}

\item{agg.params}{A named list with slots `a`, `b`, `c`, `p3`, `coef_utr`,
`coef_orf` and `keepSiteInfo` indicating the parameters for the aggregation.
Ignored if `ret!="aggregated"`. For further details see documentation of
`aggregateMatches`.}

\item{ret}{The type of data to return, either "GRanges" (default),
"data.frame", or "aggregated" (aggregates affinities/sites for each
seed-transcript pair).}

\item{BP}{Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to
enable multithreading.}

\item{verbose}{Logical; whether to print additional progress messages
(default on if not multithreading)}

\item{n_seeds}{Integer; the number of seeds that are processed in parallel to
avoid memory issues.}

\item{useTmpFiles}{Logical; whether to write results for single miRNAs in
temporary files (ignored when scanning for a single seed). Alternatively,
`useTmpFiles` can be a character vector of length 1 indicating the path to
the directory in which to write temporary files.}

\item{keepTmpFiles}{Logical; whether to keep the temporary files at the end
of the process; ignored if `useTmpFiles=FALSE`. Temporary files are removed
only upon successful completion of the function, meaning that they will not
be deleted in case of errors.}
}
\value{
A GRanges of all matches. If `seeds` is a `KdModel` or `KdModelList`,
the `log_kd` column will report the ln(Kd) multiplied by 1000, rounded and
saved as an integer. If `ret!="GRanges`, returns a data.frame.
}
\description{
`findSeedMatches` takes a set of sequences and a set of miRNAs (given either
as target seeds, mature miRNA sequences, or a \code{\link{KdModelList}}).
}
\examples{
# we create mock RNA sequences and seeds:
seqs <- getRandomSeq(n=10)
seeds <- c("AAACCAC", "AAACCUU")
findSeedMatches(seqs, seeds)
}
