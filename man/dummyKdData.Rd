% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{dummyKdData}
\alias{dummyKdData}
\title{Create dummy log_kd per 12-mer data}
\usage{
dummyKdData(mod = NULL)
}
\arguments{
\item{mod}{Optional model from which to create the dummy data}
}
\value{
A data.frame with 12-mers and log_kds
}
\description{
Create dummy log_kd per 12-mer data
}
\examples{
kd <- dummyKdData()
}
