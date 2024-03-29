% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KdModelList.R
\docType{class}
\name{KdModelList-class}
\alias{KdModelList-class}
\alias{KdModelList}
\title{KdModelList}
\usage{
KdModelList(..., description = NULL, makeUnique = FALSE)
}
\arguments{
\item{...}{Any number of \code{\link{KdModel}} objects or lists thereof.}

\item{description}{A description for the collection.}

\item{makeUnique}{Logical; whether to rename models if names are duplicated.}
}
\value{
A KdModelList
}
\description{
KdModelList
}
\examples{
data(SampleKdModel)
mods <- KdModelList(SampleKdModel, SampleKdModel, makeUnique = TRUE)
mods
}
