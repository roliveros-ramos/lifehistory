% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{get_prior}
\alias{get_prior}
\title{Get priors for lifehistory parameters}
\usage{
get_prior(taxon, par = NULL)
}
\arguments{
\item{taxon}{The taxon you want to the priors.}

\item{par}{The parameters you want, possible options are "Linf", "k", "Winf",
"tmax", "tm", "M", "Lm" and "temperature".}
}
\value{
A list containing the mean and covariance for the given parameters.
}
\description{
Get priors for lifehistory parameters
}
