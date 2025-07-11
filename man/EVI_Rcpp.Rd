% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EVI.R
\name{EVI_Rcpp}
\alias{EVI_Rcpp}
\title{Compute the posterior expected Variation of Information (EVI) for a partition}
\usage{
EVI_Rcpp(cls, cls.draw, Ks, Ks.draw)
}
\arguments{
\item{cls}{A vector containing the cluster assignments of $n$ data points, representing the partition with respect to which the EVI is computed.}

\item{cls.draw}{A matrix of the MCMC samples of partitions of $n$ data points, of dimensions $S × n$, where $S$ is the number of MCMC samples.}

\item{Ks}{The number of clusters in \code{cls}}

\item{Ks.draw}{A vector of $S$, containing the number of clusters for each of the MCMC samples.}
}
\value{
The posterior expected VI for the partition \code{cls}, where the posterior is approximated by the MCMC samples \code{cls.draw}. It corresponds to the average VI distance between \code{cls} and each MCMC sample in \code{cls.draw}.
}
\description{
Compute the posterior expected Variation of Information (EVI) for a partition
}
