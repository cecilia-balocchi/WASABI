% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{ggscatter2d}
\alias{ggscatter2d}
\title{2D cluster scatterplot}
\usage{
ggscatter2d(cls, Y)
}
\arguments{
\item{cls}{Numeric or factor vector of cluster assignments for each row.}

\item{Y}{Numeric matrix or data.frame with at least two columns.}
}
\value{
A ggplot2 plot object.
}
\description{
Plots a scatter plot of the first two columns of feature data, colored by clusters.
}
\examples{
\dontrun{
ggscatter2d(output_wvi$particles[1,], Y)
}
}
