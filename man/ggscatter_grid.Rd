% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{ggscatter_grid}
\alias{ggscatter_grid}
\title{Grid of scatterplots: clusters in multiple particles (1D)}
\usage{
ggscatter_grid(output_wvi, Y)
}
\arguments{
\item{output_wvi}{WASABI output; must contain \code{particles} and \code{part.weights}.}

\item{Y}{Numeric vector of data.}
}
\value{
A ggplot2 plot object.
}
\description{
Make a faceted scatter plot: each facet shows cluster assignments for one particle.
}
\examples{
\dontrun{
ggscatter_grid(output_wvi, Y)
}
}
