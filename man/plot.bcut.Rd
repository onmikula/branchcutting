% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.bcut.R
\name{plot.bcut}
\alias{plot.bcut}
\title{Graphical display of \code{branchcutting} results.}
\usage{
\method{plot}{bcut}(
  bc,
  which = "tree",
  edge.width = 2,
  show.tip.label = FALSE,
  tip.cex = 0.7,
  pt.cex = 2,
  show.legend = FALSE,
  ...
)
}
\arguments{
\item{bc}{an object of class "bcut".}

\item{which}{a type of plot to be created, options are "tree" (default), "species", "score" and "density",
unambiguous abbreviations are accepted.}

\item{edge.width}{an object of class "bcut".}

\item{show.tip.label}{logical, whether to show the tip labels (as in \code{plot.phylo}).}

\item{tip.cex}{factor scaling tip labels when \code{show.tip.label=TRUE}.}

\item{pt.cex}{factor scaling points when \code{which="species"}.}

\item{show.legend}{not implemented yet.}

\item{...}{Other graphical parameters.}
}
\description{
\code{plot} method for class "bcut".
}
