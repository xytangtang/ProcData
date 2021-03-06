% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtask.R
\name{segment_function}
\alias{segment_function}
\title{Segment an entropy sequence}
\usage{
segment_function(entropy_seq, lambda)
}
\arguments{
\item{entropy_seq}{a vector of entropies}

\item{lambda}{a number between 0 and 1.}
}
\value{
a vector of segment boundaries
}
\description{
\code{segment_function} segments the entropy sequence \code{entropy_seq} 
by identifying deep U-shaped curves in it.
}
