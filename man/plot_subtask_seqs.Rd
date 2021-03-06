% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtask.R
\name{plot_subtask_seqs}
\alias{plot_subtask_seqs}
\title{Plot Subtask Analysis Results for Entire Dataset}
\usage{
plot_subtask_seqs(subtask_seqs, subtasks, max_len = 5,
  col.subtask = 1:length(subtasks), plot_legend = TRUE,
  legend_pos = "topright", ...)
}
\arguments{
\item{subtask_seqs}{a list of subtask sequences}

\item{subtasks}{a vector of all subtasks}

\item{max_len}{maximum length of plotted subtasks}

\item{col.subtask}{a vector of colors for subtasks}

\item{plot_legend}{a logical value. If \code{TRUE} (default), plot the legend of subtasks.}

\item{legend_pos}{a character string or the coordinates to be used to position the legend.}

\item{...}{other arguments passed to \code{\link{legend}}}
}
\value{
this function does not return values
}
\description{
Plot Subtask Analysis Results for Entire Dataset
}
\seealso{
\code{\link{plot_subtask_seq}} for ploting results for one sequence. 
  \code{\link{plot.subtask}} for ploting an object of class \code{"subtask"}
}
