#' Class \code{"proc"} constructor
#' 
#' Create a \code{"proc"} object from given action sequences and timestamp sequences
#' 
#' An object of 
#' class \code{"proc"} is a list containing the following components:
#' \itemize{
#'   \item{action_seqs}{a list of action sequences.}
#'   \item{time_seqs}{a list of timestamp sequences.}
#' }
#' The names of the elements in \code{seqs$action_seqs} and \code{seqs$time_seqs} are 
#' process identifiers.
#'  
#' @param action_seqs a list of action sequences.
#' @param time_seqs a list of timestamp sequences.
#' @param ids ids identifiers of response processes.
#' @export
proc <- function(action_seqs, time_seqs, ids = NULL) {
  l_action <- length(action_seqs)
  ls_action <- sapply(action_seqs, length)
  if (!is.null(time_seqs)) {
    l_time <- length(time_seqs)
    if (l_action != l_time)
      stop("Number of timestamp sequences does not match the number of action sequences!\n")
    ls_time <- sapply(time_seqs, length)
    if (any(ls_time != ls_action))
      stop("Lengths of action sequences and lengths of timestamp sequences do not match!\n")
    if (any(unlist(sapply(time_seqs, diff)) < 0))
      stop("Timestamp sequences are not non-decreasing!\n")
  }  
  if (!is.null(ids)) {
    if (length(ids) != l_action) stop("Number of provided IDs does not match number of sequences!\n")
    names(action_seqs) <- ids
    if (!is.null(time_seqs)) names(time_seqs) <- ids
  }
  
  res <- list(action_seqs=action_seqs, time_seqs=time_seqs)
  
  class(res) <- "proc"
  
  res
}

#' Summary method for class \code{"proc"}
#' 
#' The summary of a "proc" object combines the summary of the action sequences and the summary
#' of the timestamp sequences.
#' 
#' @param object an object of class \code{"\link{proc}"}.
#' @param ... not used.
#' @return a list. Its components are the components returned by \link{action_seqs_summary} and
#'   \link{time_seqs_summary}.
#' @seealso \link{action_seqs_summary} and \link{time_seqs_summary}
#' @export
summary.proc <- function(object, ...) {
  res_action <- action_seqs_summary(object$action_seqs)
  
  res_time <- time_seqs_summary(object$time_seqs)
  
  c(res_action, res_time)
}

