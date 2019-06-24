#' Class \code{proc} constructor
#' 
#' Create a \code{"proc"} object from given action sequences and timestamp sequences
#' 
#' An object of 
#' class \code{"proc"} is a list containing the following components:
#'   \item{action_seqs}{a list of action sequences.}
#'   \item{time_seqs}{a list of timestamp sequences.}
#'  The names of elements in \code{action_seqs} and \code{time_seqs} are process identifiers.
#'  
#' @param action_seqs a list of action sequences
#' @param time_seqs a list of timestamp sequences
#' @param ids ids identifiers of 
#' @export
proc <- function(action_seqs, time_seqs, ids = NULL) {
  
  if (!is.null(time_seqs)) {
    l_action <- length(action_seqs)
    l_time <- length(time_seqs)
    if (l_action != l_time)
      stop("Number of timestamp sequences does not match the number of action sequences!\n")
    ls_action <- sapply(action_seqs)
    ls_time <- sapply(time_seqs)
    if (any(action_seqs != time_seqs))
      stop("Lengths of action sequences and lengths of timestamp sequences do not match!\n")
    if (any(unlist(sapply(time_seqs, diff))) < 0)
      stop("Timestamp sequences are not non-decreasing!\n")
  }  
  if (!is.null(ids)) {
    if (length(ids) != l_action) stop("Number of provided IDs does not match number of sequences!\n")
    names(action_seqs) <- ids
    names(time_seqs) <- ids
  }
  
  res <- list(action_seqs=action_seqs, time_seqs=time_seqs)
  
  class(res) <- "proc"
  
  res
}

#' Summary method for class \code{"proc"}
#' 
#' @param seqs an object of class \code{\link{"proc"}}
summary.proc <- function(seqs) {
  res_action <- action_seqs_summary(seqs$action_seqs)
  
  res_time <- time_seqs_summary(seqs$time_seqs)
  
  c(res_action, res_time)
}

