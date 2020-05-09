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
#' @return an object of class \code{"\link{proc}"} containing the provided action and
#'  timestamp sequences.
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
  
  if (is.null(object$time_seqs)) res_time <- NULL
  else
    res_time <- time_seqs_summary(object$time_seqs)
  
  c(res_action, res_time)
}

#' Print method for class \code{"proc"}
#' 
#' @param x an object of class \code{"\link{proc}"}
#' @param n number of processes to be printed.
#' @param index indice of processes to be printed.
#' @param quote logical, indicating whether or not strings should be printed with surrounding quotes.
#' @param ... not used.
#' @return \code{print.proc} invisibly returns the \code{"proc"} object it prints.
#' @export
print.proc <- function(x, n=5, index=NULL, quote=FALSE, ...) {
  n_total <- length(x$action_seqs)
  cat("'proc' object of ", n_total, " processes\n")
  cat("\n")
  if (is.null(n) & is.null(index)) n <- min(n_total, 5)
  else {
    if (!is.null(n)) n <- min(n_total, n)
    if (!is.null(index)) index <- index[index <= n_total]
    if (!is.null(n) & !is.null(index)) index <- index[1:min(n, length(index))]
  }
  if (is.null(index)) {
    cat("First ", n, " processes:\n")
    index <- 1:n
  }
  
  cat("\n")
  seq_names <- names(x$action_seqs)
  for (i in index) {
    l <- length(x$action_seqs[[i]])
    po <- data.frame(Event = x$action_seqs[[i]], stringsAsFactors = FALSE)
    if (!is.null(x$time_seqs)) po$Time=x$time_seqs[[i]]
    rownames(po) <- paste("Step", 1:l)
    cat(seq_names[i], "\n")
    print(t(po), quote=quote, ...)
    cat("\n")
  }
  invisible(x)
}