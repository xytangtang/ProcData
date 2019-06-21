#' Reading response processes from csv files
#' 
#' Reads a csv file and creates response process data.
#' 
#' \code{read.seqs} calls \code{read.csv} to read process data stored in a csv file into \code{R}.
#' The csv file to be read should at least include an identifier of distinct response processes, 
#' and action sequences. It can also include timestamp sequences.
#' 
#' The response processes (action sequences and timestamp sequences) stored in csv files can 
#' be in one of the two styles, \code{"single"} and \code{"multiple"}. In \code{"single"} style,
#' each response process occupies a single line. Actions and timestamps at different steps 
#' are separated by \code{step_sep}. In \code{"multiple"} stype, each response process occupies 
#' multiple lines with each step taking up one line.
#' 
#' @param file the name of the csv file from which the response processes are to be read.
#' @param style the style that the response processes and are stored. See 'Details'.
#' @param id_var a string giving the name of the variable storing the process identifier. 
#' @param action_var a string giving the name of the variable storing action sequences.
#' @param time_var a string giving the name of the variable storing timestamp sequences.
#' @param step_sep the step separator characters. It is only used if \code{style="single"}.
#' @param ... further arguments to be passed to \code{read.csv}.
#' @return \code{read.seqs} returns an object of class \code{"proc"}. An object of 
#' class \code{"proc"} is a list containing the following components:
#'   \item{action_seqs}{a list of action sequences.}
#'   \item{time_seqs}{a list of timestamp sequences.}
#'  The names of elements in \code{action_seqs} and \code{time_seqs} are process identifier given 
#'  by \code{id_var}.
#' @export
read.seqs <- function(file, style, id_var=NULL, action_var=NULL, time_var=NULL, step_sep = ",", ...) {
  if (!(style %in% c("multiple", "single")))
    stop("Invalid file style! Available options: multiple, and single.\n")
  
  csv_data <- read.csv(file = file, header = TRUE, stringsAsFactors = FALSE, ...)
  
  
  if (length(id_var) > 1) {
    warning("More than one variable is given in 'id_var'. Only the first one will be used\n")
    id_var <- id_var[1]
  }
  if (is.null(id_var)) stop("id_var should be provided!\n")
  if (!(id_var %in% names(csv_data))) stop("Variables in 'id_var' do not exist!\n")
  
  if (length(action_var) > 1) {
    warning("More than one variable is given in 'action_var'. Only the first one will be used\n")
    action_var <- action_var[1]
  }
  if (is.null(action_var)) stop("action_var should be provided!\n")
  if (!(action_var %in% names(csv_data))) stop("Variable in 'action_var' does not exist!\n")
  
  if (length(time_var) > 1) {
    warning("More than one variable is given in 'time_var'. Only the first one will be used\n")
    time_var <- time_var[1]
  }
  if (!(is.null(time_var)) && !(time_var %in% names(csv_data))) {
    warning("Ignoring time_var as the variable does not exist!\n")
    time_var <- NULL
  }
  
  if (style == "multiple") {
    logfiles <- split(csv_data, csv_data[, id_var])
    actions <- sapply(logfiles, function(x) x[, action_var])
    if (is.null(time_var)) times <- NULL
    else times <- sapply(logfiles, function(x) x[, time_var])
  } else if (style == "single") {
    if (length(unique(csv_data[, id_var])) < length(csv_data[, id_var])) 
      warning("Repeated IDs detected!\n")
    actions <- strsplit(csv_data[, action_var], split = step_sep, fixed=TRUE)
    names(actions) <- csv_data[, id_var]
    
    if (is.null(time_var)) times <- NULL
    else {
      times <- strsplit(csv_data[, time_var], split = step_sep, fixed=TRUE)
      names(times) <- csv_data[, id_var]
    }
  }
  
  res <- list(action_seqs = actions, time_seqs = times)
  class(res) <- "proc"
  
  res
}