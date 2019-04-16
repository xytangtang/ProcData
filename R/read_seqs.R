#' Reading action sequences 
#' 
#' Reads a csv file and creates a list of action sequences.
#' 
#' This function first calls \code{read.csv} to read in data in the csv file and then
#' organize the action sequences stored in \code{seq_var} as a list of action sequences.
#' Each action sequence in the list is a vector of actions and is named by the concatenation
#' of the corresponsing values in \code{id_vars}.
#' 
#' @param file the name of the csv file from which the action sequences are to be read.
#' @param style the style that the action sequences are stored. \code{"single"}
#'   means each sequence occupies a single line. Actions are separated by \code{seq_sep}.
#'   \code{"multiple"} means that each sequence occupies multiple lines, each action in
#'   a line.
#' @param id_vars a string or a vector of strings giving the names of the variables
#'   identify different sequences. 
#' @param seq_var a string giving the name of the variable storing action sequences.
#' @param seq_sep the action separator characters. It is only used if \code{style="single"}.
#' @param ... further arguments to be passed to \code{read.csv}.
#' @return A list containing action sequences. Each element is an action sequence.
read.seqs <- function(file, style, id_vars="ID", seq_var="Event", seq_sep = ",", ...) {
  if (!(style %in% c("multiple", "single")))
    stop("Invalid file style! Available options: multiple, and single.\n")
  csv_data <- read.csv(file = file, header = TRUE, stringsAsFactors = FALSE, ...)
  
  if (any(!(id_vars %in% names(csv_data)))) stop("Variables in 'id_vars' do not exist!\n")
  if (length(seq_var) > 1) {
    warning("More than one variable are given in 'seq_var'. Only the first one will be used\n")
    seq_var <- seq_var[1]
  }
  if (!(seq_var %in% names(csv_data))) stop("Variable in 'seq_var' does not exist!\n")
  
  if (length(id_vars) == 1) {
    id_var <- id_vars
  } else {
    id_var <- paste(id_vars, collapse = "+")
    csv_data[, id_var] <- do.call(paste, c(csv_data[, id_vars], sep=""))
  }
  if (style == "multiple") {
    logfiles <- split(csv_data, csv_data[, id_var])
    seqs <- sapply(logfiles, function(x) x[, seq_var])
  } else if (style == "single") {
    logfiles <- csv_data[, seq_var]
    seqs <- strsplit(logfiles, split = seq_sep)
    if (length(unique(csv_data[, id_var])) < length(csv_data[, id_var])) warning("Repeated IDs detected!\n")
    names(seqs) <- csv_data[, id_var]
  }
  
  seqs
}