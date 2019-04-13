#' Reading action sequences 
#' 
#' Reads a csv file and creates a list of action sequences.
#' 
#' This function first calls \code{read.csv} to read in data in the csv file and then
#' organize the action sequences stored in \code{seq_var} as a list of action sequences.
#' Each action sequence in the list is a vector of actions and is named by the concatenation
#' of the corresponsing values in \code{id_vars}.
#' 
#' @param file the name of the file from which the action sequences are to be read.
#' @param file_style the style that the action sequences are stored. \code{"oneline"}
#'   means each sequence occupies one line. Actions are separated by \code{seq_sep}.
#'   \code{"multiline"} means that each sequence occupies multiple lines. Each line
#'   represents an action in a sequence.
#' @param id_vars a string or a vector of strings giving the names of the variables
#'   identify different sequences. 
#' @param seq_var a string giving the name of the variable storing action sequences.
#' @param seq_sep the action separator characters. It is only used if \code{file_style="oneline"}.
#' @param ... further arguments to be passed to \code{read.csv}.
#' @return A list containing action sequences. Each element is an action sequence.
csv2seqs <- function(file, file_style, id_vars, seq_var, seq_sep = ",", ...) {
  if (!(file_style %in% c("multiline", "oneline")))
    stop("Invalid file style! Available options: multiline, and oneline.\n")
  csv_data <- read.csv(file = file, header = TRUE, stringsAsFactors = FALSE, ...)
  if (length(id_vars) == 1) {
    id_var <- id_vars
  } else {
    id_var <- paste(id_vars, collapse = "+")
    csv_data[, id_var] <- do.call(paste, c(csv_data[, id_vars], sep=""))
  }
  if (file_style == "multiline") {
    logfiles <- split(csv_data, csv_data[, id_var])
    seqs <- sapply(logfiles, function(x) x[, seq_var])
  } else if (file_style == "oneline") {
    logfiles <- csv_data[, seq_var]
    seqs <- strsplit(logfiles, split = seq_sep)
    names(seqs) <- csv_data[, id_var]
  }
  
  seqs
}