#' ProcData: A package for process data analysis
#'
#' General tools for exploratory process data analysis. Process data refers to
#' the data describing participants' problem solving processes in computer-based
#' assessments. It is often recorded in computer log files. This package
#' provides two action sequence generators and implements two automatic feature
#' extraction methods that compress the information stored in process data,
#' which often has a nonstandard format, into standard numerical vectors. This
#' package also provides recurrent neural network based models that relate
#' response processes with other binary or scale variables of interest. The
#' functions that involve training and evaluating neural networks are wrappers
#' of functions in keras.
#'
#' @section Sequence generator: two generators
#'
#' @section Feature extraction methods: two feature extraction methods
#'
#' @section Recurrent neural network based models: binary and scale
#' 
#' @useDynLib ProcData
#' @importFrom Rcpp sourceCpp
#' @import keras
#' @import stats
#' @import utils
"_PACKAGE"