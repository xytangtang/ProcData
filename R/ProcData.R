#' @useDynLib ProcData
#' @importFrom Rcpp sourceCpp
#' @import keras
#' @import stats
#' @import utils
NULL

#' ProcData: A package for process data analysis
#'
#' General tools for exploratory process data analysis. Process data refers to
#' the data describing participants' problem solving processes in computer-based
#' assessments. It is often recorded in computer log files. This package a process 
#' dataset and functions for reading processes from a csv file, process manipulation,
#' action sequence generators. It also implements two automatic feature
#' extraction methods that compress the information stored in process data,
#' which often has a nonstandard format, into standard numerical vectors. This
#' package also provides recurrent neural network based models that relate
#' response processes with other binary or scale variables of interest. The
#' functions that involve training and evaluating neural networks are based on
#' functions in keras.
#'
#' @section Data structure:
#' \code{ProcData} organizes response processes as an object of class \code{\link{proc}}.
#' Some functions are provided for summarizing and manipulating \code{proc} objects.
#' \itemize{
#'   \item \code{\link{summary.proc}} calculates summary statistics of a \code{proc} object.
#'   \item \code{\link{remove_action}} removes actions and the corresponding timestamps
#'   \item \code{\link{replace_action}} replaces an action by another action
#'   \item \code{\link{combine_actions}} combines consecutive action into one action.
#' }
#' @section Read sequences:
#' \itemize{
#'   \item \code{\link{read.seqs}} reads response processes from a csv file.
#' }
#' @section Sequence generators:
#' \itemize{
#'   \item \code{\link{seq_gen}} generates action sequences of an imaginery simulation-based item.
#'
#'   \item \code{\link{seq_gen2}} generates action sequences according to a given probability
#' transition matrix.
#'
#'   \item \code{\link{seq_gen3}} generates action sequences according to a recurrent neural network.
#' }
#' @section Feature extraction methods:
#' \itemize{
#'   \item \code{\link{seq2feature_mds}} extracts features from response processes by
#'     multidimensional scaling.
#' 
#'   \item \code{\link{seq2feature_seq2seq}} extracts features from response processes by
#'     autoencoder.
#' }
#' @section Sequence models:
#' \itemize{
#'   \item \code{\link{seqm}} fits a neural network model that relates response processes 
#'     with a response variable.
#'     
#'   \item \code{\link{predict.seqm}} makes predictions from the models fitted by \code{seqm}.
#' }
"_PACKAGE"
