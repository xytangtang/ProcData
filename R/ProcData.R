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
#' @section Read sequences:
#' \itemize{
#'   \item \code{\link{csv2seqs}} reads action sequences from a csv file.
#' }
#' @section Sequence generators:
#' \itemize{
#'   \item \code{\link{seq_gen}} generates action sequences of an imaginery simulation-based item.
#'
#'   \item \code{\link{seq_gen2}} generates action sequences according to a given probability
#' transition matrix.
#' }
#' @section Feature extraction methods:
#' \itemize{
#'   \item \code{\link{seq2feature_mds}} extracts \code{K} features from action sequences by
#' multidimensional scaling.
#' 
#'   \item \code{\link{seq2feature_seq2seq}} extract features from action sequences by action
#' sequence autoencoder.
#' }
#' @section Sequence models:
#' \itemize{
#'   \item \code{\link{seqm}} fits a neural network model that relates action sequences 
#'     with a response variable.
#'     
#'   \item \code{\link{predict.seqm}} makes predictions from the models fitted by \code{seqm}.
#' }
#' @useDynLib ProcData
#' @importFrom Rcpp sourceCpp
#' @import keras
#' @import stats
#' @import utils
"_PACKAGE"
