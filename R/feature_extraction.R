#' Feature extraction via multidimensional scaling
#'
#' \code{seq2feature_mds} extracts \code{K} features from response processes by
#' multidimensional scaling.
#' 
#' Since the classical MDS has a computational complexity of order \eqn{n^3} where 
#' \eqn{n} is the number of response processes, it is computational expensive to 
#' perform classical MDS when a large number of response processes is considered. 
#' In addition, storing an \eqn{n \times n} dissimilarity matrix when \eqn{n} is large
#' require a large amount of memory. In \code{seq2feature_mds}, the algorithm proposed
#' in Paradis (2018) is implemented to obtain MDS for large datasets. \code{method} 
#' specifies the algorithm to be used for obtaining MDS features. If \code{method = "small"},
#' classical MDS is used by calling \code{\link{cmdscale}}. If \code{method = "large"},
#' the algorithm for large datasets will be used. If \code{method = "auto"} (default), 
#' \code{seq2feature_mds} selects the algorithm automatically based on the sample size.
#' 
#' \code{dist_type} specifies the dissimilarity to be used for measuring the discrepancy
#' between two response processes. If \code{dist_type = "oss_action"}, the order-based 
#' sequence similarity (oss) proposed in Gomez-Alonso and Valls (2008) is used 
#' for action sequences. If \code{dist_type = "oss_both"}, both action sequences and
#' timestamp sequences are used to compute a time-weighted oss. 
#'   
#' The number of
#' features to be extracted \code{K} can be selected by cross-validation using
#' \code{\link{chooseK_mds}}.
#'
#' @family feature extraction methods
#' @param seqs a \code{"\link{proc}"} object or a square matrix. If a squared matrix is 
#'   provided, it is treated as the dissimilary matrix of a group of response processes.
#' @param K the number of features to be extracted.
#' @param method a character string specifies the algorithm used for performing MDS. See
#'   'Details'.
#' @param dist_type a character string specifies the dissimilarity measure for two
#'   response processes. See 'Details'.
#' @param pca logical. If \code{TRUE} (default), the principal components of the
#'   extracted features are returned.
#' @param subset_size,n_cand two parameters used in the large data algorithm. See 'Details'
#'   and \code{\link{seq2feature_mds_large}}.
#' @param subset_method a character string specifying the method for choosing the subset 
#'   in the large data algorithm. See 'Details' and \code{\link{seq2feature_mds_large}}.
#' @param return_dist logical. If \code{TRUE}, the dissimilarity matrix will be
#'   returned. Default is \code{FALSE}.
#' @param seed random seed.
#' @param L_set length of ngrams considered
#' @return \code{seq2feature_mds} returns a list containing 
#'   \item{theta}{a numeric matrix giving the \code{K} extracted features or principal
#'     features. Each column is a feature.} 
#'   \item{dist_mat}{the dissimilary matrix. This element exists only if 
#'     \code{return_dist=TRUE}.}
#' @seealso \code{\link{chooseK_mds}} for choosing \code{K}.
#' @references Gomez-Alonso, C. and Valls, A. (2008). A similarity measure for sequences of
#'   categorical data based on the ordering of common elements. In V. Torra & Y. Narukawa (Eds.) 
#'   \emph{Modeling Decisions for Artificial Intelligence}, (pp. 134-145). Springer Berlin Heidelberg.
#' @references Paradis, E. (2018). Multidimensional scaling with very large datasets. \emph{
#'   Journal of Computational and Graphical Statistics}, 27(4), 935-939.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' theta <- seq2feature_mds(seqs, 5)$theta
#' @export
seq2feature_mds <- function(seqs = NULL, K = 2, method = "auto", dist_type = "oss_action", 
                            pca = TRUE, subset_size = 100, subset_method = "random", n_cand = 10, 
                            return_dist = FALSE, seed = 12345, L_set = 1:3) {
  set.seed(seed)
  dist_ready <- FALSE
  if (is.null(seqs)) 
    stop("Either response processes or their dissimilarity matrix should be provided!\n")
  if (is.matrix(seqs)) {
    if (nrow(seqs) != ncol(seqs)) stop("Provided matrix is not square!\n")
    dist_mat <- seqs
    n <- nrow(dist_mat)
    dist_ready <- TRUE
  } else if (class(seqs) == "proc") {
    n <- length(seqs$action_seqs) 
  } else {
    stop("seqs should be a 'proc' object or a square matrix!\n")
  }
  
  if (n < K) stop("Not enough processes for extracting K features!\n")
  
  if (method == "auto") {
    if (n > 5000) method <- "large"
    else method <- "small"
  }
  if (method == "small") {
    if (n > 5000) warning("Using the small dataset method for large datasets!\n")
    if (!dist_ready) { # calculate dist matrix
      dist_mat <- matrix(0, n, n)
      if (is.null(seqs$time_seqs) & dist_type == "oss_both") {
        warning("Timestamp sequences are not available. 
                change dist_type from 'oss_both' to 'oss_action'!\n")
        dist_type <- "oss_action"
      }
      if (dist_type == "oss_action") {
        dist_mat <- calculate_dist_cpp(seqs$action_seqs)
      } else if (dist_type == "oss_both") {
        dist_mat <- calculate_tdist_cpp(seqs$action_seqs, seqs$time_seqs)
      } else if (dist_type == "ngram") {
        dist_mat <- calculate_ngram_dist_cpp(seqs$action_seqs, L_set)
      } else stop("Invalid dissimilarity method!\n")
    }
    
    theta <- cmdscale(dist_mat, K)
  } else { # method == "large"
    if (subset_size <= K) stop("Subset size is too small to extract K features!\n")
    if (subset_size > n) stop("Subset size is larger than the sample size!\n")
    theta <- seq2feature_mds_large(seqs = seqs, K = K, dist_type = dist_type, 
                          subset_size = subset_size, subset_method = subset_method, 
                          n_cand = n_cand, pca = pca, seed = seed, L_set = L_set)
  }
  
  if (return_dist & dist_ready) 
    res <- list(theta=theta, dist_mat=dist_mat)
  else res <- list(theta=theta)
  
  res
}




#' Choose the number of multidimensional scaling features
#' 
#' \code{chooseK_mds} choose the number of multidimensional scaling features
#'   to be extracted by cross-validation.
#'
#' @param K_cand the candidates of the number of features.
#' @param n_fold the number of folds for cross-validation.
#' @param max_epoch the maximum number of epochs for stochastic gradient
#'   descent.
#' @param step_size the step size of stochastic gradient descent.
#' @param tot the accuracy tolerance for determining convergence.
#' @inheritParams seq2feature_mds
#' @return \code{chooseK_mds} returns a list containing
#'   \item{K}{the value in \code{K_cand} producing the smallest cross-validation loss.}
#'   \item{K_cand}{the candidates of the number of features.}
#'   \item{cv_loss}{the cross-validation loss for each candidate in \code{K_cand}.}
#'   \item{dist_mat}{the dissimilary matrix. This element exists only if \code{return_dist=TRUE}.}
#' @seealso \code{\link{seq2feature_mds}} for feature extraction after choosing
#'   the number of features.
#' @references Gomez-Alonso, C. and Valls, A. (2008). A similarity measure for sequences of
#'   categorical data based on the ordering of common elements. In V. Torra & Y. Narukawa (Eds.) 
#'   \emph{Modeling Decisions for Artificial Intelligence}, (pp. 134-145). Springer Berlin Heidelberg.
#' @examples 
#' n <- 50
#' seqs <- seq_gen(n)
#' K_res <- chooseK_mds(seqs, 5:10, return_dist=TRUE)
#' theta <- seq2feature_mds(K_res$dist_mat, K_res$K)$theta
#'
#' @export 
chooseK_mds <- function(seqs=NULL, K_cand, dist_type="oss_action", n_fold=5, 
                        max_epoch=100, step_size=0.01, tot=1e-6, return_dist=FALSE,
                        seed = 12345, L_set = 1:3) {
  set.seed(seed)
  if (is.null(seqs)) 
    stop("Either response processes or their dissimilarity matrix should be provided!\n")
  if (is.matrix(seqs)) {
    if (nrow(seqs) != ncol(seqs)) stop("Provided matrix is not square!\n")
    dist_mat <- seqs
    n <- nrow(dist_mat)
    if (n > 5000) 
      warning("MDS cross-validation can take a long time due to large sample size!\n")
  } else if (class(seqs) == "proc") {
    n <- length(seqs$action_seqs)
    if (n > 5000)
      warning("MDS cross-validation can take a long time due to large sample size!\n")
    dist_mat <- matrix(0, n, n)
    if (is.null(seqs$time_seqs) & dist_type == "oss_both") {
      warning("Timestamp sequences are not available. 
              change dist_type from 'oss_both' to 'oss_action'!\n")
      dist_type <- "oss_action"
    }
    if (dist_type == "oss_action") {
      dist_mat <- calculate_dist_cpp(seqs$action_seqs)
    } else if (dist_type == "oss_both") {
      dist_mat <- calculate_tdist_cpp(seqs$action_seqs, seqs$time_seqs)
    } else if (dist_type == "ngram") {
      dist_mat <- calculate_ngram_dist_cpp(seqs$action_seqs, L_set)
    } else stop("Invalid dissimilarity method!\n")
  } else {
    stop("seqs should be a 'proc' object or a square matrix!\n")
  }
  
  n_K <- length(K_cand)
  n_pairs <- n * (n - 1) / 2
  all_pairs <- t(combn(1:n, 2)) - 1
  folds <- sample(1:n_fold, n_pairs, replace=TRUE)
  
  theta_init <- cmdscale(dist_mat, max(K_cand))
  cv_loss <- matrix(0, n_K)
  for (index_K in 1:n_K) {
    K <- K_cand[index_K]
    for (index_fold in 1:n_fold) {
      index_valid <- which(folds==index_fold)
      index_train <- which(folds!=index_fold)
      valid_set <- all_pairs[index_valid,]
      train_set <- all_pairs[index_train,]
      
      theta <- theta_init[,1:K]
      mds_res <- MDS_subset(dist_mat, theta, max_epoch, step_size, tot, train_set, valid_set)
      cv_loss[index_K] <- cv_loss[index_K] + mds_res$valid_loss
    }
  }
  
  if (return_dist) res <- list(K=K_cand[which.min(cv_loss)], K_cand=K_cand, cv_loss=cv_loss, dist_mat=dist_mat)
  else res <- list(K=K_cand[which.min(cv_loss)], K_cand=K_cand, cv_loss=cv_loss)
  
  res
}

#' Feature Extraction by MDS for Large Dataset
#' 
#' \code{seq2feature_mds_large} extracts MDS features from a large number of 
#' response processes. The algorithm proposed in Paradis (2018) is implemented with minor 
#' variations to perform MDS. The algorithm first selects a relatively small subset of 
#' response processes to perform the classical MDS. Then the coordinate of each of the 
#' other response processes are obtained by minimizing the loss function related to the target
#' response processes and the those in the subset through BFGS.
#' 
#' @family feature extraction methods
#' @inheritParams seq2feature_mds
#' @param seqs an object of class \code{"\link{proc}"}
#' @param subset_size the size of the subset on which classical MDS is performed.
#' @param subset_method a character string specifying the method for choosing the subset.
#'   It must be one of \code{"random"}, \code{"sample_avgmax"},
#'   \code{"sample_minmax"}, \code{"full_avgmax"}, and \code{"full_minmax"}.
#' @param n_cand The size of the candidate set when selecting the subset. It is only used when 
#' \code{subset_method} is \code{"sample_avgmax"} or \code{"sample_minmax"}.
#' @return \code{seq2feature_mds_large} returns an \eqn{n \times K} matrix of extracted 
#' features.
#' @references Paradis, E. (2018). Multidimensional Scaling with Very Large Datasets. 
#'   \emph{Journal of Computational and Graphical Statistics}, 27, 935--939.
#' @export
seq2feature_mds_large <- function(seqs, K, dist_type = "oss_action", subset_size, 
                                  subset_method = "random", n_cand = 10, 
                                  pca = TRUE, seed = 12345, L_set = 1:3) {
  n <- length(seqs$action_seqs)
  theta <- matrix(0, n, K)
  
  if (subset_method == "random") {
    init_obj <- sample(1:n, subset_size)
    remn_obj <- setdiff(1:n, init_obj)
  } else if (subset_method == "sample_avgmax") {
    init_obj <- sample(1:n, 1)
    remn_obj <- setdiff(1:n, init_obj)
    
    while (length(init_obj) < subset_size) {
      n_remn <- length(remn_obj)
      n_init <- length(init_obj)
      if (n_remn < n_cand) n_cand <- n_remn
      candidate_obj <- sample(remn_obj, n_cand)
      d_mat <- matrix(0, n_init, n_cand)
      for (i in 1:n_cand) {
        for (j in 1:n_init) {
          if (dist_type == "oss_action") d_mat[j, i] <- calculate_dissimilarity_cpp(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]])
          if (dist_type == "oss_both") d_mat[j, i] <- calculate_tdissimilarity_cpp(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]], seqs$time_seqs[[candidate_obj[i]]], seqs$time_seqs[[init_obj[j]]])
          if (dist_type == "ngram") d_mat[j, i] <- calculate_ngram_dissimilarity(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]], L_set)
        }
      }
      d_avg <- colSums(d_mat) / n_init
      o <- order(d_avg, decreasing = TRUE)
      current_obj <- candidate_obj[o[1]]
      init_obj <- c(init_obj, current_obj)
      remn_obj <- setdiff(remn_obj, current_obj)
    }
    
  } else if (subset_method == "sample_minmax") {
    init_obj <- sample(1:n, 1)
    remn_obj <- setdiff(1:n, init_obj)
    
    while (length(init_obj) < subset_size) {
      n_remn <- length(remn_obj)
      n_init <- length(init_obj)
      if (n_remn < n_cand) n_cand <- n_remn
      candidate_obj <- sample(remn_obj, n_cand)
      d_mat <- matrix(0, n_init, n_cand)
      for (i in 1:n_cand) {
        for (j in 1:n_init) {
          if (dist_type == "oss_action") d_mat[j, i] <- calculate_dissimilarity_cpp(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]])
          if (dist_type == "oss_both") d_mat[j, i] <- calculate_tdissimilarity_cpp(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]], seqs$time_seqs[[candidate_obj[i]]], seqs$time_seqs[[init_obj[j]]])
          if (dist_type == "ngram") d_mat[j, i] <- calculate_ngram_dissimilarity(seqs$action_seqs[[candidate_obj[i]]], seqs$action_seqs[[init_obj[j]]], L_set)
        }
      }
      d_min <- apply(d_mat, 2, min)
      o <- order(d_min, decreasing = TRUE)
      current_obj <- candidate_obj[o[1]]
      init_obj <- c(init_obj, current_obj)
      remn_obj <- setdiff(remn_obj, current_obj)
    }
  } else if (subset_method == "full_avgmax") {
    init_obj <- sample(1:n, 1)
    remn_obj <- setdiff(1:n, init_obj)
    current_obj <- init_obj
    d_mat <- numeric(0)

    while (length(init_obj) < subset_size) {
      n_remn <- length(remn_obj)
      n_init <- length(init_obj)
      d_new <- rep(0, n_remn)
      for (i in 1:n_remn) {
        if (dist_type == "oss_action") d_new[i] <- calculate_dissimilarity_cpp(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]])
        if (dist_type == "oss_both") d_new[i] <- calculate_tdissimilarity_cpp(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]], seqs$time_seqs[[remn_obj[i]]], seqs$time_seqs[[current_obj]])
        if (dist_type == "ngram") d_new[i] <- calculate_ngram_dissimilarity(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]], L_set)
      }
      d_mat <- rbind(d_mat, d_new)
      
      d_avg <- colSums(d_mat) / n_init
      o <- order(d_avg, decreasing = TRUE)
      current_obj <- remn_obj[o[1]]
      init_obj <- c(init_obj, current_obj)
      remn_obj <- setdiff(remn_obj, current_obj)
      d_mat <- d_mat[,-o[1]]
    }
  } else if (subset_method == "full_minmax") {
    init_obj <- sample(1:n, 1)
    remn_obj <- setdiff(1:n, init_obj)
    current_obj <- init_obj
    d_mat <- numeric(0)
    
    while (length(init_obj) < subset_size) {
      n_remn <- length(remn_obj)
      n_init <- length(init_obj)
      d_new <- rep(0, n_remn)
      for (i in 1:n_remn) {
        if (dist_type == "oss_action") d_new[i] <- calculate_dissimilarity_cpp(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]])
        if (dist_type == "oss_both") d_new[i] <- calculate_tdissimilarity_cpp(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]], seqs$time_seqs[[remn_obj[i]]], seqs$time_seqs[[current_obj]])
        if (dist_type == "ngram") d_new[i] <- calculate_ngram_dissimilarity(seqs$action_seqs[[remn_obj[i]]], seqs$action_seqs[[current_obj]], L_set)
      }
      d_mat <- rbind(d_mat, d_new)
      
      d_min <- apply(d_mat, 2, min) 
      o <- order(d_min, decreasing = TRUE)
      current_obj <- remn_obj[o[1]]
      init_obj <- c(init_obj, current_obj)
      remn_obj <- setdiff(remn_obj, current_obj)
      d_mat <- d_mat[,-o[1]]
    }
  } else {
    stop("Undefined subset method!\n")
  }
  
  D <- matrix(0, subset_size, subset_size)
  if (dist_type == "oss_action") D <- calculate_dist_cpp(seqs$action_seqs[init_obj])
  else if (dist_type == "oss_both") D <- calculate_tdist_cpp(seqs$action_seqs[init_obj], seqs$time_seqs[init_obj])
  else if (dist_type == "ngram") D <- calculate_ngram_dist_cpp(seqs$action_seqs[init_obj], L_set)
  
  theta[init_obj, ] <- cmdscale(D, k = K)
  
  obj_fun <- function(theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- cbind(rep(1, subset_size)) %*% t(theta_j) 
    theta_diff <- theta_m_mat - theta_j_mat
    res <- sum((d_vec - sqrt(rowSums((theta_diff)^2)))^2)
    res
  }
  grad_fun <- function(theta_j, theta_m_mat, d_vec) {
    theta_j_mat <- cbind(rep(1, subset_size)) %*% t(theta_j)
    theta_diff <- theta_m_mat - theta_j_mat
    dhat_vec <- sqrt(rowSums((theta_diff)^2))
    res <- 2 * colSums((d_vec / dhat_vec - 1) * theta_diff)
    res
  }
  theta_m_mat <- theta[init_obj, ]
  for (i in remn_obj) {
    d_vec <- rep(0, subset_size)
    for (j in 1:subset_size) {
      if (dist_type == "oss_action") d_vec[j] <- calculate_dissimilarity_cpp(seqs$action_seqs[[i]], seqs$action_seqs[[init_obj[j]]])
      else if (dist_type == "oss_both") d_vec[j] <- calculate_tdissimilarity_cpp(seqs$action_seqs[[i]], seqs$action_seqs[[init_obj[j]]], seqs$time_seqs[[i]], seqs$time_seqs[[init_obj[j]]])
      else if (dist_type == "ngram") d_vec[j] <- calculate_ngram_dissimilarity(seqs$action_seqs[[i]], seqs$action_seqs[[init_obj[j]]], L_set)
    }
    opt_res <- optim(rnorm(K), fn = obj_fun, gr = grad_fun, method = "BFGS", theta_m_mat = theta_m_mat, d_vec = d_vec)
    theta[i, ] <- opt_res$par
  }
  
  if (pca) theta <- prcomp(theta)$x
  
  theta
}

#' Feature extraction by stochastic mds
#' @param seqs a \code{"\link{proc}"} object or a square matrix. If a squared matrix is 
#'   provided, it is treated as the dissimilary matrix of a group of response processes.
#' @param K the number of features to be extracted.
#' @param dist_type a character string specifies the dissimilarity measure for two
#'   response processes. See 'Details'.
#' @param max_epoch the maximum number of epochs for stochastic gradient
#'   descent.
#' @param step_size the step size of stochastic gradient descent.
#' @param pca a logical scalar. If \code{TRUE}, the principal components of the
#'   extracted features are returned.
#' @param tot the accuracy tolerance for determining convergence.
#' @param return_dist logical. If \code{TRUE}, the dissimilarity matrix will be
#'   returned. Default is \code{FALSE}.
#' @param seed random seed.
#' @param L_set length of ngrams considered.
#' @return \code{seq2feature_mds_stochastic} returns a list containing 
#'   \item{theta}{a numeric matrix giving the \code{K} extracted features or principal
#'   features. Each column is a feature.} 
#'   \item{loss}{the value of the multidimensional scaling objective function.}
#'   \item{dist_mat}{the dissimilary matrix. This element exists only if \code{return_dist=TRUE}.}
#' @export
seq2feature_mds_stochastic <- function(seqs = NULL, K = 2, dist_type = "oss_action", 
                                       max_epoch=100, step_size=0.01, pca=TRUE, 
                                       tot=1e-6, return_dist=FALSE, seed=12345, L_set=1:3) {
  set.seed(seed)
  if (is.null(seqs)) 
    stop("Either response processes or their dissimilarity matrix should be provided!\n")
  if (is.matrix(seqs)) {
    if (nrow(seqs) != ncol(seqs)) stop("Provided matrix is not square!\n")
    dist_mat <- seqs
    n <- nrow(dist_mat)
  } else if (class(seqs) == "proc") {
    n <- length(seqs$action_seqs)
    dist_mat <- matrix(0, n, n)
    if (is.null(seqs$time_seqs) & dist_type == "oss_both") {
      warning("Timestamp sequences are not available. 
              change method from 'oss_both' to 'oss_action'!\n")
      dist_type <- "oss_action"
    }
    if (dist_type == "oss_action") {
      dist_mat <- calculate_dist_cpp(seqs$action_seqs)
    } else if (dist_type == "oss_both") {
      dist_mat <- calculate_tdist_cpp(seqs$action_seqs, seqs$time_seqs)
    } else if (dist_type == "ngram") {
      dist_mat <- calculate_ngram_dist_cpp(seqs$action_seqs, L_set)
    } else stop("Invalid dissimilarity method!\n")
  } else {
    stop("seqs should be a 'proc' object or a square matrix\n!")
  }
  
  # initialize
  theta <- cmdscale(dist_mat, K)
  
  # mds
  mds_res <- MDS(dist_mat, theta, max_epoch, step_size, tot, seed)
  if (!mds_res$convergence) warning("MDS does not converge!")
  if (pca) theta <- prcomp(theta, center=TRUE, scale=FALSE)$x
  
  if (return_dist) res <- list(theta=theta, loss=mds_res$loss, dist_mat=dist_mat)
  else res <- list(theta=theta, loss=mds_res$loss)
  
  res
}



#' Feature Extraction by autoencoder
#'
#' \code{seq2feature_seq2seq} extract features from response processes by autoencoder.
#'
#' This function wraps \code{\link{aseq2feature_seq2seq}}, 
#' \code{\link{tseq2feature_seq2seq}}, and \code{\link{atseq2feature_seq2seq}}.  
#' 
#' @family feature extraction methods
#' @param seqs an object of class \code{"\link{proc}"}.
#' @param ae_type a string specifies the type of autoencoder. The autoencoder can be an
#' action sequence autoencoder ("action"), a time sequence autoencoder ("time"), or an 
#' action-time sequence autoencoder ("both").
#' @param cumulative logical. If TRUE, the sequence of cumulative time up to each event is  
#'  used as input to the neural network. If FALSE, the sequence of inter-arrival time (gap 
#'  time between an event and the previous event) will be used as input to the neural network.
#'  Default is FALSE.
#' @param log logical. If TRUE, for the timestamp sequences, input of the neural net is
#'  the base-10 log of the original sequence of times plus 1 (i.e., log10(t+1)). If FALSE,
#'  the original sequence of times is used.
#' @param weights a vector of 2 elements for the weight of the loss of action sequences
#'  (categorical_crossentropy) and time sequences (mean squared error), respectively. 
#'  The total loss is calculated as the weighted sum of the two losses.
#' @param K the number of features to be extracted.
#' @param rnn_type the type of recurrent unit to be used for modeling
#'   response processes. \code{"lstm"} for the long-short term memory unit. 
#'   \code{"gru"} for the gated recurrent unit.
#' @param n_epoch the number of training epochs for the autoencoder.
#' @param method the method for computing features from the output of an
#'   recurrent neural network in the encoder. Available options are 
#'   \code{"last"} and \code{"avg"}.
#' @param step_size the learning rate of optimizer.
#' @param optimizer_name a character string specifying the optimizer to be used
#'   for training. Availabel options are \code{"sgd"}, \code{"rmsprop"}, 
#'   \code{"adadelta"}, and \code{"adam"}.
#' @param samples_train,samples_valid,samples_test vectors of indices specifying the
#'   training, validation and test sets for training autoencoder.
#' @param pca logical. If TRUE, the principal components of features are
#'   returned. Default is TRUE.
#' @param gpu logical. If TRUE, use gpu for training when available. Default is FALSE.
#' @param parallel logical. If TRUE, allow cpu parallel computing. Default is FALSE.
#' @param seed random seed.
#' @param verbose logical. If TRUE, training progress is printed.
#' @param return_theta logical. If TRUE, extracted features are returned.
#' @return \code{seq2feature_seq2seq} returns a list containing
#'   \item{theta}{a matrix containing \code{K} features or principal features. Each column is a feature.}
#'   \item{train_loss}{a vector of length \code{n_epoch} recording the trace of training losses.}
#'   \item{valid_loss}{a vector of length \code{n_epoch} recording the trace of validation losses.}
#'   \item{test_loss}{a vector of length \code{n_epoch} recording the trace of test losses. Exists only if \code{samples_test} is not \code{NULL}.}
#' @seealso \code{\link{chooseK_seq2seq}} for choosing \code{K} through cross-validation.
#' @examples
#' \donttest{ 
#' n <- 50
#' data(cc_data)
#' samples <- sample(1:length(cc_data$seqs$time_seqs), n)
#' seqs <- sub_seqs(cc_data$seqs, samples)
#' 
#' # action sequence autoencoder
#' K_res <- chooseK_seq2seq(seqs=seqs, ae_type="action", K_cand=c(5, 10), 
#'                          n_epoch=5, n_fold=2, valid_prop=0.2)
#' seq2seq_res <- seq2feature_seq2seq(seqs=seqs, ae_type="action", K=K_res$K, 
#'                        n_epoch=5, samples_train=1:40, samples_valid=41:50)
#' theta <- seq2seq_res$theta
#' 
#' # time sequence autoencoder
#' K_res <- chooseK_seq2seq(seqs=seqs, ae_type="time", K_cand=c(5, 10), 
#'                          n_epoch=5, n_fold=2, valid_prop=0.2)
#' seq2seq_res <- seq2feature_seq2seq(seqs=seqs, ae_type="time", K=K_res$K, 
#'                        n_epoch=5, samples_train=1:40, samples_valid=41:50)
#' theta <- seq2seq_res$theta
#' 
#' # action and time sequence autoencoder
#' K_res <- chooseK_seq2seq(seqs=seqs, ae_type="both", K_cand=c(5, 10), 
#'                          n_epoch=5, n_fold=2, valid_prop=0.2)
#' seq2seq_res <- seq2feature_seq2seq(seqs=seqs, ae_type="both", K=K_res$K, 
#'                        n_epoch=5, samples_train=1:40, samples_valid=41:50)
#' theta <- seq2seq_res$theta
#' plot(seq2seq_res$train_loss, col="blue", type="l")
#' lines(seq2seq_res$valid_loss, col="red")
#' }
#' @export
seq2feature_seq2seq <- function(seqs, ae_type="action", K, rnn_type="lstm", n_epoch=50, 
                                method="last", step_size=0.0001, optimizer_name="adam", 
                                cumulative=FALSE, log=TRUE, weights=c(1.0, 0.5), 
                                samples_train, samples_valid, samples_test=NULL, pca=TRUE, 
                                gpu=FALSE, parallel=FALSE, seed=12345, verbose=TRUE, 
                                return_theta=TRUE) {
  if (ae_type=="action")
    res <- aseq2feature_seq2seq(aseqs=seqs$action_seqs, 
                                K=K,
                                rnn_type=rnn_type,
                                n_epoch=n_epoch,
                                method=method,
                                step_size = step_size,
                                optimizer_name = optimizer_name,
                                samples_train = samples_train,
                                samples_valid = samples_valid,
                                samples_test = samples_test,
                                pca = pca,
                                gpu = gpu,
                                parallel = parallel,
                                seed = seed,
                                verbose = verbose,
                                return_theta = TRUE)
  else if (ae_type=="time")
    res <- tseq2feature_seq2seq(tseqs=seqs$time_seqs, 
                                K=K,
                                cumulative = cumulative,
                                log = log,
                                rnn_type=rnn_type,
                                n_epoch=n_epoch,
                                method=method,
                                step_size = step_size,
                                optimizer_name = optimizer_name,
                                samples_train = samples_train,
                                samples_valid = samples_valid,
                                samples_test = samples_test,
                                pca = pca,
                                gpu = gpu,
                                parallel = parallel,
                                seed = seed,
                                verbose = verbose,
                                return_theta = TRUE)
  else if (ae_type=="both")
    res <- atseq2feature_seq2seq(atseqs=seqs,
                                 K=K,
                                 weights = weights,
                                 cumulative = cumulative,
                                 log = log,
                                 rnn_type=rnn_type,
                                 n_epoch=n_epoch,
                                 method=method,
                                 step_size = step_size,
                                 optimizer_name = optimizer_name,
                                 samples_train = samples_train,
                                 samples_valid = samples_valid,
                                 samples_test = samples_test,
                                 pca = pca,
                                 gpu = gpu,
                                 parallel = parallel,
                                 seed = seed,
                                 verbose = verbose,
                                 return_theta = TRUE)
  else stop("ae_type has to be 'action', time', or 'both'!\n")
  
  res
}

#' Choose the number of autoencoder features
#'
#' \code{chooseK_seq2seq} chooses the number of features to be extracted
#'  by cross-validation.
#' 
#' @inheritParams seq2feature_seq2seq
#' @param K_cand the candidates of the number of features.
#' @param n_fold the number of folds for cross-validation.
#' @param valid_prop the proportion of validation samples in each fold.
#' @return \code{chooseK_seq2seq} returns a list containing
#'   \item{K}{the candidate in \code{K_cand} producing the smallest cross-validation loss.}
#'   \item{K_cand}{the candidates of number of features.}
#'   \item{cv_loss}{the cross-validation loss for each candidate in \code{K_cand}.}
#' @seealso \code{\link{seq2feature_seq2seq}} for feature extraction given the number of features.
#' @export
chooseK_seq2seq <- function(seqs, ae_type, K_cand, rnn_type="lstm", n_epoch=50, method="last", 
                            step_size=0.0001, optimizer_name="adam", n_fold=5, 
                            cumulative = FALSE, log = TRUE, weights = c(1., .5), 
                            valid_prop=0.1, gpu = FALSE, parallel=FALSE, seed=12345, 
                            verbose=TRUE) {
  set.seed(seed)
  n_K <- length(K_cand)
  n_seq <- length(seqs$action_seqs)
  folds <- sample(1:n_fold, n_seq, replace=TRUE)
  
  cv_loss <- matrix(0, n_K)
  for (index_K in 1:n_K) {
    K <- K_cand[index_K]
    if (verbose) cat("Candidate K:", K, "\n")
    for (index_fold in 1:n_fold) {
      index_test <- which(folds==index_fold)
      index_train_valid <- which(folds!=index_fold)
      index_valid <- sample(index_train_valid, round(length(index_train_valid)*valid_prop))
      index_train <- setdiff(index_train_valid, index_valid)
      
      seq2seq_res <- seq2feature_seq2seq(seqs = seqs, 
                                         ae_type = ae_type, 
                                         K = K, 
                                         rnn_type = rnn_type, 
                                         n_epoch = n_epoch, 
                                         method = method, 
                                         step_size = step_size, 
                                         optimizer_name = optimizer_name, 
                                         cumulative = cumulative,
                                         log = log,
                                         weights = weights,
                                         samples_train = index_train, 
                                         samples_valid = index_valid, 
                                         samples_test = index_test, 
                                         pca = FALSE, 
                                         gpu = gpu, 
                                         parallel = parallel, 
                                         seed = seed,
                                         verbose = verbose, 
                                         return_theta = FALSE)
      
      cv_loss[index_K] <- cv_loss[index_K] + seq2seq_res$test_loss[which.min(seq2seq_res$valid_loss)]
    }
  }
  
  res <- list(K=K_cand[which.min(cv_loss)], K_cand=K_cand, cv_loss=cv_loss)
}


 

