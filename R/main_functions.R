
#' Feature extraction via multidimensional scaling
#'
#' \code{seq2feature_mds} extracts \code{K} features from action sequences by
#' multidimensional scaling.
#'
#' This function minimizes the objective function by stochastic gradient
#' descent. The coordinates of the objects are extracted features. The number of
#' features to be extracted \code{K} can be selected by cross-validation using
#' \code{\link{chooseK_mds}}.
#'
#' @family feature extraction methods
#' @param seqs a list or a square matrix. If a list is provided, each of its
#'   elements should be an action sequence in the form of a vector of actions.
#'   If a squared matrix is provided, it is treated as the dissimilary matrix of
#'   a group of action sequences.
#' @param K the number of features to be extracted.
#' @param method a character string specifies the dissimilarity measure for two
#'   action sequences. Order-based sequence similarity in Gomez-Alonso and Valls
#'   (2008) (\code{"oss"}) is the only available method for now.
#' @param max_epoch the maximum number of epochs for stochastic gradient
#'   descent.
#' @param step_size the step size of stochastic gradient descent.
#' @param pca a logical scalar. If \code{TRUE}, the principal components of the
#'   extracted features are returned.
#' @param tot the accuracy tolerance for determining convergence.
#' @param return_dist logical. If \code{TRUE}, the dissimilarity matrix will be
#'   returned. Default is \code{FALSE}.
#' @return \code{seq2feature_mds} returns a list containing 
#'   \item{theta}{a numeric matrix giving the \code{K} extracted features or principal
#'   features. Each column is a feature.} 
#'   \item{loss}{the value of the multidimensional scaling objective function.}
#'   \item{dist_mat}{the dissimilary matrix. This element exists only if \code{return_dist=TRUE}.}
#' @seealso \code{\link{chooseK_mds}} for choosing \code{K}.
#' @references Gomez-Alonso, C. and Valls, A. (2008). A similarity measure for sequences of
#'   categorical data based on the ordering of common elements. In V. Torra & Y. Narukawa (Eds.) 
#'   \emph{Modeling Decisions for Artificial Intelligence}, (pp. 134-145). Springer Berlin Heidelberg.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' theta <- seq2feature_mds(seqs, 5)$theta
#' @export
seq2feature_mds <- function(seqs=NULL, K=2, method="oss", max_epoch=100, step_size=0.01, pca=TRUE, tot=1e-6, return_dist=FALSE) {
  if (is.null(seqs)) stop("Either a list of sequences or their dissimilarity matrix should be provided!\n")
  if (is.matrix(seqs)) {
    if (nrow(seqs) != ncol(seqs)) stop("Provided matrix is not square!\n")
    dist_mat <- seqs
    n <- nrow(dist_mat)
  } else if (is.list(seqs)) {
    n <- length(seqs)
    dist_mat <- matrix(0, n, n)
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        seq1 <- seqs[[i]]
        seq2 <- seqs[[j]]	
        dist_mat[i,j] <- calculate_dissimilarity(seq1, seq2, method=method) 
        dist_mat[j,i] <- dist_mat[i,j]
      }
    }
  } else {
    stop("seqs should be a list or a square matrix\n!")
  }
  
  # initialize
  theta <- cmdscale(dist_mat, K)
  
  # mds
  mds_res <- MDS(dist_mat, theta, max_epoch, step_size, tot)
  if (!mds_res$convergence) warning("MDS does not converge!")
  if (pca) theta <- prcomp(theta, center=TRUE, scale=FALSE)$x
  
  if (return_dist) res <- list(theta=theta, loss=mds_res$loss, dist_mat=dist_mat)
  else res <- list(theta=theta, loss=mds_res$loss)
  
  res
}

#' Choose the number of multidimensional scaling features
#' 
#' \code{chooseK_mds} choose the number of multidimensional scaling features
#'   to be extracted by cross-validation.
#'
#' @param K_cand the candidates of the number of features.
#' @param n_fold the number of folds for cross-validation.
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
chooseK_mds <- function(seqs=NULL, K_cand, method="oss", n_fold=5, max_epoch=100, step_size=0.01, tot=1e-6, return_dist=FALSE) {
  if (is.null(seqs)) stop("Either a list of sequences or their dissimilarity matrix should be provided!\n")
  if (is.matrix(seqs)) {
    if (nrow(seqs) != ncol(seqs)) stop("Provided matrix is not square!\n")
    dist_mat <- seqs
    n <- nrow(dist_mat)
  } else if (is.list(seqs)) {
    n <- length(seqs)
    dist_mat <- matrix(0, n, n)
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        seq1 <- seqs[[i]]
        seq2 <- seqs[[j]]	
        dist_mat[i,j] <- calculate_dissimilarity(seq1, seq2, method=method) 
        dist_mat[j,i] <- dist_mat[i,j]
      }
    }
  } else {
    stop("seqs should be a list or a square matrix!\n")
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

#' Feature Extraction by action sequence autoencoder
#'
#' \code{seq2feature_seq2seq} extract features from action sequences by action
#' sequence autoencoder.
#'
#' This function trains a sequence-to-sequence autoencoder using keras. The encoder
#' of the autoencoder consists of an embedding layer and a recurrent neural network.
#' The decoder consists of another recurrent neural network and a fully connect layer
#' with softmax activation. The outputs of the encoder are the extracted features.
#' 
#' The output of the encoder is a function of the encoder recurrent neural network.
#' It is the last output of the encoder recurrent neural network if \code{method="last"}
#' and the average of the encoder recurrent nenural network if \code{method="avg"}.
#' 
#' 
#' @family feature extraction methods
#' @param seqs a list of \code{n} action sequences. Each element is an action
#'   sequence in the form of a vector of actions.
#' @param K the number of features to be extracted.
#' @param rnn_type the type of recurrent unit to be used for modeling
#'   action sequences. \code{"lstm"} for the long-short term memory unit. 
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
#'   returned.
#' @param gpu logical. If TRUE, use gpu for training if available.
#' @param verbose logical. If TRUE, training progress is printed.
#' @param return_theta logical. If TRUE, extracted features are returned.
#' @return \code{seq2feature_seq2seq} returns a list containing
#'   \item{theta}{a matrix containing \code{K} features or principal features. Each column is a feature.}
#'   \item{train_loss}{a vector of length \code{n_epoch} recording the trace of training losses.}
#'   \item{valid_loss}{a vector of length \code{n_epoch} recording the trace of validation losses.}
#'   \item{test_loss}{a vector of length \code{n_epoch} recording the trace of test losses. Exists only if \code{samples_test} is not \code{NULL}.}
#' @seealso \code{\link{chooseK_seq2seq}} for choosing \code{K} through cross-validation.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' seq2seq_res <- seq2feature_seq2seq(seqs, 5, rnn_type="lstm", n_epoch=5, samples_train=1:40, samples_valid=41:50)
#' features <- seq2seq_res$theta
#' plot(seq2seq_res$train_loss, col="blue", type="l")
#' lines(seq2seq_res$valid_loss, col="red")
#' @export
seq2feature_seq2seq <- function(seqs, K, rnn_type="lstm", n_epoch=50, method="last", step_size=0.0001, optimizer_name="adam", samples_train, samples_valid, samples_test=NULL, pca=TRUE, gpu=FALSE, verbose=TRUE, return_theta=TRUE) {
  
  if (!(rnn_type %in% c("lstm", "gru"))) 
    stop("Invalid rnn_type! Available options: lstm, gru.\n")
  if (!(method %in% c("last", "avg"))) 
    stop("Invalid method! Available options: last, avg.\n")
  if (!(optimizer_name %in% c("sgd", "adam", "rmsprop", "adadelta"))) 
    stop("Invalid optimizer! Available options: sgd, adam, rmsprop, adadelta.\n")
  n_seq <- length(seqs)
  
  events <- unique(unlist(seqs))
  n_event <- length(events)

  # convert action sequence to one-hot vectors
  int_seqs <- list()
  target_seqs <- list()
  
  for (index_seq in 1:n_seq) {
    my_seq <- seqs[[index_seq]]
    n_l <- length(my_seq)
    onehot_mat <- matrix(0, n_l, n_event)
    tmp <- match(my_seq, events)
    int_seqs[[index_seq]] <- matrix(tmp - 1, 1, n_l)
    onehot_mat[cbind(1:n_l, tmp)] <- 1
    target_seqs[[index_seq]] <- array(onehot_mat, dim=c(1,n_l, n_event))
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # Define keras model
  # Define an input sequence and process it.
  encoder_inputs <- layer_input(shape=list(NULL))
  if (rnn_type == "lstm")
  {
    if (method=="last") {
      encoder_outputs_long <- encoder_inputs %>%
        layer_embedding(n_event, K) %>%
        layer_lstm(units=K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[3]]
      
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[3]], axis=2)
        
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_outputs_long <- encoder_inputs %>%
        layer_embedding(n_event, K) %>%
        layer_lstm(units=K, return_sequences=TRUE)
      
      fn_average <- function(x) k_mean(x, axis=2)
      encoder_outputs <- encoder_outputs_long %>% layer_lambda(fn_average)
      # Repeat 2D tensor to form a 3D tensor to be used as inputs of LSTM layer
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(k_expand_dims(x[,,1], axis = -1)) # matrix with ones, shaped as (batch, steps, 1)
        latentMatrix <- k_mean(x, axis = 2, keepdims = TRUE) # latent vars, shaped as (batch, 1, latent_dim)
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    }
    
    decoder_inputs <-  encoder_outputs_long %>% layer_lambda(fn_rp)
    decoder_outputs <- decoder_inputs %>% 
      layer_lstm(units=K, return_sequences=TRUE) %>% 
      layer_dense(n_event, activation='softmax')
  } else if (rnn_type == "gru") {
    if (method=="last") {
      encoder_outputs_long <- encoder_inputs %>%
        layer_embedding(n_event, K) %>%
        layer_gru(units=K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[2]]
      
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[2]], axis=2)
        
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_outputs_long <- encoder_inputs %>%
        layer_embedding(n_event, K) %>%
        layer_gru(units=K, return_sequences=TRUE)
      
      fn_average <- function(x) k_mean(x, axis=2)
      encoder_outputs <- encoder_outputs_long %>% layer_lambda(fn_average)
      # Repeat 2D tensor to form a 3D tensor to be used as inputs of LSTM layer
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(k_expand_dims(x[,,1], axis = -1)) # matrix with ones, shaped as (batch, steps, 1)
        latentMatrix <- k_mean(x, axis = 2, keepdims = TRUE) # latent vars, shaped as (batch, 1, latent_dim)
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    }
    
    decoder_inputs <-  encoder_outputs_long %>% layer_lambda(fn_rp)
    decoder_outputs <- decoder_inputs %>% 
      layer_gru(units=K, return_sequences=TRUE) %>% 
      layer_dense(n_event, activation='softmax')
  }
  
  autoencoder_model <- keras_model(inputs = encoder_inputs, outputs = decoder_outputs)
  encoder_model <- keras_model(inputs = encoder_inputs, outputs = encoder_outputs)
  
  # Run training
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  autoencoder_model %>% compile(optimizer=optimizer, loss='categorical_crossentropy')
  
  best_valid_loss <- Inf
  valid_loss <- rep(0, n_epoch)
  train_loss <- rep(0, n_epoch)
  if (!is.null(samples_test)) test_loss <- rep(0, n_epoch)
  
  theta <- matrix(0, n_seq, K)
  
  for (index_epoch in 1:n_epoch) {
    if (verbose) cat("Epoch", index_epoch, "\n")
    samples_train <- sample(samples_train)
    
    for (index_seq in samples_train) {
      autoencoder_model %>% train_on_batch(int_seqs[[index_seq]], target_seqs[[index_seq]])
    }
    for (index_seq in samples_train) {
      train_loss[index_epoch] <- train_loss[index_epoch] + evaluate(autoencoder_model, int_seqs[[index_seq]], target_seqs[[index_seq]], verbose=0)
    }
    for (index_seq in samples_valid) {
      valid_loss[index_epoch] <- valid_loss[index_epoch] + evaluate(autoencoder_model, int_seqs[[index_seq]], target_seqs[[index_seq]], verbose=0)
    }
    for (index_seq in samples_test) {
      test_loss[index_epoch] <- test_loss[index_epoch] + evaluate(autoencoder_model, int_seqs[[index_seq]], target_seqs[[index_seq]], verbose=0)
    }
    if (valid_loss[index_epoch] < best_valid_loss) {
      best_valid_loss <- valid_loss[index_epoch]
      if (return_theta) {
        for(index_seq in 1:n_seq) theta[index_seq,] <- predict(encoder_model, int_seqs[[index_seq]])
      }
    }
  }
  
  k_clear_session()
  res <- list(train_loss=train_loss, valid_loss=valid_loss)
  if (!is.null(samples_test)) res$test_loss <- test_loss
  if (return_theta) {
    if (pca) theta <- prcomp(theta, center=FALSE, scale=FALSE)$x
    res$theta <- theta
  }
  
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
#' @examples 
#' n <- 50
#' seqs <- seq_gen(n)
#' K_res <- chooseK_seq2seq(seqs, K_cand=c(5, 10), rnn_type="lstm", n_epoch=5, n_fold=2, valid_prop=0.2)
#' seq2seq_res <- seq2feature_seq2seq(seqs, K_res$K, rnn_type="lstm", n_epoch=10, samples_train=1:40, samples_valid=41:50)
#' theta <- seq2seq_res$theta
#' @export
chooseK_seq2seq <- function(seqs, rnn_type="lstm", K_cand, n_epoch=50, method="last", step_size=0.0001, optimizer_name="adam", n_fold=5, valid_prop=0.1, gpu = FALSE, verbose=TRUE) {
  n_K <- length(K_cand)
  n_seq <- length(seqs)
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
      
      seq2seq_res <- seq2feature_seq2seq(seqs = seqs, K = K, rnn_type = rnn_type, n_epoch = n_epoch, method = method, step_size = step_size, optimizer_name = optimizer_name, samples_train = index_train, samples_valid = index_valid, samples_test = index_test, pca = FALSE, gpu = gpu, verbose = verbose, return_theta = FALSE)
      
      cv_loss[index_K] <- cv_loss[index_K] + seq2seq_res$test_loss[which.min(seq2seq_res$valid_loss)]
    }
  }
  
  res <- list(K=K_cand[which.min(cv_loss)], K_cand=K_cand, cv_loss=cv_loss)
}


#' Fitting sequence models
#'
#' \code{seqm} is used to fit a neural network model relating action sequences
#' with a response variable.
#'
#' The model consists of an embedding layer, a recurrent layer and one or more
#' fully connected layers. The embedding layer takes a action sequence and
#' output a sequences of \code{K} dimensional numeric vectors to the recurrent
#' layer. The last output of the recurrent layer is used as the input of the
#' subsequent fully connected layers. If \code{response_type="binary"}, the last
#' layer uses the sigmoid activation to produce the probability of the response
#' being positive. If \code{response_type="scale"}, the last layer uses the linear
#' activation. The dimension of the output of other fully connected layers
#' (if any) is specified by \code{K_hidden}.
#'
#' The action sequences are re-coded into integer sequences and are padded with
#' zeros to length \code{max_len} before feeding into the model. If the provided
#' \code{max_len} is smaller than the length of the longest sequence in
#' \code{seqs}, it will be overridden.
#'
#' @inheritParams seq2feature_seq2seq
#' @param formula an object of class \code{"\link{formula}"} (or one that can be coerced
#'   to that class): a symbolic description of the model to be fitted.
#' @param response_type "binary" or "scale".
#' @param actions a character vector gives all possible actions. It is will be
#'   expanded to include all actions appear in \code{seqs} if necessary.
#' @param data a dataframe containing the variables in the model.
#' @param n_hidden the number of hidden fully-connected layers.
#' @param K_emb the latent dimension of the embedding layer.
#' @param K_rnn the latent dimension of the recurrent neural network.
#' @param K_hidden a vector of length \code{n_hidden} specifying the number of
#'   nodes in each hidden layer.
#' @param n_epoch the number of training epochs.
#' @param batch_size the batch size used in training.
#' @param valid_split proportion of sequences used as the validation set. 
#' @param index_valid a vector of indices specifying the validation set.
#' @param max_len the maximum length of input sequences.
#' @return \code{seqm} returns an object of class \code{"seqm"}, which is a list containing
#'   \item{formula}{the model formula.}
#'   \item{model_fit}{a vector of class \code{"raw"}. It is the serialized version of 
#'     the trained keras model.} 
#'   \item{actions}{all possible actions.}
#'   \item{max_len}{the maximum length of action sequences.}
#'   \item{history}{a \code{n_epoch} by 2 matrix giving the training and
#'   validation losses at the end of each epoch.} 
#'
#' @seealso \code{\link{predict.seqm}} for the \code{predict} method for \code{seqm} objects.
#' @examples
#' n <- 100
#' seqs <- seq_gen(n)
#' y1 <- sapply(seqs, function(x) "CHECK_A" %in% x)
#' y2 <- sapply(seqs, function(x) log10(length(x)))
#' x <- rnorm(n)
#' mydata <- data.frame(x=x, y1=y1, y2=y2)
#' 
#' index_test <- 91:100
#' index_train <- 1:90
#' 
#' actions <- unique(unlist(seqs))
#' 
#' res1 <- seqm(y1 ~ x, "binary", seqs[index_train], actions=actions, data=mydata[index_train, ], K_emb = 5, K_rnn = 5, valid_split=0.2, n_epoch = 5)
#' predict(res1, new_seqs = seqs[index_test], new_data=mydata[index_test, ])
#' 
#' res1_more <- seqm(y1 ~ x, "binary", seqs[index_train], actions=actions, data=mydata[index_train, ], K_emb = 5, K_rnn = 5, valid_split=0.2, n_hidden=2, K_hidden=c(10,5), n_epoch = 5)
#' predict(res1_more, new_seqs = seqs[index_test], new_data=mydata[index_test, ])
#' 
#' res2 <- seqm(y2 ~ x, "scale", seqs[index_train], actions=actions, data=mydata[index_train, ], K_emb = 5, K_rnn = 5, valid_split=0.2, n_epoch = 5)
#' predict(res2, new_seqs = seqs[index_test], new_data=mydata[index_test, ])
#' 
#' @export
seqm <- function(formula, response_type, seqs, actions = NULL, data, rnn_type = "lstm", K_emb = 20, K_rnn = 20, n_hidden = 0, K_hidden = NULL,
                 valid_split = 0, index_valid = NULL, max_len = max(sapply(seqs, length)),
                 n_epoch = 20, batch_size = 16, optimizer_name = "rmsprop", step_size = 0.001, gpu = TRUE)
{
  n_person <- length(seqs)
  if (is.null(actions)) events <- unique(unlist(seqs))
  else {
    events <- union(unique(unlist(seqs)), actions)
    if (length(events) > length(actions)) warning("Action set is expanded to include all actions in seqs.\n")
  }
  n_event <- length(events)
  
  max_len0 <- max(sapply(seqs, length))
  if (max_len < max_len0) {
    warning("max_len is set as the max length in seqs!\n")
    max_len <- max_len0
  }

  if (!is.null(index_valid)) {
    index_train <- setdiff(1:n_person, index_valid)
    if (valid_split != 0) warning("Both valid_split and index_valid are set. Use index_valid as validation set.\n")
  } else if (valid_split !=0) {
    index_valid <- sample(1:n_person, round(n_person*valid_split))
    index_train <- setdiff(1:n_person, index_valid)
  } else {
    stop("Validation set is empty! Specify either valid_split or index_valid.\n")
  }

  int_seqs <- matrix(0, n_person, max_len)
  
  for (index_seq in 1:n_person) {
    my_seq <- seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:n_l] <- tmp
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # organize covariates
  covariates <- model.matrix(formula, data)[,-1,drop=FALSE]
  response <- model.extract(model.frame(formula, data), "response")
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, K_emb, mask_zero=TRUE)
  
  K_cov <- ncol(covariates)
  cov_inputs <- layer_input(shape=list(K_cov))
  cov_expand <- cov_inputs %>% layer_repeat_vector(n=max_len)
  seq_plus_cov <- layer_concatenate(inputs=c(seq_emb, cov_expand), axis=-1)
  if (rnn_type == "lstm") seq_feature <- seq_plus_cov %>% layer_lstm(units=K_rnn)
  else if (rnn_type == "gru") seq_feature <- seq_plus_cov %>% layer_gru(units=K_rnn)
  
  outputs <- seq_feature
  n_hidden <- min(n_hidden, length(K_hidden))
  
  if (n_hidden > 0) {
    for (index_hidden in 1:n_hidden) outputs <- outputs %>% layer_dense(units=K_hidden[index_hidden], activation='tanh')
  }
  
  if (response_type == "binary") outputs <- outputs %>% layer_dense(units=1, activation='sigmoid')
  else if (response_type == "scale") outputs <- outputs %>% layer_dense(units=1, activation='linear')

  seq_model <- keras_model(c(seq_inputs, cov_inputs), outputs)
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  if (response_type == "binary") seq_model %>% compile(optimizer = optimizer, 
                                                       loss='binary_crossentropy')
  else if (response_type == "scale") seq_model %>% compile(optimizer = optimizer, 
                                                           loss='mean_squared_error')
  best_valid <- Inf
  trace_res <- matrix(0, n_epoch, 2)
  colnames(trace_res) <- c("train", "valid")
  for (index_epoch in 1:n_epoch) {
    model_res <- seq_model %>% fit(list(int_seqs[index_train, ], covariates[index_train,]), response[index_train], 
                                 epochs=1, batch_size=batch_size, verbose=FALSE, 
                                 validation_data=list(list(int_seqs[index_valid,], covariates[index_valid,]), response[index_valid]))
    trace_res[index_epoch, 1] <- model_res$metrics$loss
    trace_res[index_epoch, 2] <- model_res$metrics$val_loss
    if (model_res$metrics$val_loss < best_valid) {
      best_valid <- model_res$metrics$val_loss
      model_save <- serialize_model(seq_model)
    }
  }
  
  k_clear_session()
  
  res <- list(formula = formula, model_fit = model_save, actions = events, max_len = max_len, history = trace_res) 
  class(res) <- "seqm"

  res  
}

#' Predict method for sequence models
#' 
#' Obtains predictions from a fitted sequence model object.
#' 
#' It unserialize object$model to obtain a keras model of class 
#' \code{"keras.engin.training.Model"} and then calls \code{predict} for this class
#' to obtain predictions.
#' 
#' @param object a fitted object of class \code{"seqm"} from \code{seqm}.
#' @param new_seqs a list of action sequences with which to predict.
#' @param new_data a dataframe in which to look for variables with which to predict.
#' @param ... further arguments to be passed to \code{predict.keras.engine.training.Model}.
#' 
#' @return a vector of predictions. If \code{response_type="binary"}, predictions are
#'   probabilities of the response variable being 1.
#' @seealso \code{\link{seqm}} for fitting sequence models.
#' @export
predict.seqm <- function(object, new_seqs, new_data, ...) {
  model <- unserialize_model(object$model)
  max_len <- object$max_len
  events <- object$actions
  ff <- object$formula
  
  n <- length(new_seqs)
  new_int_seqs <- matrix(0, n, max_len)
  
  for (index_seq in 1:n) {
    my_seq <- new_seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    new_int_seqs[index_seq, 1:n_l] <- tmp
  }
  new_covariates <- model.matrix(ff, new_data)[,-1,drop=FALSE]
  
  pred_res <- predict(model, list(new_int_seqs, new_covariates), ...)
  
  k_clear_session()
  
  pred_res
}

