
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

#' Fitting sequence models with numeric responses
#' 
#' \code{seq2scale} is used to fit a neural network model relating 
#' action sequences with numeric responses.
#'
#' The model fitted in this function is the same as that in \code{seq2binary} except
#' that the last fully-connected layer has a linear activation.
#' 
#' @inheritParams seq2feature_seq2seq
#' @inheritParams seq2binary
#' @param response the continuous response variable.
#' @param index_train,index_valid,index_test vectors of indices specifying the training,
#'   validation, and test sets.
#' @return \code{seq2scale} returns a list containing
#'   \item{model}{the trained keras model}
#'   \item{summary}{a 3 by 2 matrix summarizing the prediction performance. 
#'     The rows correspond to the training, validation, and test sets. 
#'     The two columns gives the mean squared error and the R-square.}
#'   \item{trace}{a \code{n_epoch} by 2 matrix giving the trace of training 
#'     and validation losses in the training process.}
#'   \item{pred_train}{the fitted probabilities on the training set.}
#'   \item{pred_valid}{the predicted probabilities for the validation set.}
#'   \item{pred_test}{the predicted probabilities for the test set.}
#'   \item{events}{the set of all possible actions.}
#'   \item{max_len}{the length of padded sequnces.}
#' @family sequence models
#' @seealso \code{\link{seq_predict}} for prediction from a sequence model.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' y <- log10(sapply(seqs, length))
#' index_train <- sample(1:n, round(0.8*n))
#' index_valid <- sample(setdiff(1:n, index_train), round(0.1*n))
#' index_test <- setdiff(1:n, c(index_train, index_valid))
#' res <- seq2scale(seqs, y, index_train = index_train, index_valid = index_valid, index_test = index_test)
#' @export
seq2scale <- function(seqs, response, rnn_type = "lstm", n_hidden = 0, K = 20, K_hidden = NULL, n_epoch=20, batch_size = 16, optimizer_name="rmsprop", step_size=0.001, index_train, index_valid, index_test, gpu=TRUE, return_model=TRUE, model_output = "seq2scale_model.h5", max_len = max(sapply(seqs, length)))
{
  n_person <- length(seqs)
  events <- unique(unlist(seqs))
  n_event <- length(events)
  max_len <- max(sapply(seqs, length))
  
  int_seqs <- matrix(0, n_person, max_len)
  
  for (index_seq in 1:n_person) {
    my_seq <- seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:n_l] <- tmp
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, K, mask_zero=TRUE)
  if (rnn_type == "lstm") seq_feature <- seq_emb %>% layer_lstm(units=K)
  else if (rnn_type == "gru") seq_feature <- seq_emb %>% layer_gru(units=K)
  
  n_hidden <- min(n_hidden, length(K_hidden))
  
  ff_string <- "scale_outputs <- seq_feature %>% "
  if (n_hidden > 0)
  {
    for (index_hidden in 1:n_hidden) ff_string <- paste(ff_string, "layer_dense(units=", K_hidden[index_hidden], ", activation='tanh') %>% ", sep="")
  }
  ff_string <- paste(ff_string, "layer_dense(units=1, activation='linear')", sep="")
  eval(parse(text=ff_string))
  
  seq2scale_model <- keras_model(seq_inputs, scale_outputs)
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  seq2scale_model %>% compile(optimizer = optimizer, loss='mean_squared_error')
  
  model_res <- seq2scale_model %>% fit(int_seqs[index_train, ], response[index_train], epochs=n_epoch, batch_size=batch_size, verbose=FALSE, validation_data=list(int_seqs[index_valid,], response[index_valid]), callbacks=list(callback_model_checkpoint(model_output, monitor="val_loss", save_best_only = TRUE)))
  
  trace_res <- cbind(model_res$metrics$loss, model_res$metrics$val_loss)
  colnames(trace_res) <- c("train", "validation")
  
  seq2scale_model <- load_model_hdf5(model_output)
  
  pred_train <- predict(seq2scale_model, int_seqs[index_train,], batch_size = batch_size)
  pred_valid <- predict(seq2scale_model, int_seqs[index_valid,], batch_size = batch_size)
  pred_test <- predict(seq2scale_model, int_seqs[index_test,], batch_size = batch_size)
  
  summary_res <- matrix(0, 3, 2)
  summary_res[1,1] <- mean((pred_train - response[index_train])^2)
  summary_res[2,1] <- mean((pred_valid - response[index_valid])^2)
  summary_res[3,1] <- mean((pred_test - response[index_test])^2)
  summary_res[1,2] <- cor(pred_train, response[index_train])^2
  summary_res[2,2] <- cor(pred_valid, response[index_valid])^2
  summary_res[3,2] <- cor(pred_test, response[index_test])^2
  colnames(summary_res) <- c("mse", "rsq")
  rownames(summary_res) <- c("train", "valid", "test")
  
  if (return_model) {
    res <- list(model = seq2scale_model, summary = summary_res, trace = trace_res, pred = pred_test, events=events, max_len = max_len)
  } else {
      res <- list(summary = summary_res, trace = trace_res, pred = pred_test)
      k_clear_session()
  }
  
  res  
    
}


#' Fitting sequence models with binary responses
#' 
#' \code{seq2binary} is used to fit a neural network model relating 
#' action sequences with binary responses.
#' 
#' The model consists of an embedding layer, a recurrent layer and one or more 
#' fully connected layers. The embedding layer takes a action sequence and output a
#' sequences of \code{K} dimensional numeric vectors to the recurrent layer. The last
#' output of the recurrent layer is used as the input of the subsequent fully connected
#' layers. The last layer uses the sigmoid activation to produce a probability of the
#' response being positive. The dimension of the output of other fully connected layers
#' (if any) is specified by \code{K_hidden}. 
#' 
#' The action sequences are re-coded into integer sequences and padded with zeros so that
#' every sequene is of length \code{max_len}. If the provided \code{max_len} is smaller
#' than the length of the longest sequence in \code{seqs}, it will be overridden.
#' 
#' @inheritParams seq2feature_seq2seq
#' @param response the binary response variable.
#' @param n_hidden the number of hidden fully-connected layers.
#' @param K the latent dimension of the embedding layer and the recurrent layer. 
#' @param K_hidden a vector of length \code{n_hidden} specifying
#'   the number of nodes in each hidden layer.
#' @param n_epoch the number of training epochs.
#' @param batch_size the batch size used in training.
#' @param index_train,index_valid,index_test vectors of indices specifying the training,
#'   validation, and test sets.
#' @param return_model logical. If TRUE, the trained keras model is returned.
#' @param model_output a character string specifying the filename for saving the trained keras model.
#' @param max_len the maximum length of sequences.
#' @return \code{seq2binary} returns a list containing
#'   \item{model}{the trained keras model.}
#'   \item{summary}{a vector of length 3 summarizing the prediction accuracy on the training, validation, and test sets.}
#'   \item{trace}{a \code{n_epoch} by 2 matrix giving the trace of training and validation losses in the training process.}
#'   \item{pred_train}{the fitted probabilities on the training set.}
#'   \item{pred_valid}{the predicted probabilities for the validation set.}
#'   \item{pred_test}{the predicted probabilities for the test set.}
#'   \item{events}{the set of all possible actions.}
#'   \item{max_len}{the length of padded sequnces.}
#' @family sequence models
#' @seealso \code{\link{seq_predict}} for prediction from a sequence model.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' y <- sapply(seqs, function(x) "CHECK_A" %in% x)
#' index_train <- sample(1:n, round(0.8*n))
#' index_valid <- sample(setdiff(1:n, index_train), round(0.1*n))
#' index_test <- setdiff(1:n, c(index_train, index_valid))
#' res <- seq2binary(seqs, y, index_train = index_train, index_valid = index_valid, index_test = index_test)
#' @export
seq2binary <- function(seqs, response, rnn_type = "lstm", n_hidden = 0, K = 20, K_hidden = NULL, n_epoch=20, batch_size = 16, optimizer_name="rmsprop", step_size=0.001, index_train, index_valid, index_test, gpu = TRUE, return_model = TRUE, model_output = "seq2binary_model.h5", max_len = max(sapply(seqs, length)))
{
  n_person <- length(seqs)
  events <- unique(unlist(seqs))
  n_event <- length(events)
  max_len0 <- max(sapply(seqs, length))
  if (max_len < max_len0) {
    warning("max_len is set as the max length in seqs!\n")
    max_len <- max_len0
  }
  int_seqs <- matrix(0, n_person, max_len)
  
  for (index_seq in 1:n_person) {
    my_seq <- seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:n_l] <- tmp
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, K, mask_zero=TRUE)
  if (rnn_type == "lstm") seq_feature <- seq_emb %>% layer_lstm(units=K)
  else if (rnn_type == "gru") seq_feature <- seq_emb %>% layer_gru(units=K)
  
  n_hidden <- min(n_hidden, length(K_hidden))
  
  ff_string <- "prob_outputs <- seq_feature %>% "
  if (n_hidden > 0)
  {
    for (index_hidden in 1:n_hidden) ff_string <- paste(ff_string, "layer_dense(units=", K_hidden[index_hidden], ", activation='tanh') %>% ", sep="")
  }
  ff_string <- paste(ff_string, "layer_dense(units=1, activation='sigmoid')", sep="")
  eval(parse(text=ff_string))
  seq2binary_model <- keras_model(seq_inputs, prob_outputs)
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  seq2binary_model %>% compile(optimizer = optimizer, loss='binary_crossentropy', metrics='accuracy')
  
  model_res <- seq2binary_model %>% fit(int_seqs[index_train, ], response[index_train], epochs=n_epoch, batch_size=batch_size, verbose=FALSE, metrics='accuracy', validation_data=list(int_seqs[index_valid,], response[index_valid]), callbacks=list(callback_model_checkpoint(model_output, monitor="val_acc", save_best_only = TRUE)))
  
  trace_res <- cbind(model_res$metrics$acc, model_res$metrics$val_acc)
  colnames(trace_res) <- c("train", "validation")
  
  seq2binary_model <- load_model_hdf5(model_output)
  
  pred_train <- predict(seq2binary_model, int_seqs[index_train,], batch_size = batch_size)
  pred_valid <- predict(seq2binary_model, int_seqs[index_valid,], batch_size = batch_size)
  pred_test <- predict(seq2binary_model, int_seqs[index_test,], batch_size = batch_size)
  
  summary_res <- rep(0, 3)
  summary_res[1] <- mean(as.numeric(pred_train > 0.5) == response[index_train])
  summary_res[2] <- mean(as.numeric(pred_valid > 0.5) == response[index_valid])
  summary_res[3] <- mean(as.numeric(pred_test > 0.5) == response[index_test])
  
  names(summary_res) <- c("train", "valid", "test")
  
  if (return_model) {
    res <- list(model = seq2binary_model, summary = summary_res, trace = trace_res, pred = pred_test, events = events, max_len = max_len)
  } else {
      res <- list(summary = summary_res, trace = trace_res, pred = pred_test)
      k_clear_session()
  }
  
  res  
} 

#' Sequence model predictions
#'
#' \code{seq_predict} is a function for predictions from the models fitted by 
#'   \code{seq2binary} or \code{seq2scale}.
#'
#' @param fit_res an object returned by \code{seq2binary} or \code{seq2scale}
#'   with \code{return_model=TRUE}.
#' @param new_seqs action sequences with which to predict.
#' @param batch_size batch size for computing predicted values
#' 
#' @return \code{seq_predict} produces a vector of predicted values.
#' @seealso The model fitting functions \code{\link{seq2binary}} and \code{\link{seq2scale}}.
#' @examples
#' n <- 100
#' seqs <- seq_gen(n)
#' y1 <- sapply(seqs, function(x) "CHECK_A" %in% x)
#' y2 <- log10(sapply(seqs, length))
#' index_train <- 1:70
#' index_valid <- 71:85
#' index_test <- 86:100
#' res1 <- seq2binary(seqs, y1, index_train = index_train, index_valid = index_valid, index_test = index_test)
#' res2 <- seq2scale(seqs, y2, index_train = index_train, index_valid = index_valid, index_test = index_test)
#' res1_pred <- seq_predict(res1, seqs[index_test])
#' res2_pred <- seq_predict(res2, seqs[index_test])
#' all.equal(res1$pred, res1_pred) 
#' all.equal(res2$pred, res2_pred)
#'
#' @export
seq_predict <- function(fit_res, new_seqs, batch_size=16) {
  model = fit_res$model
  if (is.null(model)) stop("No model is available for prediction! Set return_model=TRUE in seq2binary or seq2scale.\n")
  max_len <- fit_res$max_len
  events <- fit_res$events

  n <- length(new_seqs)
  new_int_seqs <- matrix(0, n, max_len)
  
  for (index_seq in 1:n) {
    my_seq <- new_seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    new_int_seqs[index_seq, 1:n_l] <- tmp
  }

  predict(model, new_int_seqs, batch_size = batch_size)
  
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
#' @param response the binary response variable.
#' @param response_type "binary" or "scale".
#' @param n_hidden the number of hidden fully-connected layers.
#' @param K the latent dimension of the embedding layer and the recurrent layer.
#' @param K_hidden a vector of length \code{n_hidden} specifying the number of
#'   nodes in each hidden layer.
#' @param n_epoch the number of training epochs.
#' @param batch_size the batch size used in training.
#' @param index_train,index_valid,index_test vectors of indices specifying the
#'   training, validation, and test sets.
#' @param return_model logical. If TRUE, the trained keras model is returned.
#' @param model_output a character string specifying the filename for saving the
#'   trained keras model.
#' @param max_len the maximum length of sequences.
#' @return \code{seq2binary} returns a list containing \item{model}{the trained
#'   keras model.} \item{summary}{a vector of length 3 summarizing the
#'   prediction accuracy on the training, validation, and test sets.}
#'   \item{trace}{a \code{n_epoch} by 2 matrix giving the trace of training and
#'   validation losses in the training process.} \item{pred_train}{the fitted
#'   probabilities on the training set.} \item{pred_valid}{the predicted
#'   probabilities for the validation set.} \item{pred_test}{the predicted
#'   probabilities for the test set.} \item{events}{the set of all possible
#'   actions.} \item{max_len}{the length of padded sequnces.}
#' @family sequence models
#' @seealso \code{\link{seq_predict}} for prediction from a sequence model.
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' y <- sapply(seqs, function(x) "CHECK_A" %in% x)
#' index_train <- sample(1:n, round(0.8*n))
#' index_valid <- sample(setdiff(1:n, index_train), round(0.1*n))
#' index_test <- setdiff(1:n, c(index_train, index_valid))
#' res <- seq2binary(seqs, y, index_train = index_train, index_valid = index_valid, index_test = index_test)
#' @export
seqm <- function(seqs, response, response_type = "scale", rnn_type = "lstm", K = 20, n_hidden = 0, K_hidden = NULL,
                 index_train=NULL, index_valid=NULL, max_len = max(sapply(seqs, length)),
                 n_epoch=20, batch_size = 16, optimizer_name="rmsprop", step_size=0.001, gpu = TRUE, model_output = "seq_model.h5")
{
  n_person <- length(seqs)
  events <- unique(unlist(seqs))
  n_event <- length(events)
  max_len0 <- max(sapply(seqs, length))
  if (max_len < max_len0) {
    warning("max_len is set as the max length in seqs!\n")
    max_len <- max_len0
  }
  int_seqs <- matrix(0, n_person, max_len)
  
  for (index_seq in 1:n_person) {
    my_seq <- seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:n_l] <- tmp
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, K, mask_zero=TRUE)
  if (rnn_type == "lstm") seq_feature <- seq_emb %>% layer_lstm(units=K)
  else if (rnn_type == "gru") seq_feature <- seq_emb %>% layer_gru(units=K)
  
  n_hidden <- min(n_hidden, length(K_hidden))
  
  ff_string <- "prob_outputs <- seq_feature %>% "
  if (n_hidden > 0)
  {
    for (index_hidden in 1:n_hidden) ff_string <- paste(ff_string, "layer_dense(units=", K_hidden[index_hidden], ", activation='tanh') %>% ", sep="")
  }
  if (response_type == "binary") last_act_fun <- "'sigmoid'"
  else if (response_type == "scale") last_act_fun <- "'linear'"
  
  ff_string <- paste(ff_string, "layer_dense(units=1, activation=", last_act_fun, ")", sep="")
  eval(parse(text=ff_string))
  seq_model <- keras_model(seq_inputs, prob_outputs)
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  if (response_type == "binary") seq_model %>% compile(optimizer = optimizer, 
                                                       loss='binary_crossentropy')
  else if (response_type == "scale") seq_model %>% compile(optimizer = optimizer, 
                                                           loss='mean_squared_error')
  
  model_res <- seq_model %>% fit(int_seqs[index_train, ], response[index_train], 
                                 epochs=n_epoch, batch_size=batch_size, verbose=FALSE, 
                                 validation_data=list(int_seqs[index_valid,], response[index_valid]), 
                                 callbacks=list(callback_model_checkpoint(model_output, monitor="val_loss", save_best_only = TRUE)))
  
  trace_res <- cbind(model_res$metrics$loss, model_res$metrics$val_loss)
  colnames(trace_res) <- c("train", "validation")
  
  seq_model <- load_model_hdf5(model_output)
  
  res <- list(model = seq_model, actions = events, max_len = max_len, history = trace_res) 
  
  res  
}