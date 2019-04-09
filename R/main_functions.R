#' @useDynLib ProcData
#' @importFrom Rcpp sourceCpp
#' @import keras
#' @import stats
#' @import utils
NULL

#' Extract features from action sequences by multidimensional scaling
#' 
#' @param seqs a list or a square matrix. If a list, each element is an action sequence in the form of a vector of actions. Alternatively, their dissimilarity matrix can be provided
#' @param K the number of features to be extracted.
#' @param method the dissimilarity measure to be used to calculate distance matrix. Currently available methods: \code{oss}. It implements the order-based sequence similarity in Gomez-Alonso and Valls (2008).
#' @param max_epoch the maximum number of epochs in stochastic gradient descent.
#' @param step_size the step size to be used in stochastic gradient descent
#' @param pca logical. If TRUE, the principal components of features are returned.
#' @param tot the tolerance of accuracy
#' @param return_dist logical. If the dissimilarity matrix should be returned. Default is \code{FALSE}.
#' @return a list
#'   \item{theta}{an \code{n} by \code{K} matrix of \code{K} raw features or principal features for \code{n} sequences.}
#'   \item{loss}{discrepancy between estimated and true dissimilarity}
#'   \item{dist_mat}{dissimilary matrix}
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

#' Choose the number of MDS features by cross validation
#' 
#' @param seqs a list of \code{n} action sequences or their dissimilarity matrix.
#' @param K_cand candidates of number of features.
#' @param method the dissimilarity measure to be used to calculate distance matrix. Currently available methods: \code{oss}. It implements the order-based sequence similarity in Gomez-Alonso and Valls (2008).
#' @param max_epoch the maximum number of epochs in stochastic gradient descent.
#' @param n_fold number of folds for cross-validation
#' @param step_size the step size to be used in stochastic gradient descent
#' @param tot the tolerance of accuracy
#' @param return_dist logical. If dissimilarity matrix should be returned. Default is \code{FALSE}.
#' @return a list
#'   \item{K}{number of features with smallest cross-validation loss}
#'   \item{K_cand}{candidates of number of features}
#'   \item{cv_loss}{cross-validation loss for each candidate K}
#'   \item{dist_mat}{dissimilarity matrix of sequences}
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


#' Extract features from action sequences by sequence-to-sequence autoencoder
#' 
#' @param seqs a list of \code{n} action sequences. Each element is an action sequence in the form of a vector of actions.
#' @param K the number of features to be extracted.
#' @param rnn_type the type of recurrent neural network to be used for modeling action sequences. \code{"lstm"} for long-short term memory. \code{"gru"} for gated recurrent unit.
#' @param n_epoch the number of epochs to be run when minimizing the loss function
#' @param method the method to be used to compute features from the output of an RNN; With \code{method="last"}, the features are the last output of the encoder RNN; with \code{method="avg"}, the features are the average of the outputs of the encoder RNN. keep the output of the last step; "average": average of the outputs of all time steps
#' @param step_size (baseline) learning rate of optimizer
#' @param optimizer_name a character string specifying the optimizer to be used for training: sgd, rmsprop, adadelta, adam
#' @param samples_train specify train set;
#' @param samples_valid specify validation set; 
#' @param samples_test specify test set;
#' @param pca logical. If TRUE, the principal components of features are returned.
#' @param gpu logical. If TRUE, use gpu (if available) for training
#' @param verbose logical. If TRUE, training progress is printed.
#' @param return_theta logical. If TRUE, extracted features are returned.
#' @return a list
#'   \item{features}{an \code{n} by \code{K} matrix which gives \code{K} raw features or principal features for \code{n} sequences.}
#'   \item{train_loss}{a vector of length \code{n_epoch}; training loss trace}
#'   \item{valid_loss}{a vector of length \code{n_epoch}; validation loss trace}
#'   \item{test_loss}{a vector of length \code{n_epoch}; test loss trace}
#' @examples 
#' n <- 50
#' seqs <- seq_gen(n)
#' seq2seq_res <- seq2feature_seq2seq(seqs, 5, n_epoch=10, samples_train=1:40, samples_valid=41:50)
#' features <- seq2seq_res$theta
#' plot(seq2seq_res$train_loss, col="blue", type="l")
#' lines(seq2seq_res$valid_loss, col="red")
#' @export
seq2feature_seq2seq <- function(seqs, K, rnn_type="lstm", n_epoch=50, method="last", step_size=0.0001, optimizer_name="adam", samples_train, samples_valid, samples_test=NULL, pca=TRUE, gpu=FALSE, verbose=TRUE, return_theta=TRUE) {
  
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
      layer_lstm(units=K, return_sequences=TRUE) %>% 
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
  
  # if (length(valid_set) > 1) samples_valid <- valid_set
  # else 
  # {
  #   if (valid_set < 1) n_valid <- floor(n_seq*valid_set)
  #   else n_valid <- floor(valid_set)
  #   samples_valid <- sample(1:n_seq, n_valid)
  # }
  # samples_train <- setdiff(1:n_seq, samples_valid)
  
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

#' Choose the number of seq2seq features by cross validation
#' 
#' @param seqs a list of \code{n} action sequences. Each element is an action sequence in the form of a vector of actions.
#' @param K_cand candidates of number of features.
#' @param rnn_type the type of recurrent neural network to be used for modeling action sequences. \code{"lstm"} for long-short term memory. \code{"gru"} for gated recurrent unit.
#' @param n_epoch the number of epochs to be run when minimizing the loss function
#' @param method the method to be used to compute features from the output of an RNN; With \code{method="last"}, the features are the last output of the encoder RNN; with \code{method="avg"}, the features are the average of the outputs of the encoder RNN. keep the output of the last step; "average": average of the outputs of all time steps
#' @param step_size (baseline) learning rate of optimizer
#' @param optimizer_name a character string specifying the optimizer to be used for training: sgd, rmsprop, adadelta, adam
#' @param n_fold number of folds for cross-validation
#' @param valid_prop proportion of validation samples in each fold.
#' @param gpu logical. If TRUE, use gpu (if available) for training
#' @param verbose logical. If TRUE, training progress is printed.
#' @return a list
#'   \item{K}{number of features with smallest cross-validation loss}
#'   \item{K_cand}{candidates of number of features}
#'   \item{cv_loss}{cross-validation loss for each candidate K}
#' @examples 
#' n <- 50
#' seqs <- seq_gen(n)
#' K_res <- chooseK_seq2seq(seqs, 5:10, n_epoch=10, n_fold=2, valid_prop=0.2)
#' seq2seq_res <- seq2feature_seq2seq(seqs, K_res$K, n_epoch=10, samples_train=1:40, samples_valid=41:50)
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

#' Predict continuous variable from action sequences
#' 
#' @param seqs a list of \code{n} action sequences. Each element is an action sequence in the form of a vector of actions.
#' @param response the continuous response variable
#' @param model_output a character string giving the name of the file to store keras model.
#' @param rnn_type the type of recurrent neural network to be used for modeling action sequences. \code{"lstm"} for long-short term memory. \code{"gru"} for gated recurrent unit.
#' @param n_hidden the number of hidden layers in the feed-forward neural network.
#' @param K the latent dimension in the recurrent neural network.
#' @param K_hidden a vector of length \code{n_hidden} specifying the number of nodes in each hidden layers.
#' @param n_epoch the number of epochs to be run in training.
#' @param batch_size batch size in training.
#' @param optimizer a character string specifies the optimizer to be used for training.
#' @param index_train a vector of indices specifying the training set
#' @param index_valid a vector of indices specifying the validation set
#' @param index_test a vector of indices specifying the test set
#' @param gpu logical. If TRUE, use gpu (if available) for training
#' @return a list
#'   \item{summary}{a 3 by 2 matrix summarizing the prediction performance. Rows of the matrix indicate training set, validation set and test set. Columns of the matrix indicate mean squared error and R-square.}
#'   \item{trace}{a \code{n_epoch} by 2 matrix giving the trace of model training. The two columns gives the trace for the training set and the validation set, respectively.}
#'   \item{pred_train}{fitted values for the training set.}
#'   \item{pred_valid}{predicted values for the validation set.}
#'   \item{pred_test}{predicted values for the test set.}
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' y <- log10(sapply(seqs, length))
#' index_train <- sample(1:n, round(0.8*n))
#' index_valid <- sample(setdiff(1:n, index_train), round(0.1*n))
#' index_test <- setdiff(1:n, c(index_train, index_valid))
#' res <- seq2scale(seqs, y, index_train = index_train, index_valid = index_valid, index_test = index_valid)
#' @export
seq2scale <- function(seqs, response, model_output = "seq2scale_model.h5", rnn_type = "lstm", n_hidden = 0, K = 20, K_hidden = NULL, n_epoch=20, batch_size = 16, optimizer="rmsprop", index_train, index_valid, index_test, gpu=TRUE)
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
  
  k_clear_session()
  
  return(list(summary = summary_res, trace = trace_res, pred = pred_test))
    
}


#' Predict binary variable from action sequences
#' 
#' @param seqs a list of \code{n} action sequences. Each element is an action sequence in the form of a vector of actions.
#' @param response the binary response variable
#' @param model_output a character string giving the name of the file to store keras model.
#' @param rnn_type the type of recurrent neural network to be used for modeling action sequences. \code{"lstm"} for long-short term memory. \code{"gru"} for gated recurrent unit.
#' @param n_hidden the number of hidden layers in the feed-forward neural network.
#' @param K the latent dimension in the recurrent neural network. 
#' @param K_hidden a vector of length \code{n_hidden} specifying the number of nodes in each hidden layers.
#' @param n_epoch the number of epochs to be run in training.
#' @param batch_size batch size in training.
#' @param optimizer a character string specifies the optimizer to be used for training.
#' @param index_train a vector of indices specifying the training set
#' @param index_valid a vector of indices specifying the validation set
#' @param index_test a vector of indices specifying the test set
#' @param gpu logical. If TRUE, use gpu (if available) for training
#' @return a list
#'   \item{summary}{a vector of length 3 summarizing the prediction accuracy in the training, validation and test sets.}
#'   \item{trace}{a \code{n_epoch} by 2 matrix giving the trace of model training. The two columns gives the trace for the training set and the validation set, respectively.}
#'   \item{pred_train}{fitted probabilities for the training set.}
#'   \item{pred_valid}{predicted probabilities for the validation set.}
#'   \item{pred_test}{predicted probabilities for the test set.}
#' @examples
#' n <- 50
#' seqs <- seq_gen(n)
#' y <- sapply(seqs, function(x) "CHECK_A" %in% x)
#' index_train <- sample(1:n, round(0.8*n))
#' index_valid <- sample(setdiff(1:n, index_train), round(0.1*n))
#' index_test <- setdiff(1:n, c(index_train, index_valid))
#' res <- seq2binary(seqs, y, index_train = index_train, index_valid = index_valid, index_test = index_valid)
#' @export
seq2binary <- function(seqs, response, model_output = "seq2binary_model.h5", rnn_type = "lstm", n_hidden = 0, K = 20, K_hidden = NULL, n_epoch=20, batch_size = 16, optimizer="rmsprop", index_train, index_valid, index_test, gpu = TRUE)
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
  
  ff_string <- "prob_outputs <- seq_feature %>% "
  if (n_hidden > 0)
  {
    for (index_hidden in 1:n_hidden) ff_string <- paste(ff_string, "layer_dense(units=", K_hidden[index_hidden], ", activation='tanh') %>% ", sep="")
  }
  ff_string <- paste(ff_string, "layer_dense(units=1, activation='sigmoid')", sep="")
  eval(parse(text=ff_string))
  seq2binary_model <- keras_model(seq_inputs, prob_outputs)
  
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
  
  k_clear_session()
  
  return(list(summary = summary_res, trace = trace_res, pred = pred_test))
} 
  
  
