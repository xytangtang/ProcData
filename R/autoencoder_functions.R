#' Feature Extraction by action sequence autoencoder
#'
#' \code{aseq2feature_seq2seq} extract features from action sequences by action
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
#' @inheritParams seq2feature_seq2seq
#' @param aseqs a list of \code{n} action sequences. Each element is an action
#'   sequence in the form of a vector of actions.
#' @return \code{aseq2feature_seq2seq} returns a list containing
#'   \item{theta}{a matrix containing \code{K} features or principal features. Each column is a feature.}
#'   \item{train_loss}{a vector of length \code{n_epoch} recording the trace of training losses.}
#'   \item{valid_loss}{a vector of length \code{n_epoch} recording the trace of validation losses.}
#'   \item{test_loss}{a vector of length \code{n_epoch} recording the trace of test losses. Exists only if \code{samples_test} is not \code{NULL}.}
#' @seealso \code{\link{chooseK_seq2seq}} for choosing \code{K} through cross-validation.
#' @examples
#' \donttest{
#' n <- 50
#' seqs <- seq_gen(n)
#' seq2seq_res <- aseq2feature_seq2seq(seqs$action_seqs, 5, rnn_type="lstm", n_epoch=5, 
#'                                    samples_train=1:40, samples_valid=41:50)
#' features <- seq2seq_res$theta
#' plot(seq2seq_res$train_loss, col="blue", type="l")
#' lines(seq2seq_res$valid_loss, col="red")
#' }
#' @export
aseq2feature_seq2seq <- function(aseqs, K, rnn_type="lstm", n_epoch=50, method="last", 
                                step_size=0.0001, optimizer_name="adam", 
                                samples_train, samples_valid, samples_test=NULL, 
                                pca=TRUE, gpu=FALSE, parallel=FALSE, seed=12345,
                                verbose=TRUE, return_theta=TRUE) {
  use_session_with_seed(seed, disable_gpu = !gpu, disable_parallel_cpu = !parallel)
  #tensorflow::use_compat("v1")
  if (!(rnn_type %in% c("lstm", "gru"))) 
    stop("Invalid rnn_type! Available options: lstm, gru.\n")
  if (!(method %in% c("last", "avg"))) 
    stop("Invalid method! Available options: last, avg.\n")
  if (!(optimizer_name %in% c("sgd", "adam", "rmsprop", "adadelta"))) 
    stop("Invalid optimizer! Available options: sgd, adam, rmsprop, adadelta.\n")
  n_seq <- length(aseqs)
  
  events <- unique(unlist(aseqs))
  n_event <- length(events)
  
  # convert action sequence to one-hot vectors
  int_seqs <- list()
  target_seqs <- list()
  
  for (index_seq in 1:n_seq) {
    my_seq <- aseqs[[index_seq]]
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
    train_loss[index_epoch] <- train_loss[index_epoch] / length(samples_train)
    valid_loss[index_epoch] <- valid_loss[index_epoch] / length(samples_valid)
    if (!is.null(samples_test)) 
      test_loss[index_epoch] <- test_loss[index_epoch] / length(samples_test)
    
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


#' Feature Extraction by time sequence autoencoder
#'
#' \code{tseq2feature_seq2seq} extract features from timestamps of action sequences by a 
#' sequence autoencoder.
#'
#' This function trains a sequence-to-sequence autoencoder using keras. The encoder
#' of the autoencoder consists of a recurrent neural network.
#' The decoder consists of another recurrent neural network and a fully connected layer
#' with ReLU activation. The outputs of the encoder are the extracted features.
#' 
#' The output of the encoder is a function of the encoder recurrent neural network.
#' It is the last latent state of the encoder recurrent neural network if \code{method="last"}
#' and the average of the encoder recurrent neural network latent states if \code{method="avg"}.
#' 
#' 
#' @family feature extraction methods
#' @inheritParams seq2feature_seq2seq
#' @param tseqs a list of \code{n} timestamp sequences. Each element is a numeric
#'   sequence in the form of a vector of timestamps associated with actions, with
#'   the timestamp of the first event (e.g., "start") of 0.
#' @return \code{tseq2feature_seq2seq} returns a list containing
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
#' tseqs <- cc_data$seqs$time_seqs[samples]
#' time_seq2seq_res <- tseq2feature_seq2seq(tseqs, 5, rnn_type="lstm", n_epoch=5, 
#'                                    samples_train=1:40, samples_valid=41:50)
#' features <- time_seq2seq_res$theta
#' plot(time_seq2seq_res$train_loss, col="blue", type="l",
#'      ylim = range(c(time_seq2seq_res$train_loss, time_seq2seq_res$valid_loss)))
#' lines(time_seq2seq_res$valid_loss, col="red", type = 'l')
#' }
#' @export
tseq2feature_seq2seq <- function(tseqs, K, cumulative = FALSE, log = TRUE, rnn_type="lstm",
                                 n_epoch=50, method="last", step_size=0.0001, 
                                 optimizer_name="rmsprop", samples_train, samples_valid, 
                                 samples_test=NULL, pca=TRUE, gpu=FALSE, parallel=FALSE, 
                                 seed=12345, verbose=TRUE, return_theta=TRUE) {
  
  use_session_with_seed(seed, disable_gpu = !gpu, disable_parallel_cpu = !parallel)
  if (!(rnn_type %in% c("lstm", "gru"))) 
    stop("Invalid rnn_type! Available options: lstm, gru.\n")
  if (!(method %in% c("last", "avg"))) 
    stop("Invalid method! Available options: last, avg.\n")
  if (!(optimizer_name %in% c("sgd", "adam", "rmsprop", "adadelta"))) 
    stop("Invalid optimizer! Available options: sgd, adam, rmsprop, adadelta.\n")
  
  n_seq <- length(tseqs)
  
  # Get input/output sequences
  if(cumulative){
    if(log){
      times_input <- sapply(tseqs, function(x) array(log10(x+1), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(log10(x+1), dim = c(1,length(tseq2interval(x)),1)))
    }else{
      times_input <- sapply(tseqs, function(x) array((x), dim = c(1,length(x))))
      times_output <- sapply(tseqs, function(x) array((x), dim = c(1,length(x),1)))
    }
  }else{
    if(log){
      times_input <- sapply(tseqs, function(x) array(log10(tseq2interval(x)+1), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(log10(tseq2interval(x)+1), dim = c(1,length(tseq2interval(x)),1)))
    }else{
      times_input <- sapply(tseqs, function(x) array(tseq2interval(x), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(tseq2interval(x), dim = c(1,length(tseq2interval(x)),1)))
    }
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # Define keras model
  # Define an input sequence and process it.
  encoder_inputs <- layer_input(shape=list(NULL),dtype = 'float32')
  expand <- function(x){
    return(k_expand_dims(x, axis = 3))
  }
  
  if (rnn_type == "lstm")
  {
    if (method=="last") {
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand) %>%
        layer_lstm(units = K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[3]]
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[3]], axis=2)
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand) %>%
        layer_lstm(units = K, return_sequences=TRUE)
      
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
      layer_dense(1, activation='relu')
  } else if (rnn_type == "gru") {
    if (method=="last") {
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand) %>%
        layer_gru(units = K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[2]]
      
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[2]], axis=2)
        
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand) %>%
        layer_gru(units = K, return_sequences=TRUE)
      
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
      layer_dense(1, activation='relu')
  }
  
  autoencoder_model <- keras_model(inputs = encoder_inputs, outputs = decoder_outputs)
  encoder_model <- keras_model(inputs = encoder_inputs, outputs = encoder_outputs)
  
  # Run training
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  autoencoder_model %>% compile(optimizer=optimizer, loss='mean_squared_error')
  
  best_valid_loss <- Inf
  valid_loss <- rep(0, n_epoch)
  train_loss <- rep(0, n_epoch)
  if (!is.null(samples_test)) test_loss <- rep(0, n_epoch)
  
  theta <- matrix(0, n_seq, K)
  
  for (index_epoch in 1:n_epoch) {
    if (verbose) cat("Epoch", index_epoch, "\n")
    samples_train <- sample(samples_train)
    
    for (index_seq in samples_train) {
      autoencoder_model %>% train_on_batch(times_input[[index_seq]],times_output[[index_seq]])
    }
    for (index_seq in samples_train) {
      train_loss[index_epoch] <- train_loss[index_epoch] + evaluate(autoencoder_model, times_input[[index_seq]], times_output[[index_seq]], verbose=0)
    }
    for (index_seq in samples_valid) {
      valid_loss[index_epoch] <- valid_loss[index_epoch] + evaluate(autoencoder_model, times_input[[index_seq]], times_output[[index_seq]], verbose=0)
    }
    for (index_seq in samples_test) {
      test_loss[index_epoch] <- test_loss[index_epoch] + evaluate(autoencoder_model, times_input[[index_seq]], times_output[[index_seq]], verbose=0)
    }
    train_loss[index_epoch] <- train_loss[index_epoch] / length(samples_train)
    valid_loss[index_epoch] <- valid_loss[index_epoch] / length(samples_valid)
    if (!is.null(samples_test)) 
      test_loss[index_epoch] <- test_loss[index_epoch] / length(samples_test)
    
    if (valid_loss[index_epoch] < best_valid_loss) {
      best_valid_loss <- valid_loss[index_epoch]
      if (return_theta) {
        for(index_seq in 1:n_seq) theta[index_seq,] <- as.numeric(predict(encoder_model, times_input[[index_seq]]))
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

#' Feature Extraction by action and time sequence autoencoder
#'
#' \code{atseq2feature_seq2seq} extract features from action and timestamp sequences by a 
#' sequence autoencoder.
#'
#' This function trains a sequence-to-sequence autoencoder using keras. The encoder
#' of the autoencoder consists of a recurrent neural network.
#' The decoder consists of another recurrent neural network followed by a fully connected layer
#' with softmax activation for actions and another fully connected layer with ReLU activation 
#' for times. The outputs of the encoder are the extracted features.
#' 
#' The output of the encoder is a function of the encoder recurrent neural network.
#' It is the last latent state of the encoder recurrent neural network if \code{method="last"}
#' and the average of the encoder recurrent neural network latent states if \code{method="avg"}.
#' 
#' 
#' @family feature extraction methods
#' @inheritParams seq2feature_seq2seq
#' @param atseqs a list of two elements, first element is the list of \code{n} action sequences, Each element 
#'   is an action sequence in the form of a vector of actions. The second element is the list of \code{n} 
#'   timestamp sequences corresponding to the action sequences. Each element is a numeric sequence in the form 
#'   of a vector of timestamps associated with actions, with the timestamp of the first event (e.g., "start") of 0.
#' @return \code{tseq2feature_seq2seq} returns a list containing
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
#' atseqs <- sub_seqs(cc_data$seqs, samples)
#' action_and_time_seq2seq_res <- atseq2feature_seq2seq(atseqs, 5, rnn_type="lstm", n_epoch=5, 
#'                                    samples_train=1:40, samples_valid=41:50)
#' features <- action_and_time_seq2seq_res$theta
#' plot(action_and_time_seq2seq_res$train_loss, col="blue", type="l",
#'      ylim = range(c(action_and_time_seq2seq_res$train_loss, 
#'                     action_and_time_seq2seq_res$valid_loss)))
#' lines(action_and_time_seq2seq_res$valid_loss, col="red", type = 'l')
#' }
#' @export
atseq2feature_seq2seq <- function(atseqs, K, weights = c(1, .5), cumulative = FALSE, 
                                  log = TRUE, rnn_type="lstm", n_epoch=50, method="last", 
                                  step_size=0.0001, optimizer_name="rmsprop", 
                                  samples_train, samples_valid, samples_test=NULL, 
                                  pca=TRUE, gpu=FALSE, parallel=FALSE, seed=12345,
                                  verbose=TRUE, return_theta=TRUE) {
  
  use_session_with_seed(seed, disable_gpu = !gpu, disable_parallel_cpu = !parallel)
  if (!(rnn_type %in% c("lstm", "gru"))) 
    stop("Invalid rnn_type! Available options: lstm, gru.\n")
  if (!(method %in% c("last", "avg"))) 
    stop("Invalid method! Available options: last, avg.\n")
  if (!(optimizer_name %in% c("sgd", "adam", "rmsprop", "adadelta"))) 
    stop("Invalid optimizer! Available options: sgd, adam, rmsprop, adadelta.\n")
  
  aseqs <- atseqs[[1]]
  tseqs <- atseqs[[2]]
  if(length(tseqs) != length(aseqs))
    stop("Number of timestamp sequences does not match the number of action sequences!\n")
  n_seq <- length(tseqs)
  
  n_actions_aseqs <- sapply(aseqs, length)
  n_actions_tseqs <- sapply(tseqs, length)
  if(!all(n_actions_aseqs==n_actions_tseqs))
    stop("Lengths of action sequences and lengths of timestamp sequences do not match!\n")
  
  events <- unique(unlist(aseqs))
  n_event <- length(events)
  
  # convert action sequence to one-hot vectors
  int_seqs <- list()
  target_seqs <- list()
  
  for (index_seq in 1:n_seq) {
    my_seq <- aseqs[[index_seq]]
    n_l <- length(my_seq)
    onehot_mat <- matrix(0, n_l, n_event)
    tmp <- match(my_seq, events)
    int_seqs[[index_seq]] <- matrix(tmp - 1, 1, n_l)
    onehot_mat[cbind(1:n_l, tmp)] <- 1
    target_seqs[[index_seq]] <- array(onehot_mat, dim=c(1,n_l, n_event))
  }
  
  # Get input/output time sequences
  if(cumulative){
    if(log){
      times_input <- sapply(tseqs, function(x) array(log10(x+1), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(log10(x+1), dim = c(1,length(tseq2interval(x)),1)))
    }else{
      times_input <- sapply(tseqs, function(x) array((x), dim = c(1,length(x))))
      times_output <- sapply(tseqs, function(x) array((x), dim = c(1,length(x),1)))
    }
  }else{
    if(log){
      times_input <- sapply(tseqs, function(x) array(log10(tseq2interval(x)+1), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(log10(tseq2interval(x)+1), dim = c(1,length(tseq2interval(x)),1)))
    }else{
      times_input <- sapply(tseqs, function(x) array(tseq2interval(x), dim = c(1,length(tseq2interval(x)))))
      times_output <- sapply(tseqs, function(x) array(tseq2interval(x), dim = c(1,length(tseq2interval(x)),1)))
    }
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # Define keras model
  # Define inputs and process it.
  encoder_inputs_aseq <- layer_input(shape=list(NULL))
  encoder_inputs_tseq <- layer_input(shape=list(NULL), dtype = 'float32')
  # expand RT sequence dimesion and concatenate with action embedding sequence
  expand_and_concatenate <- function(x){
    tmp <- k_expand_dims(x[[2]],axis = 3)
    return(k_concatenate(list(x[[1]],tmp),axis = 3))
  }
  
  if (rnn_type == "lstm")
  {
    if (method=="last") {
      encoder_aseq_embeddings <- encoder_inputs_aseq %>%
        layer_embedding(n_event, K)
      encoder_inputs <- list(encoder_aseq_embeddings,encoder_inputs_tseq)
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand_and_concatenate) %>%
        layer_lstm(units = K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[3]]
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[3]], axis=2)
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_aseq_embeddings <- encoder_inputs_aseq %>%
        layer_embedding(n_event, K)
      encoder_inputs <- list(encoder_aseq_embeddings,encoder_inputs_tseq)
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand_and_concatenate) %>%
        layer_lstm(units = K, return_sequences=TRUE)
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
    decoder_states <- decoder_inputs %>% 
      layer_lstm(units=K, return_sequences=TRUE)
    decoder_outputs_tseq <- decoder_states %>%
      layer_dense(1, activation='relu',name = 'tseq_output')
    decoder_outputs_aseq <- decoder_states %>%
      layer_dense(n_event, activation='softmax',name = 'aseq_output')
  } else if (rnn_type == "gru") {
    if (method=="last") {
      encoder_aseq_embeddings <- encoder_inputs_aseq %>%
        layer_embedding(n_event, K)
      encoder_inputs <- list(encoder_aseq_embeddings,encoder_inputs_tseq)
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand_and_concatenate) %>%
        layer_gru(units = K, return_sequences=TRUE, return_state=TRUE)
      encoder_outputs <- encoder_outputs_long[[2]]
      fn_rp <- function(x) {
        stepMatrix <- k_ones_like(x[[1]][,,1, drop=FALSE])
        latentMatrix <- k_expand_dims(x[[2]], axis=2)
        return(k_permute_dimensions(k_batch_dot(latentMatrix,stepMatrix, axes=list(2, 3)), list(1,3,2)))
      }
    } else {
      encoder_aseq_embeddings <- encoder_inputs_aseq %>%
        layer_embedding(n_event, K)
      encoder_inputs <- list(encoder_aseq_embeddings,encoder_inputs_tseq)
      encoder_outputs_long <- encoder_inputs %>%
        layer_lambda(expand_and_concatenate) %>%
        layer_gru(units = K, return_sequences=TRUE)
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
    decoder_states <- decoder_inputs %>% 
      layer_gru(units=K, return_sequences=TRUE)
    decoder_outputs_tseq <- decoder_states %>%
      layer_dense(1, activation='relu', name = 'tseq_output')
    decoder_outputs_aseq <- decoder_states %>%
      layer_dense(n_event, activation='softmax', name = 'aseq_output')
  }
  
  autoencoder_model <- keras_model(inputs = list(encoder_inputs_aseq,encoder_inputs_tseq),
                                   outputs = list(decoder_outputs_aseq,decoder_outputs_tseq))
  encoder_model <- keras_model(inputs = list(encoder_inputs_aseq,encoder_inputs_tseq), 
                               outputs = encoder_outputs)
  
  # Run training
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  autoencoder_model %>% compile(optimizer=optimizer, 
                                loss = list('aseq_output' = 'categorical_crossentropy','tseq_output'='mean_squared_error'),
                                loss_weights= list('aseq_output' = weights[1], 'tseq_output' = weights[2]))
  
  best_valid_loss <- Inf
  valid_loss <- rep(0, n_epoch)
  train_loss <- rep(0, n_epoch)
  if (!is.null(samples_test)) test_loss <- rep(0, n_epoch)
  
  theta <- matrix(0, n_seq, K)
  
  for (index_epoch in 1:n_epoch) {
    if (verbose) cat("Epoch", index_epoch, "\n")
    samples_train <- sample(samples_train)
    
    for (index_seq in samples_train) {
      autoencoder_model %>% train_on_batch(list(int_seqs[[index_seq]],times_input[[index_seq]]),
                                           list(target_seqs[[index_seq]],times_output[[index_seq]]))
    }
    for (index_seq in samples_train) {
      train_loss[index_epoch] <- train_loss[index_epoch] + evaluate(autoencoder_model,
                                                                    list(int_seqs[[index_seq]],times_input[[index_seq]]), 
                                                                    list(target_seqs[[index_seq]],times_output[[index_seq]]), verbose=0)[[1]]
    }
    for (index_seq in samples_valid) {
      valid_loss[index_epoch] <- valid_loss[index_epoch] + evaluate(autoencoder_model, 
                                                                    list(int_seqs[[index_seq]],times_input[[index_seq]]), 
                                                                    list(target_seqs[[index_seq]],times_output[[index_seq]]), verbose=0)[[1]]
    }
    for (index_seq in samples_test) {
      test_loss[index_epoch] <- test_loss[index_epoch] + evaluate(autoencoder_model, 
                                                                  list(int_seqs[[index_seq]],times_input[[index_seq]]), 
                                                                  list(target_seqs[[index_seq]],times_output[[index_seq]]), verbose=0)[[1]]
    }
    train_loss[index_epoch] <- train_loss[index_epoch] / length(samples_train)
    valid_loss[index_epoch] <- valid_loss[index_epoch] / length(samples_valid)
    if (!is.null(samples_test)) 
      test_loss[index_epoch] <- test_loss[index_epoch] / length(samples_test)
    
    if (valid_loss[index_epoch] < best_valid_loss) {
      best_valid_loss <- valid_loss[index_epoch]
      if (return_theta) {
        for(index_seq in 1:n_seq) theta[index_seq,] <- as.numeric(predict(encoder_model, list(int_seqs[[index_seq]],times_input[[index_seq]])))
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


