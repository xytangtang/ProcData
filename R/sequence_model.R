K2string <- function(K_emb, K_rnn, K_hidden = NULL, include_time, rnn_type) {
  if (include_time) input_string <- "Action+Time"
  else input_string <- "Action"
  emb_string <- paste("Embed(K=", K_emb, ")", sep="")
  if (rnn_type == "lstm") rnn_string <- paste("LSTM(K=", K_rnn, ")", sep="") 
  else if (rnn_type == "gru") rnn_string <- paste("GRU(K=", K_rnn, ")", sep="") 
  
  dense_strings <- character(0)
  if (!is.null(K_hidden)) dense_strings <- c(dense_strings, paste("Dense(K=", K_hidden, ")", sep=""))
  dense_strings <- c(dense_strings, "Dense(K=1)")
  model_string <- paste(input_string, emb_string, rnn_string, paste(dense_strings, collapse=","), sep=",")
  
  model_string 
}

#' Fitting sequence models
#'
#' \code{seqm} is used to fit a neural network model relating a response process
#' with a variable.
#'
#' The model consists of an embedding layer, a recurrent layer and one or more
#' fully connected layers. The embedding layer takes an action sequence and
#' output a sequences of \code{K} dimensional numeric vectors to the recurrent
#' layer. If \code{include_time = TRUE}, the embedding sequence is combined with
#' the timestamp sequence in the response process as the input the recurrent
#' layer. The last output of the recurrent layer and the covariates specified in 
#' \code{covariates} are used as the input of the subsequent fully connected layer.
#' If \code{response_type="binary"}, the last layer uses the sigmoid activation
#' to produce the probability of the response being one. If
#' \code{response_type="scale"}, the last layer uses the linear activation. The
#' dimension of the output of other fully connected layers (if any) is specified
#' by \code{K_hidden}.
#'
#' The action sequences are re-coded into integer sequences and are padded with
#' zeros to length \code{max_len} before feeding into the model. If the provided
#' \code{max_len} is smaller than the length of the longest sequence in
#' \code{seqs}, it will be overridden.
#'
#' @inheritParams seq2feature_seq2seq
#' @param response response variable.
#' @param covariates covariate matrix.
#' @param response_type "binary" or "scale".
#' @param actions a character vector gives all possible actions. It is will be
#'   expanded to include all actions appear in \code{seqs} if necessary.
#' @param include_time logical. If the timestamp sequence should be included in the model.
#' @param time_interval logical. If the timestamp sequence is included as a sequence of 
#'   inter-arrival time.
#' @param log_time logical. If take the logarithm of the time sequence.
#' @param n_hidden the number of hidden fully-connected layers.
#' @param K_emb the latent dimension of the embedding layer.
#' @param K_rnn the latent dimension of the recurrent neural network.
#' @param K_hidden a vector of length \code{n_hidden} specifying the number of
#'   nodes in each hidden layer.
#' @param n_epoch the number of training epochs.
#' @param batch_size the batch size used in training.
#' @param index_valid proportion of sequences used as the validation set or a vector 
#'   of indices specifying the validation set.
#' @param max_len the maximum length of response processes.
#' @return \code{seqm} returns an object of class \code{"seqm"}, which is a list containing
#'   \item{structure}{a string describing the neural network structure.}
#'   \item{coefficients}{a list of fitted coefficients. The length of the list is 
#'     6 + 2 * \code{n_hidden}. The first element gives the action embedding. 
#'     Elements 2-4 are parameters in the recurrent unit. The rest of the elements are 
#'     for the fully connected layers. Elements 4 + (2 * i - 1) and 4 + 2 * i give the parameters
#'     for the i-th fully connected layer.}
#'   \item{model_fit}{a vector of class \code{"raw"}. It is the serialized version of 
#'     the trained keras model.} 
#'   \item{feature_model}{a vector of class \code{"raw"}. It is the serialized version of the
#'     keras model for obtaining the rnn outputs.}
#'   \item{include_time}{if the timestamp sequence is included in the model.}
#'   \item{time_interval}{if inter-arrival time is used.}
#'   \item{log_time}{if the logarithm time is used.}
#'   \item{actions}{all possible actions.}
#'   \item{max_len}{the maximum length of action sequences.}
#'   \item{history}{a \code{n_epoch} by 2 matrix giving the training and
#'   validation losses at the end of each epoch.} 
#'
#' @seealso \code{\link{predict.seqm}} for the \code{predict} method for \code{seqm} objects.
#' @examples
#' \donttest{
#' n <- 100
#' data(cc_data)
#' samples <- sample(1:length(cc_data$responses), n)
#' seqs <- sub_seqs(cc_data$seqs, samples)
#' 
#' y <- cc_data$responses[samples]
#' x <- matrix(rnorm(n*2), ncol=2)
#' 
#' index_test <- 91:100
#' index_train <- 1:90
#' seqs_train <- sub_seqs(seqs, index_train)
#' seqs_test <- sub_seqs(seqs, index_test)
#' 
#' actions <- unique(unlist(seqs$action_seqs))
#' 
#' ## no covariate is used
#' res1 <- seqm(seqs = seqs_train, response = y[index_train], 
#'              response_type = "binary", actions=actions, K_emb = 5, K_rnn = 5, 
#'              n_epoch = 5)
#' pred_res1 <- predict(res1, new_seqs = seqs_test)
#' 
#' mean(as.numeric(pred_res1 > 0.5) == y[index_test])
#' 
#' ## add more fully connected layers after the recurrent layer.
#' res2 <- seqm(seqs = seqs_train, response = y[index_train],
#'              response_type = "binary", actions=actions, K_emb = 5, K_rnn = 5, 
#'              n_hidden=2, K_hidden=c(10,5), n_epoch = 5)
#' pred_res2 <- predict(res2, new_seqs = seqs_test)
#' mean(as.numeric(pred_res2 > 0.5) == y[index_test])
#' 
#' ## add covariates
#' res3 <- seqm(seqs = seqs_train, response = y[index_train], 
#'              covariates = x[index_train, ],
#'              response_type = "binary", actions=actions, 
#'              K_emb = 5, K_rnn = 5, n_epoch = 5)
#' pred_res3 <- predict(res3, new_seqs = seqs_test, 
#'                      new_covariates=x[index_test, ])
#'                      
#' ## include time sequences
#' res4 <- seqm(seqs = seqs_train, response = y[index_train], 
#'              response_type = "binary", actions=actions,
#'              include_time=TRUE, K_emb=5, K_rnn=5, n_epoch=5)
#' pred_res4 <- predict(res4, new_seqs = seqs_test)
#' }
#' @export
seqm <- function(seqs, response, covariates = NULL, response_type, 
                 actions = unique(unlist(seqs$action_seqs)), 
                 rnn_type = "lstm", include_time = FALSE, time_interval = TRUE, 
                 log_time = TRUE, K_emb = 20, K_rnn = 20, n_hidden = 0, 
                 K_hidden = NULL, index_valid = 0.2, 
                 verbose = FALSE, 
                 max_len = NULL, n_epoch = 20, batch_size = 16, optimizer_name = "rmsprop", 
                 step_size = 0.001, gpu = FALSE, parallel = FALSE, seed = 12345) {
  use_session_with_seed(seed, disable_gpu = !gpu, disable_parallel_cpu = !parallel)
  n_person <- length(seqs$action_seqs)
  if (is.null(actions)) events <- unique(unlist(seqs$action_seqs))
  else {
    events <- union(unique(unlist(seqs$action_seqs)), actions)
    if (length(events) > length(actions)) 
      warning("Action set is expanded to include all actions in action sequences.\n")
  }
  n_event <- length(events)
  
  if (!include_time) {
    seqs$time_seqs <- NULL
    time_interval <- NA
    log_time <- NA
  } else if (is.null(seqs$time_seqs)) {
    stop("Timestamp sequences are not available!\n")
  } else {
    if (time_interval)
      seqs$time_seqs <- sapply(seqs$time_seqs, tseq2interval)
    if (log_time)
      seqs$time_seqs <- sapply(seqs$time_seqs, function(x) log10(1.0 + x))
  }
  
  max_len0 <- max(sapply(seqs$action_seqs, length))
  if (is.null(max_len)){
    max_len <- max_len0
  } else if (max_len < max_len0) {
    warning("max_len is set as the max length in seqs!\n")
    max_len <- max_len0
  }
  
  if (is.null(index_valid)) {
    stop("Validation set is empty! Set proportion or indices of validation set 
         through index_valid.\n")
  } else if (length(index_valid) == 1) {
    if (index_valid < 1 && index_valid > 0) {
      index_valid <- sample(1:n_person, round(n_person*index_valid))
      index_train <- setdiff(1:n_person, index_valid)
    } else {
      stop("Invalid validation set! Set index_valid as a number between zero and 1 or as 
           a vector of indices.\n")
    }
  } else {
    if (!all(index_valid == round(index_valid)) || !all(index_valid > 0)) {
      stop("Invalid validation set! Set index_valid as a number between zero and 1 or as 
           a vector of indices.\n")
    } else {
      index_valid <- index_valid[index_valid <= n_person]
      index_train <- setdiff(1:n_person, index_valid)
    }
  }
  
  int_seqs <- matrix(0, n_person, max_len)
  
  for (index_seq in 1:n_person) {
    my_seq <- seqs$action_seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:n_l] <- tmp
  }
  
  # pad time sequence with -1
  if (include_time) {
    time_seqs <- array(-1.0, dim=c(n_person, max_len, 1))
    for (index_seq in 1:n_person) {
      tseq <- seqs$time_seqs[[index_seq]]
      n_l <- length(tseq)
      time_seqs[index_seq, 1:n_l, 1] <- tseq
    }
  }
  
  if (!gpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, K_emb, mask_zero=TRUE)
  if (include_time) {
    time_inputs <- layer_input(shape=list(max_len, 1))
    time_masked <- time_inputs %>% layer_masking(mask_value = -1.0)
    seq_emb <- layer_concatenate(list(seq_emb, time_masked), axis=-1)
  }
  
  if (is.null(covariates)) {
    K_cov <- 0
    covariates <- matrix(numeric(0), nrow = length(response), ncol=0)
  }
  else K_cov <- ncol(covariates)
  cov_inputs <- layer_input(shape=list(K_cov))
  
  if (rnn_type == "lstm") seq_feature <- seq_emb %>% layer_lstm(units=K_rnn)
  else if (rnn_type == "gru") seq_feature <- seq_emb %>% layer_gru(units=K_rnn)
  
  outputs <- layer_concatenate(inputs=c(seq_feature, cov_inputs), axis=-1)
  n_hidden <- min(n_hidden, length(K_hidden))
  
  if (n_hidden > 0) {
    for (index_hidden in 1:n_hidden) 
      outputs <- outputs %>% layer_dense(units=K_hidden[index_hidden], activation='tanh')
  }
  
  if (response_type == "binary") 
    outputs <- outputs %>% layer_dense(units=1, activation='sigmoid')
  else if (response_type == "scale") 
    outputs <- outputs %>% layer_dense(units=1, activation='linear')
  
  if (include_time) {
    seq_model <- keras_model(c(seq_inputs, time_inputs, cov_inputs), outputs)
    feature_model <- keras_model(c(seq_inputs, time_inputs), seq_feature)
  }
  else {
    seq_model <- keras_model(c(seq_inputs, cov_inputs), outputs)
    feature_model <- keras_model(seq_inputs, seq_feature)
  }
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  if (response_type == "binary") 
    seq_model %>% compile(optimizer = optimizer, loss='binary_crossentropy')
  else if (response_type == "scale") 
    seq_model %>% compile(optimizer = optimizer, loss='mean_squared_error')
  
  best_valid <- Inf
  trace_res <- matrix(0, n_epoch, 2)
  colnames(trace_res) <- c("train", "valid")
  for (index_epoch in 1:n_epoch) {
    if (include_time)
      model_res <- seq_model %>% fit(list(int_seqs[index_train, , drop = FALSE], 
                                          time_seqs[index_train, , , drop = FALSE], 
                                          covariates[index_train, , drop = FALSE]), 
                                     response[index_train], 
                                     epochs=1, batch_size=batch_size, verbose = verbose, 
                                     validation_data=list(list(int_seqs[index_valid,], 
                                                               time_seqs[index_valid,,,drop=FALSE],
                                                               covariates[index_valid,]), 
                                                          response[index_valid]))
    else  
      model_res <- seq_model %>% fit(list(int_seqs[index_train, ], 
                                          covariates[index_train,]), 
                                     response[index_train], 
                                     epochs = 1, batch_size = batch_size, verbose = verbose, 
                                     validation_data=list(list(int_seqs[index_valid, , drop=FALSE], 
                                                               covariates[index_valid, , drop=FALSE]), 
                                                          response[index_valid]))
    trace_res[index_epoch, 1] <- model_res$metrics$loss
    trace_res[index_epoch, 2] <- model_res$metrics$val_loss
    if (model_res$metrics$val_loss < best_valid) {
      best_valid <- model_res$metrics$val_loss
      model_save <- serialize_model(seq_model)
      feature_model_save <- serialize_model(feature_model)
    }
  }
  weights <- get_weights(seq_model)
  k_clear_session()
  model_string <- K2string(K_emb, K_rnn, K_hidden[1:n_hidden], rnn_type, include_time = include_time)  
  res <- list(structure = model_string, coefficients = weights, n_cov = K_cov,
              model_fit = model_save, feature_model = feature_model_save, 
              include_time = include_time, time_interval = time_interval,
              log_time = log_time, actions = events, max_len = max_len, history = trace_res) 
  class(res) <- "seqm"
  
  res  
}

#' Predict method for sequence models
#' 
#' Obtains predictions from a fitted sequence model object.
#' 
#' It unserialize \code{object$model_fit} to obtain a keras model of class 
#' \code{"keras.engin.training.Model"} and then calls \code{predict}
#' to obtain predictions.
#' 
#' @param object a fitted object of class \code{"seqm"} from \code{seqm}.
#' @param new_seqs an object of class \code{"\link{proc}"} with which to predict.
#' @param new_covariates a new covariate matrix with which to predict.
#' @param type a string specifying whether to predict responses (\code{"response"}) 
#'   or features (\code{"feature"}) or both (\code{"both"}).
#' @param ... further arguments to be passed to \code{predict.keras.engine.training.Model}.
#' 
#' @return If \code{type="response"}, a vector of predictions. The vector gives the 
#'   probabilities of the response variable being one if \code{response_type="binary"}.
#'   If \code{type="feature"}, a matrix of rnn outputs. If \code{type="both"}, a list 
#'   containing both the vector of response variable prediction and the rnn output matrix.
#' @seealso \code{\link{seqm}} for fitting sequence models.
#' @export
predict.seqm <- function(object, new_seqs, new_covariates = NULL, type="response", ...) {
  
  max_len <- object$max_len
  events <- object$actions
  n_cov <- object$n_cov
  n_seq <- length(new_seqs$action_seqs)
  
  if (is.null(new_covariates)) {
    new_covariates <- matrix(numeric(0), nrow = n_seq, ncol=0)
    if (n_cov > 0) 
      stop("No covariate is provided!\n")
  } else {
    if (n_cov == 0) {
      warning("Covariates are not used!\n")
      new_covariates <- matrix(numeric(0), nrow = n_seq, ncol=0)
    }
    else if (n_cov != ncol(new_covariates)) {
      stop("New covariates does not match covariates in the model!\n")
    }
  }
  if (type=="response" || type=="both") seq_model <- unserialize_model(object$model_fit)
  if (type=="feature" || type=="both") feature_model <- unserialize_model(object$feature_model)
  
  n <- length(new_seqs$action_seqs)
  new_int_seqs <- matrix(0, n, max_len)
  for (index_seq in 1:n) {
    my_seq <- new_seqs$action_seqs[[index_seq]]
    n_l <- min(length(my_seq), max_len)
    tmp <- match(my_seq, events)
    new_int_seqs[index_seq, 1:n_l] <- tmp[1:n_l]
  }
  
  if (object$include_time && !is.null(new_seqs$time_seqs)) {
    new_time_seqs <- array(-1, dim=c(n, max_len, 1))
    for (index_seq in 1:n) {
      tseq <- new_seqs$time_seqs[[index_seq]]
      n_l <- min(length(tseq), max_len)
      if (object$time_interval) tseq <- tseq2interval(tseq)
      if (object$log_time) tseq <- log(1.0 + tseq)
      new_time_seqs[index_seq, 1:n_l, 1] <- tseq[1:n_l]
    }
  } else if (object$include_time && is.null(new_seqs$time_seqs)) {
    stop("Timestamp sequences are not available!\n")
  }
  
  if (object$include_time) {
    if (type=="response" || type=="both") 
      response_res <- predict(seq_model, list(new_int_seqs, new_time_seqs, new_covariates), ...)
    if (type=="feature" || type=="both") 
      feature_res <- predict(feature_model, list(new_int_seqs, new_time_seqs), ...)
  }
  else {
    if (type=="response" || type=="both") 
      response_res <- predict(seq_model, list(new_int_seqs, new_covariates), ...)
    if (type=="feature" || type=="both")
      feature_res <- predict(feature_model, new_int_seqs, ...)
  }
  k_clear_session()
  
  if (type=="response") res <- response_res
  if (type=="feature") res <- feature_res
  if (type=="both") res <- list(response=response_res, feature=feature_res)
  
  res
}
