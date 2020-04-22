#' Action sequence generator
#'
#' \code{seq_gen} generates action sequences of an imaginary simulation-based item.
#' 
#' The format of the generated sequences resembles that of the response processes of 
#' simulation-based items. In these items, participants are asked to answer
#' a question by running simulated experiments in which two conditions can 
#' be controlled. A simulated experiment can be run by setting the two conditions
#' at one of the given choices and click "Run" button.
#' 
#' The possible actions are "Start", "End", "Run", and the elements in \code{action_set1}, 
#' \code{action_set2}, and \code{answer_set}. The generated sequences begin with "Start"
#' and continue with groups of three actions. Each group of three actions, representing 
#' one experiment, consists of an action chosen from \code{action_set1} according to 
#' \code{p1}, an action chosen from \code{action_set2} according to \code{p2}, and "Run".
#' The probability of performing an experiment after "Start" or one experiment is 
#' \code{p_continue}. After the experiment process, with probability \code{p_choose}, an 
#' answer will be chosen. The chosen answer is randomly sampled from \code{answer_set} 
#' according to \code{p_answer}. All generated sequences end with "End".
#' 
#' @param n An integer. The number of action sequences to be generated.
#' @param action_set1,action_set2 Character vectors giving the choices for 
#'   the first and the second conditions.
#' @param answer_set A character vector giving the choices for the answer. 
#' @param p1,p2 Nonnegative numeric vectors. They are the weights for sampling 
#'   from \code{action_set1} and \code{action_set2}.
#' @param p_answer A nonnegative numeric vector giving the weights for sampling
#'   from \code{answer_set}.
#' @param p_continue Probability of running an/another experiment.
#' @param p_choose Probability of choosing an answer.
#' @param include_time logical. Indicate if timestamp sequences should be generated. Default is
#'   FALSE.
#' @param time_intv_dist A list specifying the distribution of the inter-arrival time.
#' @param seed Random seed.
#' @return An object of class \code{"\link{proc}"} with \code{time_seqs = NULL}.
#' @family sequence generators
#' @export
seq_gen <- function(n, action_set1 = c("OPT1_1", "OPT1_2", "OPT1_3"), 
                    action_set2 = c("OPT2_1", "OPT2_2"), 
                    answer_set = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"), 
                    p1 = rep(1,length(action_set1)), p2 = rep(1, length(action_set2)), 
                    p_answer = rep(1, length(answer_set)), p_continue = 0.5, p_choose = 0.5, 
                    include_time = FALSE, time_intv_dist = list("exp", 1), seed = 12345) {
  set.seed(seed)
  seqs <- list()
  for (i in 1:n) {
    cur_seq <- c("Start")
    while(runif(1) < p_continue) {
      cur_seq <- c(cur_seq, sample(action_set1, size = 1, prob = p1), sample(action_set2, size = 1, prob = p2), "RUN")
    }
    if (runif(1) < p_choose) cur_seq <- c(cur_seq, sample(answer_set, size = 1, prob=p_answer))
    cur_seq <- c(cur_seq, "End")
    seqs[[i]] <- cur_seq
  }
  
  if (include_time) {
    time_seqs <- list()
    dist_name <- time_intv_dist[[1]]
    
    for (i in 1:n) {
      l <- length(seqs[[i]])
      if (dist_name == "exp") {
        theta <- time_intv_dist[[2]]
        tseq <- rexp(l - 1, rate = theta)
        tseq <- c(0, cumsum(tseq))
      } else if (dist_name == "lognorm") {
        theta1 <- time_intv_dist[[2]]
        theta2 <- time_intv_dist[[3]]
        tseq <- rnorm(l - 1, mean = theta1, sd = theta2)
        tseq <- c(0, cumsum(exp(tseq)))
      }
      time_seqs[[i]] <- tseq
    }
  } else {
    time_seqs <- NULL
  }
  
  proc(action_seqs=seqs, time_seqs=time_seqs, ids=1:n)
}

#' Markov action sequence generator
#'
#' \code{seq_gen2} generates action sequences according to a given probability
#' transition matrix.
#'
#' This function generates \code{n} action sequences according \code{Pmat}. The
#' set of possible actions is \code{events}. All generated sequences start with
#' \code{events[start_index]} and end with \code{events[end_index]}. If
#' \code{Pmat} is not supplied, actions is uniformly drawn from
#' \code{events[-start_index]} until \code{events[end_index]} appears.
#'
#' @param n An integer. The number of action sequences to be generated.
#' @param events A character vector specifying the set of \code{N} possible
#'   actions. Default is \code{letters}.
#' @param Pmat An \code{N} by \code{N} probability transition matrix.
#' @param start_index Index of the action indicating the start of an item in
#'   \code{events}.
#' @param end_index Index of the action indicating the end of an item in
#'   \code{events}.
#' @param max_len Maximum length of generated sequences.
#' @param include_time logical. Indicate if timestamp sequences should be generated. Default is
#'   FALSE.
#' @param time_intv_dist A list specifying the distribution of the inter-arrival time.
#' @param seed random generator seed.
#' @return An object of class \code{"\link{proc}"} with \code{time_seqs = NULL}.
#' @family sequence generators
#' @export
seq_gen2 <- function(n, Pmat = NULL, events = letters, 
                     start_index=1, end_index=length(events), 
                     max_len=200, include_time = FALSE, 
                     time_intv_dist = list("exp", 1), seed=12345) {
  set.seed(seed)
  n_event <- length(events)
  if (is.null(Pmat)) {
    Pmat <- matrix(0, n_event, n_event)
    Pmat[-end_index, -start_index] <- 1 / (n_event - 1)
    Pmat[end_index, end_index] <- 1
  }
  
  seqs <- list()
  
  for (i in 1:n)
  {
    int_seq <- start_index
    event_index <- start_index
    while (event_index != end_index & length(int_seq) < max_len) {
      event_index <- sample(1:n_event, 1, prob = Pmat[event_index, ])
      int_seq <- c(int_seq, event_index)
    }
    if (tail(int_seq, 1) != end_index) int_seq <- c(int_seq, end_index)
    seqs[[i]] <- events[int_seq]
  }
  
  if (include_time) {
    time_seqs <- list()
    dist_name <- time_intv_dist[[1]]
    
    for (i in 1:n) {
      l <- length(seqs[[i]])
      if (dist_name == "exp") {
        theta <- time_intv_dist[[2]]
        tseq <- rexp(l - 1, rate = theta)
        tseq <- c(0, cumsum(tseq))
      } else if (dist_name == "lognorm") {
        theta1 <- time_intv_dist[[2]]
        theta2 <- time_intv_dist[[3]]
        tseq <- rnorm(l - 1, mean = theta1, sd = theta2)
        tseq <- c(0, cumsum(exp(tseq)))
      }
      time_seqs[[i]] <- tseq
    }
  } else {
    time_seqs <- NULL
  }
  
  proc(action_seqs=seqs, time_seqs=time_seqs, ids=1:n)
}

# Check if two lists have the same shape
same_shape <- function(target, current) {
  if (length(current) != length(target)) return(FALSE)
  n <- length(target)
  for (i in 1:n) {
    target_i <- target[[i]]
    current_i <- current[[i]]
    if (is.array(target_i)) {
      if (any(dim(current_i) != dim(target_i))) return(FALSE)
    } else {
      if (length(current_i) != length(target_i)) return(FALSE)
    }
  }
  
  return(TRUE)
}

#' RNN action sequence generator
#' 
#' \code{seq_gen3} generates action sequences according to a recurrent neural network
#' 
#' @inheritParams seq_gen2
#' @inheritParams seq2feature_seq2seq
#' @param rnn_type the type of recurrent unit to be used for generating sequences. 
#'   \code{"lstm"} for the long-short term memory unit. \code{"gru"} for the gated
#'   recurrent unit.
#' @param K the latent dimension of the recurrent unit.
#' @param weights a list containing the weights in the embedding layer, the recurrent 
#'   unit, the fully connected layer. If not (properly) specified, randomly generated 
#'   weights are used.
#' @param initial_state a list containing the initial state of the recurrent neural 
#'   network. If \code{rnn_type="lstm"}, it contains two 1 by \code{K} matrices. If
#'   \code{rnn_type="gru"}, it contains one 1 by \code{K} matrix. If not specified, 
#'   all the elements are set to zero.
#' @return A list containing the following elements
#'     \item{seqs}{an object of class \code{"\link{proc}"} with \code{time_seqs=NULL}.}
#'     \item{weights}{a list containing the weights used for generating sequences.}
#' @family sequence generators
#' @export
seq_gen3 <- function(n, events = letters, rnn_type = "lstm", K = 10, weights=NULL, 
                     max_len = 100, initial_state = NULL, start_index=1, 
                     end_index=length(events), include_time = FALSE,
                     time_intv_dist = list("exp", 1), gpu=FALSE, parallel=FALSE, 
                     seed=12345) {
  use_session_with_seed(seed, disable_gpu = !gpu, disable_parallel_cpu = !parallel)
  n_event <- length(events)
  
  if (!(rnn_type) %in% c("lstm", "gru")) 
    stop("Undefined type of RNN! Available options: lstm and gru.\n")
  
  if (rnn_type == "lstm") {
    state_c_inputs <- layer_input(shape=list(K))
    state_h_inputs <- layer_input(shape=list(K))
    state_inputs <- c(state_c_inputs, state_h_inputs)
  } else if (rnn_type == "gru") {
    state_inputs <- layer_input(shape=list(K))
  }
  seq_inputs <- layer_input(shape=list(NULL))
  emb <- seq_inputs %>% layer_embedding(n_event, K)
  if (rnn_type == "lstm") rnn_unit <- layer_lstm(units = K, return_sequences = FALSE, return_state=TRUE)
  else if (rnn_type == "gru") rnn_unit <- layer_gru(units = K, return_sequences = FALSE, return_state=TRUE)
  rnn_outputs <- rnn_unit(emb, initial_state=state_inputs)
  prob_outputs <- rnn_outputs[[1]] %>% layer_dense(n_event, activation="softmax", bias_initializer=initializer_random_uniform(minval=-5, maxval=5))
  if (rnn_type == "lstm") state_outputs <- rnn_outputs[2:3]
  else if (rnn_type == "gru") state_outputs <- rnn_outputs[2]
  seq_pred_model <- keras_model(c(seq_inputs, state_inputs), c(prob_outputs, state_outputs))
  
  if (!is.null(weights)) {
    curr_weights <- get_weights(seq_pred_model)
    if (same_shape(curr_weights, weights)) set_weights(seq_pred_model, weights)
    else warning("Provided weights does not match the model! Use randomly generated weights.\n")
  }
  
  seqs <- list()
  
  for (i in 1:n) {
    
    seq_index <- start_index - 1
    int_seqs <- c(seq_index + 1)
    if (is.null(initial_state)) {
      if (rnn_type == "lstm") initial_state <- list(array(0, dim=c(1, K)), array(0, dim=c(1,K)))
      else if (rnn_type == "gru") initial_state <- list(array(0, dim=c(1,K)))
    }
    while (seq_index != end_index - 1 && length(int_seqs) < max_len) {
      pred_res <- predict(seq_pred_model, c(list(matrix(seq_index, 1,1)), initial_state))
      prob_vec <- pred_res[[1]]
      if (rnn_type == "lstm") initial_state <- pred_res[2:3]
      else if (rnn_type == "gru") initial_state <- pred_res[2]
      seq_index <- sample(c(0:(n_event-1))[-start_index], 1, prob=prob_vec[-start_index])
      int_seqs <- c(int_seqs, seq_index + 1)
    }
    
    if (seq_index != end_index - 1) int_seqs <- c(int_seqs, end_index)
    
    seqs[[i]] <- events[int_seqs]
  } 
  weights = get_weights(seq_pred_model)
  k_clear_session()
  
  if (include_time) {
    time_seqs <- list()
    dist_name <- time_intv_dist[[1]]
    
    for (i in 1:n) {
      l <- length(seqs[[i]])
      if (dist_name == "exp") {
        theta <- time_intv_dist[[2]]
        tseq <- rexp(l - 1, rate = theta)
        tseq <- c(0, cumsum(tseq))
      } else if (dist_name == "lognorm") {
        theta1 <- time_intv_dist[[2]]
        theta2 <- time_intv_dist[[3]]
        tseq <- rnorm(l - 1, mean = theta1, sd = theta2)
        tseq <- c(0, cumsum(exp(tseq)))
      }
      time_seqs[[i]] <- tseq
    }
  } else {
    time_seqs <- NULL
  }
  
  seqs_res <- proc(action_seqs=seqs, time_seqs=time_seqs, ids=1:n)
  list(seqs=seqs_res, weights = weights)
  
}