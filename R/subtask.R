
calculate_entropy <- function(p_vec) {
  -sum(p_vec * log(p_vec))
}


# step 1: fit action prediction model and compute entropy processes
#' Step 1 of Subtask Analysis: obtaining entropy sequences of action sequences
#' 
#' \code{action2entropy} fit a recurrent-neural-network-based action prediction 
#' model to a set of action sequences \code{action_seqs}.
#' 
#' @inheritParams seqm
#' @param action_seqs a list of action sequences
#' @param rnn_dim latent dimension of RNN
#' @return \code{action2entropy} returns a list containing
#' \item{entropy_seqs}{a list of entropy sequences. The length of each entropy 
#'   sequence is one less than that of the corresponding action sequence.}
#' \item{loss_history}{a \code{n_epoch} by 2 matrix. The two columns contain
#'   the loss at the end of each epoch on the training set and the validation 
#'   set, respectively.}
#' \item{rnn_dim}{the latent dimension of the recurrent neural network}
#' \item{model_fit}{a vector of class \code{"raw"}. It is the serialized version of 
#'     the trained keras model.}
#' \item{actions}{a vector of the actions in \code{action_seqs}.}
#' \item{max_len}{maximum length of the action sequences.}
#' @seealso \code{\link{entropy2segment}} and \code{\link{segment2subtask}} for 
#'   steps 2 and 3 of the subtask analysis procedure; \code{\link{subtask_analysis}} 
#'   for the complete procedure.
#' @references Wang, Z., Tang, X., Liu, J., and Ying, Z. (2020) Subtask analysis of process data
#'   through a predictive model. \link{https://arxiv.org/abs/2009.00717}
#' @export
action2entropy <- function(action_seqs, 
                           rnn_dim = 20, 
                           n_epoch = 50, 
                           step_size = 0.001, 
                           batch_size = 1,
                           optimizer_name = "rmsprop",
                           index_valid = 0.2,
                           verbose = FALSE) {
  
  n_person <- length(action_seqs)
  
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
  
  events <- unique(unlist(action_seqs))
  n_event <- length(events)
  
  max_len <- max(sapply(action_seqs, length))
  
  int_seqs <- matrix(0, n_person, max_len-1)
  target_seqs <- array(0, dim=c(n_person, max_len-1, n_event))
  
  for (index_seq in 1:n_person) {
    my_seq <- action_seqs[[index_seq]]
    n_l <- length(my_seq)
    tmp <- match(my_seq, events)
    int_seqs[index_seq, 1:(n_l-1)] <- tmp[1:(n_l-1)]
    for (index_t in 1:(n_l - 1)) target_seqs[index_seq, index_t, tmp[index_t + 1]] <- 1
  }
  
  Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  
  # build keras model
  seq_inputs <- layer_input(shape=list(max_len - 1))
  seq_emb <- seq_inputs %>% layer_embedding(n_event + 1, rnn_dim, mask_zero=TRUE)
  hidden_states <- seq_emb %>% layer_gru(units=rnn_dim, return_sequences = TRUE)
  out_probs <- hidden_states %>% layer_dense(units=n_event, activation = "softmax")

  action_pred_model <- keras_model(seq_inputs, out_probs)
  
  if (optimizer_name == "sgd") optimizer <- optimizer_sgd(lr=step_size)
  else if (optimizer_name == "rmsprop") optimizer <- optimizer_rmsprop(lr=step_size)
  else if (optimizer_name == "adadelta") optimizer <- optimizer_adadelta(lr=step_size)
  else if (optimizer_name == "adam") optimizer <- optimizer_adam(lr=step_size)
  
  action_pred_model %>% compile(optimizer = optimizer, loss='categorical_crossentropy')
  
  best_valid <- Inf
  trace_res <- matrix(0, n_epoch, 2)
  colnames(trace_res) <- c("train", "valid")
  if (verbose) cat("Training action prediction model...\n")
  for (index_epoch in 1:n_epoch) {
    if (verbose) cat("Epoch ", index_epoch, "\n")
      model_res <- action_pred_model %>% fit(int_seqs[index_train, ], 
                                             target_seqs[index_train,,], 
                                             epochs = 1, 
                                             batch_size = batch_size, 
                                             verbose = 0, 
                                             validation_data=list(int_seqs[index_valid, ], 
                                                                  target_seqs[index_valid,,]))
    trace_res[index_epoch, 1] <- model_res$metrics$loss
    trace_res[index_epoch, 2] <- model_res$metrics$val_loss
    if (model_res$metrics$val_loss < best_valid) {
      best_valid <- model_res$metrics$val_loss
      model_save <- serialize_model(action_pred_model)
    }
  }
  if (verbose) cat("Done\n")
  
  if (verbose) cat("Computing predictive distributions...  ")
  action_pred_model <- unserialize_model(model_save)
  pred_res <- predict(action_pred_model, x = int_seqs)
  k_clear_session()
  if (verbose) cat("Done\n")
  
  if (verbose) cat("Computing entropy processes...  ")
  entropy_seqs <- list()
  for (index_seq in 1:n_person) {
    n_l <- length(action_seqs[[index_seq]])
    pred_dist <- pred_res[index_seq, 1:(n_l - 1), ]
    entropy_seqs[[index_seq]] <- apply(pred_dist, 1, calculate_entropy)
  }
  if (verbose) cat("Done\n")
  
  
  res <- list(entropy_seqs = entropy_seqs,
              loss_history = trace_res,
              rnn_dim = rnn_dim,
              model_fit = model_save,
              actions = events, 
              max_len = max_len) 

  res
  
}

# take one entropy sequence and tunnig parameter lambda and return (internal) segment points
filter_function <- function(entropy_seq, lambda) {
  l <- length(entropy_seq)
  if (l < 3) return(numeric(0))
  threshold <- lambda * (max(entropy_seq) - min(entropy_seq))
  seq_ext <- entropy_seq
  seq_ext[0] <- Inf
  seq_ext[l] <- Inf
  
  # collect local maxima
  local_max <- numeric(0)
  for (i in 2:(l-1)) {
    if (seq_ext[i] >= max(seq_ext[i-1], seq_ext[i+1])) 
      local_max <- c(local_max, i)
  }
  
  # filter U curve and output locations for segmentation
  bds <- numeric(0)
  index_begin <- 1
  for (index_end in local_max) {
    v_low <- min(seq_ext[index_begin:index_end])
    v_high <- min(seq_ext[c(index_begin, index_end)])
    if (v_high - v_low > threshold) {
      bds <- c(bds, index_end)
      index_begin <- index_end
    }
  }
  
  bds
  
}

#' Segment an entropy sequence
#' 
#' \code{segment_function} segments the entropy sequence \code{entropy_seq} 
#' by identifying deep U-shaped curves in it.  
#' 
#' @param entropy_seq a vector of entropies
#' @param lambda a number between 0 and 1.
#' @return a vector of segment boundaries
#' @export
segment_function <- function(entropy_seq, lambda) {
  l <- length(entropy_seq)
  segs_l2r <- filter_function(entropy_seq, lambda)
  segs_r2l <- filter_function(rev(entropy_seq), lambda)
  segs_r2l <- l - segs_r2l + 1 # plus one or minus one?
  segs <- sort(unique(c(segs_l2r, segs_r2l, l)))
  
  segs
}

#' Step 2 of Subtask Analysis: Segmenting Entropy Sequences
#' 
#' \code{entropy2segment} segments the entropy sequences in \code{entropy_seqs} 
#' using \code{segment_function}.
#' 
#' @param entropy_seqs a list of entropy sequences
#' @param lambda a number between 0 and 1
#' @param verbose print progress if TRUE. default is FALSE
#' @return a list containg the segment boundaries of each entropy sequence.
#' @seealso \code{\link{action2entropy}} and \code{\link{segment2subtask}} for 
#'   steps 1 and 3 of the subtask analysis procedure; \code{\link{subtask_analysis}} 
#'   for the complete procedure.
#' @references Wang, Z., Tang, X., Liu, J., and Ying, Z. (2020) Subtask analysis of process data
#'   through a predictive model. \link{https://arxiv.org/abs/2009.00717}
#' @export
entropy2segment <- function(entropy_seqs, lambda = 0.3, verbose = FALSE) {
  if (verbose) cat("Segmenting entropy sequences...  ")
  res <- sapply(entropy_seqs, segment_function, lambda=lambda)
  if (verbose) cat("Done\n")
  
  list(seg_seqs = res)
  
}

get_seg_profile <- function(action_seq, segment_bds, actions) {
  profiles <- numeric()
  n_event <- length(actions)
  index_begin <- 1
  for (index_end in segment_bds) {
    freqs <- numeric(n_event)
    for (index_event in 1:n_event) {
      freqs[index_event] <- sum(action_seq[index_begin:index_end] == actions[index_event])
    }
    profiles <- rbind(profiles, freqs/sum(freqs))
    index_begin <- index_end
  }
  colnames(profiles) <- actions
  profiles
  
}
  
elbow_cluster_selection <- function(x, k) {
  n <- length(x)
  slopes <- -diff(x)/diff(k)
  slopes_pct_change <- numeric(n)
  for (i in 1:(n-2)) {
    slopes_pct_change[i] <- (slopes[i] - slopes[i+1])/slopes[i]
  }
  
  selected_index <- which.max(slopes_pct_change) + 1 # check add 1 or 2
  
  selected_index
}

derive_subtask_seqs <- function(action_seqs, seg_seqs, cluster_ids) {
  
  n_person <- length(action_seqs)
  
  subtask_seqs <- list()
  index_profile <- 1
  for (index_person in 1:n_person) {
    subtask_seq <- rep(0, length(action_seqs[[index_person]]))
    index_begin <- 1
    for (index_end in seg_seqs[[index_person]]) {
      if (length(seg_seqs[[index_person]]) == 1) subtask_seq[] <- cluster_ids[index_profile]
      else if (index_end == seg_seqs[[index_person]][1]) subtask_seq[(index_begin):(index_end)] <- cluster_ids[index_profile]
      else if (index_end < max(seg_seqs[[index_person]])) subtask_seq[(index_begin+1):(index_end)] <- cluster_ids[index_profile]
      else subtask_seq[(index_begin+1):(index_end+1)] <- cluster_ids[index_profile]
      index_profile <- index_profile + 1
      index_begin <- index_end
    }
    subtask_seqs[[index_person]] <- subtask_seq
  }
  
  subtask_seqs
}

compute_relative_cluster_profile <- function(action_seqs, subtask_seqs, actions) {
  all_elements <- unlist(action_seqs)
  all_subtasks <- unlist(subtask_seqs)
  
  n_subtask <- length(unique(all_subtasks))
  n_event <- length(actions)
  
  relative_cluster_profiles <- matrix(0, n_subtask, n_event)
  
  for (index_cluster in 1:n_subtask) {
    select_elements <- all_subtasks == index_cluster
    for (index_action in 1:n_event) {
      relative_cluster_profiles[index_cluster, index_action] <- sum(all_elements[select_elements] == actions[index_action])
    }
  }
  overall_profile <- colSums(relative_cluster_profiles)
  overall_profile <- overall_profile / sum(overall_profile)
  relative_cluster_profiles <- relative_cluster_profiles/rowSums(relative_cluster_profiles)
  for (index_cluster in 1:n_subtask) relative_cluster_profiles[index_cluster, ] <- relative_cluster_profiles[index_cluster, ]/overall_profile
  
  colnames(relative_cluster_profiles) <- actions
  relative_cluster_profiles
}

#' Step 3 of Subtask Analysis: Grouping Segments
#' 
#' \code{segment2subtask} clustering action sequence segments according to their action
#' frequency profiles. Each cluster forms a subtask.
#' 
#' @inheritParams action2entropy
#' @param seg_seqs a list of segment locations
#' @param n_subtask the desired number of subtasks or a vector of candidate number of subtasks
#' @param actions a set of actions
#' @param ... additional arguments passed to \code{kmeans}
#' @return a list containing 
#' \item{n_subtask}{the number of subtasks}
#' \item{subtasks}{a vector of subtasks}
#' \item{subtask_seqs}{a list of subtask sequences}
#' \item{tot.withinss}{a vector of total within cluster sum of squares}
#' \item{relative_cluster_profiles}{a \code{n_subtask} by \code{length(actions)} matrix. Each
#'   row gives the action frequency profile of each subtask relative to the overall action 
#'   frequency profile}
#' @seealso \code{\link{action2entropy}} and \code{\link{segment2subtask}} for 
#'   steps 1 and 3 of the subtask analysis procedure; \code{\link{subtask_analysis}} 
#'   for the complete procedure. 
#' @references Wang, Z., Tang, X., Liu, J., and Ying, Z. (2020) Subtask analysis of process data
#'   through a predictive model. \link{https://arxiv.org/abs/2009.00717}
#' @export
segment2subtask <- function(action_seqs, seg_seqs, n_subtask, actions, verbose = FALSE, ...) {
  
  n_person <- length(action_seqs)
  
  if (length(action_seqs) != length(seg_seqs)) stop("action sequences and segment sequences do not match!\n")
  if (any(sapply(action_seqs, length) - 1 < sapply(seg_seqs, max))) stop("segment location is out of bounds!\n")
  
  if (length(n_subtask) == 1) {
    if (n_subtask < 1) stop("n_subtask should be an integer greater than one!\n")
    else {
      select_num_cluster <- FALSE
      n_subtask <- round(n_subtask)
    }
  } else {
    if (any(n_subtask < 1)) stop("all elements of n_subtask should be greater than one!\n")
    else {
      select_num_cluster <- TRUE
      n_subtask_candidate <- sort(unique(n_subtask))
      n_subtask_candidate <- unique(round(n_subtask_candidate))
    }
  }
  

  if (verbose) cat("Computing action frequency profiles...  ")
  all_profiles <- numeric()
  for (index_person in 1:n_person) {
    profiles <- get_seg_profile(action_seqs[[index_person]], seg_seqs[[index_person]], actions)
    all_profiles <- rbind(all_profiles, profiles)
  }
  if (verbose) cat("Done\n")
  
  if (verbose) cat("Clustering segments...  ")
  if (select_num_cluster) {
    n_cand <- length(n_subtask_candidate)
    cluster_loss <- rep(0, n_cand)
    cluster_ids_all <- matrix(0, nrow(all_profiles), n_cand)
    cat("\n")
    for (index_cand in 1:n_cand) {
      if (verbose) cat("Candidate ", index_cand, "\n")
      cluster_res <- kmeans(all_profiles, n_subtask_candidate[index_cand], ...)
      cluster_loss[index_cand] <- cluster_res$tot.withinss
      cluster_ids_all[, index_cand] <- cluster_res$cluster
    }
    
    # select number of cluster by elbow method
    selected_index <- elbow_cluster_selection(cluster_loss, n_subtask_candidate)
    n_subtask <- n_subtask_candidate[selected_index]
    cluster_ids <- cluster_ids_all[, selected_index]
    
  } else {
    cluster_res <- kmeans(x=all_profiles, centers=n_subtask, ...)
    cluster_loss <- cluster_res$tot.withinss
    cluster_ids <- cluster_res$cluster
  }
  if (verbose) cat("Done\n")
  
  if (verbose) cat("Deriving subtask sequences...  ")
  subtask_seqs <- derive_subtask_seqs(action_seqs, seg_seqs, cluster_ids)
  if (verbose) cat("Done\n")
  
  if (verbose) cat("Computing relative action frequency profile for each cluster...  ")
  relative_cluster_profiles <- compute_relative_cluster_profile(action_seqs, subtask_seqs, actions)
  if (verbose) cat("Done\n")
  
  res <- list(n_subtask = n_subtask, subtasks=1:n_subtask, subtask_seqs=subtask_seqs, tot.withinss=cluster_loss, relative_cluster_profiles = relative_cluster_profiles)
  
  res
  
}

#' Subtask Analysis
#' 
#' \code{subtask_analysis} performs subtask identification procedure.
#' 
#' @inheritParams action2entropy
#' @inheritParams entropy2segment
#' @inheritParams segment2subtask
#' @return an object of class \code{"subtask"}. It is a list containing
#' \item{action_seqs}{a list of action sequences}
#' \item{entropy_seqs}{a list of entropy sequences}
#' \item{seg_seqs}{a list of segment boundaries}
#' \item{subtask_seqs}{a list of subtask sequences}
#' \item{subtasks}{a vector of subtasks}
#' \item{n_subtask}{the number of subtasks}
#' \item{tot.withinss}{a vector of total within cluster sum of squares}
#' \item{relative_cluster_profiles}{a \code{n_subtask} by \code{length(actions)} matrix. Each
#'   row gives the action frequency profile of each subtask relative to the overall action 
#'   frequency profile}
#' \item{loss_history}{a \code{n_epoch} by 2 matrix. The two columns contain
#'   the loss at the end of each epoch on the training set and the validation 
#'   set, respectively.}
#' \item{rnn_dim}{the latent dimension of the recurrent neural network}
#' \item{model_fit}{a vector of class \code{"raw"}. It is the serialized version of 
#'     the trained action prediction model.}
#' \item{actions}{a vector of the actions in \code{action_seqs}.}
#' \item{max_len}{maximum length of the action sequences.}
#' @seealso \code{\link{action2entropy}}, \code{\link{entropy2segment}}, and 
#'   \code{\link{segment2subtask}} for the three steps of subtask analysis.
#' @references Wang, Z., Tang, X., Liu, J., and Ying, Z. (2020) Subtask analysis of process data
#'   through a predictive model. \link{https://arxiv.org/abs/2009.00717}
#' @export
subtask_analysis <- function(action_seqs, lambda=0.3, n_subtask, 
                             rnn_dim = 20,
                             n_epoch = 20, 
                             step_size = 0.001, 
                             batch_size = 1,
                             optimizer_name = "rmsprop",
                             index_valid = 0.2,
                             verbose = FALSE, ...) {
  
  step1_res <- action2entropy(action_seqs = action_seqs,
                              rnn_dim = rnn_dim,
                              n_epoch = n_epoch,
                              step_size = step_size,
                              batch_size = batch_size,
                              optimizer_name = optimizer_name,
                              index_valid = index_valid,
                              verbose = verbose)
 
  step2_res <- entropy2segment(entropy_seqs = step1_res$entropy_seqs, lambda = lambda, verbose = verbose)
  
  step3_res <- segment2subtask(action_seqs = action_seqs, seg_seqs = step2_res$seg_seqs, n_subtask = n_subtask, actions = step1_res$actions, verbose = verbose, ...)
  
  res <- c(list(action_seqs = action_seqs), step1_res, step2_res, step3_res)
  class(res) <- "subtask"
  
  res
}

subtask2seg <- function(subtask_seq) {
  x <- rep(0, length(subtask_seq))
  for (i in 1:(length(x)-1)) {
    if (subtask_seq[i] != subtask_seq[i+1]) x[i] <- 1
  }
  which(x == 1)
}


#' Plot Subtask Analysis Results for One Sequence
#' 
#' @param action_seq an action sequence
#' @param entropy_seq an entropy sequence
#' @param subtask_seq a subtask sequence
#' @param subtasks a vector of all subtasks
#' @param col.subtask a vector of colors for subtasks 
#' @param lty line types 
#' @param pch point characters
#' @param plot_legend a logical value. If \code{TRUE} (default), plot the legend of subtasks.
#' @param legend_pos a character string or the coordinates to be used to position the legend.
#' @param ... other arguments passed to \code{\link{legend}}
#' @return this function does not return values
#' @seealso \code{\link{plot_subtask_seqs}} for plotting results for all sequences. 
#'   \code{\link{plot.subtask}} for the plot method of \code{"subtask"} object.
#' @export
plot_subtask_seq <- function(action_seq, entropy_seq, subtask_seq , subtasks, col.subtask = 1:length(subtasks), cex.action = 0.5, lty = 1, pch = 16, srt = -90, plot_legend = TRUE, legend_pos = "topleft", ...) {
  l_action <- length(action_seq)
  l_entropy <- length(entropy_seq)
  l_subtask <- length(subtask_seq)
  n_subtask <- length(subtasks)
  names(col.subtask) <- subtasks
  
  yrange <- max(entropy_seq) - min(entropy_seq)
  ylim <- range(entropy_seq) + c(-0.1*yrange, 0.1*yrange)
  plot(entropy_seq, type = "p", ylim=ylim, col = col.subtask[subtask_seq[-l_subtask]], xlab="Time Steps", ylab="Entropy", pch=pch)
  segments(1:(l_entropy-1), entropy_seq[1:(l_entropy-1)], 2:l_entropy, entropy_seq[2:l_entropy], col = col.subtask[subtask_seq[2:(l_subtask-1)]], lty=lty)
  text(entropy_seq, action_seq[-l_action], col = col.subtask[subtask_seq[-l_subtask]], cex = cex.action, srt = srt, pos=1, offset = 1)
  if (plot_legend) legend(legend_pos, legend=subtasks, col=col.subtask, lty = lty, pch=16, ...)
  
  invisible()
}

simplify_subtask_seq <- function(subtask_seq, max_len) {
  x <- subtask_seq
  l <- length(x)
  x1 <- x[-1]
  x2 <- x[-l]
  rm_id <- which(x1 == x2)
  if (length(rm_id) > 0) x <- x[-rm_id]
  
  if (length(x) > max_len) x <- x[1:max_len]
  
  x
}

#' Plot Subtask Analysis Results for Entire Dataset
#' 
#' @inheritParams plot_subtask_seq
#' @param subtask_seqs a list of subtask sequences
#' @param max_len maximum length of plotted subtasks
#' @return this function does not return values
#' @seealso \code{\link{plot_subtask_seq}} for ploting results for one sequence. 
#'   \code{\link{plot.subtask}} for ploting an object of class \code{"subtask"}
#' @export
plot_subtask_seqs <- function(subtask_seqs, subtasks, max_len=5, col.subtask = 1:length(subtasks), plot_legend = TRUE, legend_pos="topright", ...) {
  op <- par(mar=c(1,1,1,1), oma=c(2,2,2,0))
  
  plot_seqs <- sapply(subtask_seqs, simplify_subtask_seq, max_len = max_len)
  max_subtask <- max(sapply(plot_seqs, length))
  max_len <- min(max_subtask, max_len)
  n_seq <- length(plot_seqs)
  o <- order(sapply(plot_seqs, paste, collapse=""))
  plot(0, xlim=c(0, max_len), ylim=c(0, n_seq), type="n", xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=F)
  mtext("Subtasks", side=1, line=0.5, outer=F, at=max_len/2)
  mtext("Respondents", side=2, line=0.5, outer=F, at=n_seq/2)
  for (i in 1:n_seq) {
    seqi <- plot_seqs[[o[i]]]
    li <- length(seqi)
    for (l in 1:li) {
      rect(xleft=l-1, xright=l, ybottom = n_seq - i, ytop=n_seq - i + 1, col=col.subtask[seqi[l]], border=NA)
    }
  }
  if (plot_legend) legend(legend_pos, legend=subtasks, xpd=TRUE, fill=col.subtask, border=col.subtask, ...)

  on.exit(par(op), add=TRUE)
  invisible()
}

#' Plot an subtask Object
#' 
#' Plot the subtask analysis results for either the entire dataset or individual sequences.
#' 
#' @inheritParams plot_subtask_seq
#' @inheritParams plot_subtask_seqs
#' @param object an object of class \code{"subtask"}
#' @param type \code{"all"} or \code{"individual"}
#' @param index a vector of indices of sequences to plot
#' @return this function does not return values
#' @seealso \code{\link{plot_subtask_seq}}, \code{\link{plot_subtask_seqs}}.
#' @export
plot.subtask <- function(object, type="all", index = NULL, max_len = 5, col.subtask = 1:length(object$subtasks), cex.action = 0.5, lty = 1, pch = 16, srt = -90, plot_legend=TRUE, legend_pos = "topright", ...) {

  if (type=="all") plot_subtask_seqs(object$subtask_seqs, object$subtasks, max_len = max_len, col.subtask = col.subtask, plot_legend = plot_legend, legend_pos = legend_pos, ...)
  else if (type=="individual") {
    if (is.null(index)) stop("index is empty!\n")
    n_plot_seq <- length(index) 
    for (index_seq in index) {
      plot_subtask_seq(object$action_seqs[[index_seq]], 
                       object$entropy_seqs[[index_seq]], 
                       object$subtask_seqs[[index_seq]], 
                       object$subtasks, 
                       col.subtask = col.subtask, 
                       cex.action = cex.action,
                       lty = lty, 
                       pch = pch, 
                       srt = srt,
                       plot_legend = plot_legend, 
                       legend_pos = legend_pos, ...)
    }
  }
  invisible()
}

