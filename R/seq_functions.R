
#' Count action appearances
#' 
#' This function counts the appearances of each action in \code{actions} in 
#' action sequence \code{x}.
#' 
#' @param x an action sequence.
#' @param actions a set of actions whose number of appearances will be count.
#' @return an integer vector of counts.
#' @export
count_actions <- function(x, actions) 
{
	sapply(actions, function(a) sum(x==a))
}

# action transitions each row is an action transition (bigram).
aseq2atranseqs <- function(x) {
  l <- length(x)
  
  rbind(x[-l], x[-1])
}

#' Summarize action sequences
#' 
#' @param action_seqs a list of action sequences.
#' @return a list containing the following objects:
#'   \item{n_seq}{the number of action sequences}
#'   \item{n_action}{the number of distinct actions}
#'   \item{action}{the action set}
#'   \item{seq_length}{sequence lengths}
#'   \item{action_freq}{action counts}
#'   \item{action_seqfreq}{the number of sequences that each action appears}
#'   \item{trans_count}{a \code{length(action)} by \code{length(action)} matrix whose
#'     element in the i-th row and j-th column is the counts of transition from 
#'     \code{action[i]} to \code{action[j]}.}
#' @seealso \code{\link{time_seqs_summary}} for summarizing timestamp sequences.
#' @export
action_seqs_summary <- function(action_seqs)
{
  # number of sequences
	n_seq <- length(action_seqs)
	# sequence length
	seq_length <- sapply(action_seqs, length)
	# action set
	actions <- sort(unique(unlist(action_seqs)))
	# number of actions
	n_action <- length(actions)
	
	# action frequency
	action_freq_by_seq <- t(sapply(action_seqs, count_actions, actions=actions))
	action_freq <- colSums(action_freq_by_seq)
	# action sequence frequency
	action_seq_counts_by_seq <- array(as.numeric(action_freq_by_seq > 0), 
	                                  dim=dim(action_freq_by_seq))
	action_seq_freq <- colSums(action_seq_counts_by_seq)
	names(action_seq_freq) <- actions
	
	# tfidf weighted action frequency
	# inv_seq_freq <- log(n_seq / action_seq_freq)
	# term_freq <- log(action_freq_by_seq)
	# term_freq[term_freq==-Inf] <- 0
	# term_freq <- 1 + term_freq
	# tfidf <- term_freq * (rep(1, n_seq) %*% rbind(inv_seq_freq))
	# weighted_action_freq <- colSums(tfidf)
	
	# action transtion counts 
	trans_counts <- matrix(0, n_action, n_action)
	colnames(trans_counts) <- actions
	rownames(trans_counts) <- actions
	
	action_tran_seqs <- sapply(action_seqs, aseq2atranseqs)
	all_pairs <- matrix(unlist(action_tran_seqs), nrow=2)
	n_pair <- ncol(all_pairs)
	for (i in 1:n_pair) 
	  trans_counts[all_pairs[1,i], all_pairs[2,i]] <- trans_counts[all_pairs[1,i], all_pairs[2,i]] + 1 
	
	list(n_seq=n_seq, n_action = n_action,
	     actions=actions, seq_length=seq_length, 
	     action_freq=action_freq, action_seqfreq = action_seq_freq, 
	     # weighted_action_freq = weighted_action_freq, 
	     trans_count = trans_counts)
}

#' Transform a timestamp sequence into a inter-arrival time sequence
#' 
#' @param x a timestamp sequence
#' @return a numeric vector of the same length as \code{x}. The first element in 
#'   the returned vector is 0. The t-th returned element is \code{x[t] - x[t-1]}.
#' @export
tseq2interval <- function(x) {
  c(0, diff(x))
}

#' Summarize timestamp sequences
#' 
#' @param time_seqs a list of timestamp sequences
#' @return a list containing the following objects
#'   \item{total_time}{a summary of response time of \code{n_seq} response processes}
#'   \item{mean_react_time}{a summary of mean reaction time of \code{n_seq} response processes}
#' @export
time_seqs_summary <- function(time_seqs) {
  # total time
  total_time <- summary(sapply(time_seqs, max))
  # time per action
  mean_react_time <- summary(sapply(time_seqs, max) / sapply(time_seqs, length))
  # inter-arrival time seqs
  #time_interval <- summary(unlist(sapply(time_seqs, tseq2interval)))
  
  list(total_time=total_time, mean_react_time=mean_react_time)
}



