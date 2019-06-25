
#' count appearances of actions in \code{actions} in an action sequence
#' 
#' @param x an action sequence
#' @param actions a set of actions whose number of appearances will be count
count_actions <- function(x, actions) 
{
	sapply(actions, function(a) sum(x==a))
}

# action transitions
aseq2atranseqs <- function(x) {
  l <- length(x)
  
  rbind(x[-l], x[-1])
}

#' Summarize action sequences
#' 
#' @param action_seqs a list of action sequences
#' @return a list containing the following objects:
#'   \item{n}{number of action sequences}
#'   \item{action}{action set}
#'   \item{seq_length}{sequence lengths}
#'   \item{count}{action counts}
#'   \item{count_by_seq}{action counts in each sequence}
#'   \item{seq_count}{counts of number of sequences with a given action}
#'   \item{seq_count_by_seq}{if actions appears in a sequence}
#'   \item{trans_count}{action transition counts}
#' @export
action_seqs_summary <- function(action_seqs)
{
	n_seq <- length(action_seqs)
	
	seq_length <- sapply(action_seqs, length)
	
	actions <- sort(unique(unlist(action_seqs)))
	
	action_counts_by_seq <- t(sapply(action_seqs, count_actions, actions=actions))
	action_counts <- colSums(action_counts_by_seq)
	
	action_seq_counts_by_seq <- array(as.numeric(action_counts_by_seq > 0), 
	                                  dim=dim(action_counts_by_seq))
	action_seq_counts <- colSums(action_seq_counts_by_seq)
	
	n_action <- length(actions)
	trans_counts <- matrix(0, n_action, n_action)
	colnames(trans_counts) <- actions
	rownames(trans_counts) <- actions
	
	action_tran_seqs <- sapply(action_seqs, aseq2atranseqs)
	all_pairs <- matrix(unlist(action_tran_seqs), nrow=2)
	n_pair <- ncol(all_pairs)
	for (i in 1:n_pair) 
	  trans_counts[all_pairs[1,i], all_pairs[2,i]] <- trans_counts[all_pairs[1,i], all_pairs[2,i]] + 1 
	
	list(n=n_seq, action=actions, seq_length=seq_length, 
	     count=action_counts, count_by_seq=seq_freqs, 
	     seq_count=action_seq_counts, seq_count_by_seq=action_seq_counts_by_seq,
	     trans_count = trans_counts)
}

#' Transform a timestamp sequence to a time interval sequence
#' @param x a timestamp sequence
#' @export
tseq2interval <- function(x) {
  c(0, diff(x))
}

#' Summarize timestamp sequences
#' 
#' @param time_seqs a list of timestamp sequences
#' @return a list containing the following objects
#'   \item{total_time}{total time elapsed}
#'   \item{time_per_action}{average time between two consecutive actions}
#'   \item{time_interval_seqs}{time interval sequences}
#' @export
time_seqs_summary <- function(time_seqs) {
  # total time
  total_time <- sapply(time_seqs, max)
  # time per action
  time_per_action <- total_time / sapply(time_seqs, length)
  # inter-arrival time seqs
  time_interval_seqs <- sapply(time_seqs, tseq2interval)
  
  list(total_time=total_time, time_per_action=time_per_action, 
       time_interval_seqs=time_interval_seqs)
}



