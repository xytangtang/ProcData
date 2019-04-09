#' Sequence generator
#'
#' \code{seq_gen} generates action sequences similar to PISA data: choose the levels of two conditions to run simulations and select an answer
#'
#' @param n number of sequences to be generated
#' @param action_set1 choices for the first condition
#' @param action_set2 choices for the second condition
#' @param answer_set choices for the answer 
#' @param p1 probability weights to sample from \code{action_set1}
#' @param p2 probability weights to sample from \code{action_set2}
#' @param p_answer probability weights to sample from \code{answer_set}
#' @param p_continue probability of running another simulation
#' @param p_choose probability of choose an answer
#' @return a list of \code{n} elements. Each element is a sequence of actions.
#' @export
seq_gen <- function(n, action_set1 = c("OPT1_1", "OPT1_2", "OPT1_3"), action_set2 = c("OPT2_1", "OPT2_2"), answer_set = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"), p1 = rep(1,length(action_set1)), p2 = rep(1, length(action_set2)), p_answer = rep(1, length(answer_set)), p_continue = 0.5, p_choose = 0.5)
{
  seqs <- list()
  for (i in 1:n)
  {
    cur_seq <- c("Start")
    while(runif(1) < p_continue)
    {
      cur_seq <- c(cur_seq, sample(action_set1, size = 1, prob = p1), sample(action_set2, size = 1, prob = p2), "RUN")
    }
    if (runif(1) < p_choose) cur_seq <- c(cur_seq, sample(answer_set, size = 1, prob=p_answer))
    cur_seq <- c(cur_seq, "End")
    seqs[[i]] <- cur_seq
  }
  
  seqs
}

#' Action sequence generator
#'
#' \code{seq_gen2} generates action sequences according to a given probability transition matrix
#'
#' @param n number of sequences to be generated
#' @param events a set of \code{N} possible actions including start and end
#' @param Pmat \code{N} by \code{N} probability transition matrix; if \code{Pmat} not supplied, each action is randomly drawn from \code{events} (excluding start).
#' @param start_index index of start in \code{events}
#' @param end_index index of end in \code{events}
#' @param max_len maximum length of generated sequences
#' @return a list of \code{n} elements. Each element is a sequence of actions.
#' @export
seq_gen2 <- function(n, Pmat = NULL, events = letters, start_index=1, end_index=length(events), max_len=200)
{
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
    while(event_index != end_index)
    {
      event_index <- sample(1:n_event, 1, prob=Pmat[event_index,])
      int_seq <- c(int_seq, event_index)
    }
    seqs[[i]] <- events[int_seq]
  }
  seqs
}

#' Action counts in a sequence
#' 
#' @param x an action sequences
#' @param actions a set of actions whose number of appearances will be count

count_actions <- function(x, actions) 
{
	sapply(actions, function(a) sum(x==a))
}

#' Find all actions appears in a list of action sequences
#' 
#' @param seqs a list of action sequences
#' @param count logical value. If TRUE, unigram frequence will also be output
find_unigram <- function(seqs, count)
{
	n_seq <- length(seqs)
	unigrams <- unique(seqs[[1]])
	
	for (i in 2:n_seq) unigrams <- union(unigrams, unique(seqs[[i]]))
	
	if (!count) return(unigrams)
	
	unigram_counts <- count_actions(seqs[[1]], unigrams)
	
	for (i in 2:n_seq) unigram_counts <- unigram_counts + count_actions(seqs[[i]], unigrams)
	
	list(unigram=unigrams, count=unigram_counts)
}

seq2bseq <- function(myseq, sep="\t") # given a unigram sequence, output its bigram version
{
	seq_len <- length(myseq)
	
	paste(myseq[1:(seq_len-1)], myseq[2:seq_len], sep=sep)
}

count_bigram <- function(seqs) # given a list of sequences, find all appeared bigrams and count frequencies
{	
	n_seq <- length(seqs)
	
	bigrams <- unique(seq2bseq(seqs[[1]]))
	
	for (i in 2:n_seq) bigrams <- union(bigrams, unique(seq2bseq(seqs[[i]])))
		
	bigram_counts <- count_actions(seq2bseq(seqs[[1]]), bigrams)
	for (i in 2:n_seq) bigram_counts <- bigram_counts + count_actions(seq2bseq(seqs[[i]]), bigrams)
	
	list(bigram=bigrams, count=bigram_counts)
}

seq2position <- function(myseq, actions) # given a sequence and a set of actions (unigrams or bigrams or combination), find the position these actions appear in the sequence
{
	action_counts <- count_actions(myseq, actions)
	
	action_positions <- list()
	
	n_action <- length(actions)
	
	for (i in 1:n_action) action_positions[[i]] <- which(myseq == actions[i])
	
	list(counts = action_counts, positions = action_positions)
}

