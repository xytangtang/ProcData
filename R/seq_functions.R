#' Action sequence generator
#'
#' \code{seq_gen} generates action sequences of an imaginery simulation-based item.
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
#' @param n An integer. The number of sequences to be generated.
#' @param action_set1,action_set2 Character vectors giving the choices for 
#'   the first and the second conditions.
#' @param answer_set A character vector giving the choices for the answer. 
#' @param p1,p2 Nonnegative numeric vectors. They are the weights for sampling 
#'   from \code{action_set1} and \code{action_set2}.
#' @param p_answer A nonnegative numeric vector giving the weights for sampling
#'   from \code{answer_set}.
#' @param p_continue Probability of running an/another experiment.
#' @param p_choose Probability of choosing an answer.
#' @return A list of \code{n} elements. Each element is a generated action sequence.
#' @examples 
#' seqs <- seg_gen(100)
#' 
#' @family sequence generators
#' @export
seq_gen <- function(n, action_set1 = c("OPT1_1", "OPT1_2", "OPT1_3"), 
                    action_set2 = c("OPT2_1", "OPT2_2"), 
                    answer_set = c("CHECK_A", "CHECK_B", "CHECK_C", "CHECK_D"), 
                    p1 = rep(1,length(action_set1)), p2 = rep(1, length(action_set2)), 
                    p_answer = rep(1, length(answer_set)), p_continue = 0.5, p_choose = 0.5)
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
#' @param n An integer. The number of sequences to be generated.
#' @param events A character vector specifying the set of \code{N} possible
#'   actions. Default is \code{letters}.
#' @param Pmat An \code{N} by \code{N} probability transition matrix.
#' @param start_index Index of the action indicating the start of an item in
#'   \code{events}.
#' @param end_index Index of the action indicating the end of an item in
#'   \code{events}.
#' @param max_len Maximum length of generated sequences.
#' @return A list of \code{n} elements. Each element is a generated action
#'   sequence.
#' @examples 
#' seqs <- seq_gen2(100)
#' 
#' @family sequence generators
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

