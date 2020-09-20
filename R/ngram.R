
collect_ngrams <- function(x, n, sep="\t") {
  if (n == 1) return(x)
  
  L <- length(x)
  if (L < n) return(NULL)
  
  ngrams <- x[1:(L - (n - 1))]
  for (i in 1:(n-1)) {
    ngrams <- paste(ngrams, x[i + 1:(L - (n - 1))], sep = sep)
  }
  ngrams
}

log_plus_one <- function(x) {
  ifelse(x == 0, 0, 1 + log(x))
}

#' ngram feature extraction
#' 
#' @family feature extraction methods
#' 
#' @param seqs an object of class \code{"\link{proc}"}
#' @param level an integer specifying the max length of ngrams
#' @param type a character string (\code{"binary"}, \code{"freq"}, or \code{"weighted"}) specifying the type of ngram features.
#' @return a matrix of ngram features
#' @examples 
#' seqs <- seq_gen(100)
#' theta <- seq2feature_ngram(seqs)
#' @export
seq2feature_ngram <- function(seqs, level = 2, type = "binary", sep="\t") {
  
  if (!(type %in% c("binary", "freq", "weighted"))) stop("Undefined ngram feature type!\n")
  if (class(seqs) != "proc") stop("seqs should be a proc object!\n")
  level <- round(level)
  
  theta <- numeric(0)
  action_seqs <- seqs$action_seqs
  n_seq <- length(action_seqs)
  
  for (index_level in 1:level) {

    lgram_seqs <- sapply(action_seqs, collect_ngrams, n=index_level, sep=sep)
    lgram_vec <- unlist(lgram_seqs)
    lgrams <- unique(lgram_vec)
    
    n_lgram <- length(lgrams)
    
    lgram_tf <- matrix(0, n_seq, n_lgram)
    colnames(lgram_tf) <- lgrams
    
    for (index_seq in 1:n_seq) {
      lgrams_seq <- lgram_seqs[[index_seq]]
      for (lgram in lgrams_seq) {
        lgram_tf[index_seq, lgram] <- lgram_tf[index_seq, lgram] + 1
      }
    }
    
    if (type == "freq") theta <- cbind(theta, lgram_tf)
    else if (type == "binary") {
      lgram_tf_binary <- (lgram_tf > 0) + 0
      theta <- cbind(theta, lgram_tf_binary)
    }
    else {
      lgram_sf <- numeric(n_lgram)
      names(lgram_sf) <- lgrams
      
      for (lgram in lgrams) {
        lgram_sf[lgram] <- sum(sapply(lgram_seqs, function(x) lgram %in% x))
      }
      
      lgram_weight <- t(log(n_seq / lgram_sf) * t(log_plus_one(lgram_tf)))
      
      theta <- cbind(theta, lgram_weight)
    }
    
  }
  
  return(theta)
}

