#' @export
sub_seqs <- function(seqs, ids) {
  proc(seqs$action_seqs[ids], seqs$time_seqs[ids])
}

# remove
remove_action <- function(seqs, actions) {
  
  n_seq <- length(seqs$action_seqs)
  
  for (i in 1:n_seq) {
    aseq <- seqs$action_seqs[[i]]
    tseq <- seqs$time_seqs[[i]]
    rm_pos <- which(aseq %in% actions)
    if (length(rm_pos) > 0) {
      seqs$action_seqs[[i]] <- aseq[-rm_pos]
      seqs$time_seqs[[i]] <- tseq[-rm_pos]
    }
  }
  
  seqs
  
}

# replace
replace_action <- function(seqs, old_action, new_action) {
  n_seq <- length(seqs$action_seqs)
  
  action_seqs <- seqs$action_seqs
  
  action_seqs <- sapply(action_seqs, function(x) {x[x==old_action] <- new_action; x})
  
  seqs$action_seqs <- action_seqs
  
  seqs
}

# combine
combine_action <- function(seqs, old_actions, new_action) {
  
}