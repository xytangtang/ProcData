#' Subset response processes
#' 
#' @param seqs an object of class \code{"\link{proc}"}
#' @param ids a vector of indices 
#' @return an object of class \code{"\link{proc}"}
#' @examples 
#' data(cc_data)
#' seqs <- sub_seqs(cc_data$seqs, 1:10)
#' @export
sub_seqs <- function(seqs, ids) {
  proc(seqs$action_seqs[ids], seqs$time_seqs[ids])
}

#' Remove actions from response processes
#' 
#' Remove actions in \code{actions} and the corresponding timestamps
#'  in response processes \code{seqs}.
#'
#' @param seqs an object of class \code{"\link{proc}"}
#' @param actions a character vector. Each element is an action to be removed.
#' @return an object of class \code{"\link{proc}"} with actions in \code{actions}
#'   and the corresponding timestamps removed.
#' @examples 
#'   seqs <- seq_gen(10)
#'   new_seqs <- remove_action(seqs, c("RUN", "Start"))
#' @export
remove_action <- function(seqs, actions) {
  
  n_seq <- length(seqs$action_seqs)
  
  for (i in 1:n_seq) {
    aseq <- seqs$action_seqs[[i]]
    rm_pos <- which(aseq %in% actions)
    if (length(rm_pos) > 0) {
      seqs$action_seqs[[i]] <- aseq[-rm_pos]
      if (!is.null(seqs$time_seqs)) {
        tseq <- seqs$time_seqs[[i]]
        seqs$time_seqs[[i]] <- tseq[-rm_pos]
      }
    }
  }
  
  seqs
  
}

#' Replace actions in response processes
#' 
#' Replace \code{old_action} with \code{new_action} in \code{seqs}. Timestamp
#' sequences are not affected.
#' 
#' @param seqs an object of class \code{"\link{proc}"}
#' @param old_action a string giving the action to be replaced.
#' @param new_action a string giving the action replacing \code{old_action}
#' @return an object of class \code{"\link{proc}"}
#' @examples 
#' seqs <- seq_gen(10)
#' new_seqs <- replace_action(seqs, "Start", "Begin")
#' @export
replace_action <- function(seqs, old_action, new_action) {
  n_seq <- length(seqs$action_seqs)
  
  action_seqs <- seqs$action_seqs
  
  action_seqs <- sapply(action_seqs, function(x) {x[x==old_action] <- new_action; x})
  
  seqs$action_seqs <- action_seqs
  
  seqs
}

# check if pattern exist in action sequence x
check_pattern <- function(x, pattern)
{
  n <- length(x)
  l <- length(pattern)
  if (n < l) return(FALSE)
  for (i in 1:(n-l+1)) if (all(x[i:(i+l-1)] == pattern)) return(TRUE)
  
  return(FALSE)
}

#' Remove repeated actions
#' 
#' @param seqs an object of class \code{"\link{proc}"}
#' @param ignore repeated actions in ignore will not be deleted.
#' @return an object of class \code{"\link{proc}"}
#' @export
remove_repeat <- function(seqs, ignore=NULL) {
   aseqs <- seqs$action_seqs
   tseqs <- seqs$time_seqs
   
   n <- length(aseqs)
   for (i in 1:n) {
     aseq <- aseqs[[i]]
     tseq <- tseqs[[i]]
     l <- length(aseq)
     rm_idx <- which(aseq[1:(l-1)] == aseq[2:l] & !(aseq[2:l] %in% ignore))
     if (length(rm_idx) > 0) {
       rm_idx <- rm_idx + 1
       aseqs[[i]] <- aseq[-rm_idx]
       tseqs[[i]] <- tseq[-rm_idx]
     }
   }
   
   proc(action_seqs=aseqs, time_seqs=tseqs)
}

#' Combine consecutive actions into a single action
#' 
#' Combine the action pattern described in \code{old_actions} into a single action 
#' \code{new_action}. The timestamp of the combined action can be the timestamp of the
#' first action in the action pattern, the timestamp of the last action in the action
#' pattern, or the average of the two timestamps.
#' 
#' @param seqs an object of class \code{"\link{proc}"}
#' @param old_actions a character vector giving consecutive actions to be replaced.
#' @param new_action a string giving the combined action
#' @param timestamp "first", "last", or "avg", specifying how the timestamp of the combined
#'   action should be derived.
#' @return an object of class \code{"\link{proc}"}
#' @examples 
#' seqs <- seq_gen(100)
#' new_seqs <- combine_actions(seqs, 
#'                             old_actions=c("OPT1_3", "OPT2_2", "RUN"), 
#'                             new_action="KEY_ACTION")
#' @export
combine_actions <- function(seqs, old_actions, new_action, timestamp="first") {
  
  l <- length(old_actions)
  n_seq <- length(seqs$action_seqs)
  if (is.null(seqs$time_seqs)) {
    for (index_seq in 1:n_seq) {
      aseq <- seqs$action_seqs[[index_seq]]
      while(check_pattern(x=aseq, pattern=old_actions)) {
        n <- length(aseq)
        for (i in 1:(n-l+1)) {
          if (all(aseq[i:(i+l-1)] == old_actions)) break
        }
        aseq <- c(aseq[1:i], aseq[(i+l):n])
        aseq[i] <- new_action
      }
      seqs$action_seqs[[index_seq]] <- aseq
    }
  } else {
    for (index_seq in 1:n_seq) {
      aseq <- seqs$action_seqs[[index_seq]]
      tseq <- seqs$time_seqs[[index_seq]]
      while(check_pattern(x=aseq, pattern=old_actions)) {
        n <- length(aseq)
        for (i in 1:(n-l+1)) {
          if (all(aseq[i:(i+l-1)] == old_actions)) break
        }
        aseq <- c(aseq[1:i], aseq[(i+l):n])
        aseq[i] <- new_action
        
        if (timestamp=="last") tseq[i] <- tseq[i+l-1]
        else if (timestamp == "avg") tseq[i] <- (tseq[i+l-1] + tseq[i]) / 2
        tseq <- c(tseq[1:i], tseq[(i+l):n])
      }
      seqs$action_seqs[[index_seq]] <- aseq
      seqs$time_seqs[[index_seq]] <- tseq
    }
  }
  
  seqs
}