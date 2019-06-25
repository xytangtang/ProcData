rm(list=ls())

cc_data_raw <- read.csv("~/Documents/Research/Process_Data/data/PISA_2012_released/climateproc.csv", stringsAsFactors=F)

cc_data_log <- cc_data_raw[cc_data_raw$Event!="diagram", c("ID", "time", "Event")]

cc_data_log$Event <- gsub(" ", "", cc_data_log$Event, fixed = TRUE)

cc_logfile <- split(cc_data_log, f=cc_data_log$ID)

cc_data_outcome <- cc_data_raw[,c("ID", "Correctness")]

outcome_split <- split(cc_data_outcome$Correctness, cc_data_outcome$ID)


cc_response <- sapply(outcome_split, function(x) x[1])
cc_action_seqs <- sapply(cc_logfile, function(x) x$Event)
cc_time_seqs <- sapply(cc_logfile, function(x) x$time)
cc_time_interval_seqs <- sapply(cc_time_seqs, diff)
rm_obs <- which(sapply(cc_time_interval_seqs, function(x) any(x < 0)))
cc_time_seqs <- cc_time_seqs[-rm_obs]
cc_action_seqs <- cc_action_seqs[-rm_obs]

cc_seqs <- proc(action_seqs = cc_action_seqs, time_seqs = cc_time_seqs)
cc_data <- list(seqs = cc_seqs, responses = cc_response)
  
# devtools::use_data(cc_data, overwrite = TRUE)
usethis::use_data(cc_data, overwrite = TRUE)
