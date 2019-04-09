
calculate_f1 <- function(seq1, seq2)
{
	common_events <- intersect(seq1, seq2)
	n_common_event <- length(common_events)
	f_score <- 0
	for (event in common_events)
	{
		set1 <- which(seq1 == event)
		set2 <- which(seq2 == event)

		n_common <- min(length(set1), length(set2))

		f_score <- f_score + sum(abs(set1[1:n_common] - set2[1:n_common]))
	}

	f_score <- f_score / max(length(seq1), length(seq2))

	f_score
}

calculate_time_f1 <- function(seq1, seq2, ts1, ts2) # t1, t2 are time stamps, not time intervals
{
	common_events <- intersect(seq1, seq2)
	n_common_event <- length(common_events)
	f_score <- 0
	for (event in common_events)
	{
		set1 <- which(seq1 == event)
		set2 <- which(seq2 == event)

		n_common <- min(length(set1), length(set2))

		f_score <- f_score + sum(abs(ts1[set1[1:n_common]] - ts2[set2[1:n_common]]))
	}

	n1 <- length(seq1)
	n2 <- length(seq2)
	f_score <- f_score / max(ts1[n1], ts2[n2])

	f_score
}

calculate_f2 <- function(seq1, seq2)
{
	common_events <- intersect(seq1, seq2)
	n_common_event <- length(common_events)

	f_score <- 0

	if (n_common_event == 0) return(f_score)

	for (event in common_events)
	{
		set1 <- which(seq1 == event)
		set2 <- which(seq2 == event)

		n_common <- min(length(set1), length(set2))

		f_score <- f_score + sum(abs(set1[1:n_common] - set2[1:n_common])) + abs(length(set1) - length(set2))
	}

	f_score <- f_score / n_common_event

	f_score
}

calculate_g <- function(seq1, seq2)
{
	sum(!(seq1 %in% seq2)) + sum(!(seq2 %in% seq1))
}

calculate_time_g <- function(seq1, seq2, t1, t2)
{
	n1 <- length(seq1)
	n2 <- length(seq2)
	t1_pre <- c(0, t1[-n1])
	t2_pre <- c(0, t2[-n2])
	ti1 <- t1 - t1_pre
	ti2 <- t2 - t2_pre

	sum(ti1[!(seq1 %in% seq2)]) + sum(ti2[!(seq2 %in% seq1)])
}

calculate_tdissimilarity <- function(seq1, seq2, t1, t2)
{
	n1 <- length(seq1)
	n2 <- length(seq2)
	f_score <- calculate_time_f1(seq1, seq2, t1, t2)
	g_score <- calculate_time_g(seq1, seq2, t1, t2)
	res <- (f_score + g_score) / (t1[n1] + t2[n2])

	return(res)
}

calculate_dissimilarity <- function(seq1, seq2, method)
{
# method: f1
#		  f2
#		  g
#		  oss
#		  osa
#		  edit
#		  lendiff

	if (method == "f1") res <- calculate_f1(seq1, seq2)
	else if (method == "f2") res <- calculate_f2(seq1, seq2)
	else if (method == "g") res <- calculate_g(seq1, seq2)
	else if (method == "oss") {
		n1 <- length(seq1)
		n2 <- length(seq2)
		f_score <- calculate_f1(seq1, seq2)
		g_score <- calculate_g(seq1, seq2)
		res <- (f_score + g_score) / (n1 + n2)
	}
	else if (method == "osa") {
		n1 <- length(seq1)
		n2 <- length(seq2)
		f_score <- calculate_f2(seq1, seq2)
		g_score <- calculate_g(seq1, seq2)
		res <- (f_score + g_score) / (n1 + n2)
	}
	else if (method == "edit") stop("Edit distance should be calculated in python!\n")
	else if (method == "lendiff") res <- abs(length(seq1) - length(seq2))
	else stop("Invalid distance method!\n")

	return(res)
}


