library(tracerer)

if(file.exists("beast.log")) {
	# Parse the log file
	beast_log_full <- parse_beast_tracelog_file("beast.log")

	# Remove the burn-in
	beast_log <- remove_burn_ins(
	  beast_log_full,
	  burn_in_fraction = 0.1
	)

	# Calculates the effective sample sizes of all parameter estimates
	esses <- calc_esses(beast_log, sample_interval = 5000)
	names(esses) <- NULL
	rownames(esses) <- NULL
	v <- unlist(esses)
	write(v, file="v.out")
	cat(any(v<300))
} else {
	cat(TRUE)
}
