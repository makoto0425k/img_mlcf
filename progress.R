
progress.start <- function(total) {
	progress.time <<- Sys.time()
	progress.total <<- total
	message(format(progress.time, "%H:%M:%S"))	
}

progress.end <- function(skip = 0) {
	now <- Sys.time()
	elapsed <- floor(as.numeric(difftime(now, progress.time, units="secs")))
	elapsed.s <- floor(elapsed %% 60)
	elapsed.m <- floor(floor(elapsed / 60) %% 60)
	elapsed.h <- floor(elapsed / 3600)
	progress.report(progress.total, comment = sprintf("TOTAL %2.2d:%2.2d:%2.2d", elapsed.h, elapsed.m, elapsed.s), skip = skip)
}

progress.report <- function(current, skip = 0, comment = "", step = 1) {
	now <- Sys.time()
	elapsed <- floor(as.numeric(difftime(now, progress.time, units="secs")))
	time.report <- ""
	
	if (current != progress.total && step > 1 && current %% step) {
		return()
	}
	
	if (current - skip > 0) {
		estimate <- floor((elapsed) / (current - skip) * (progress.total - skip))
		rest <- floor(estimate / 60) - floor(elapsed / 60)
		rest.h <- floor(rest / 60)
		rest.m <- rest %% 60
		time.report <- sprintf(" %d-%d=%d %dh%dm", floor(estimate / 60), floor(elapsed / 60), rest, rest.h, rest.m)
	}
	report <- sprintf(" %d/%d(+%d) %d%%", 
										as.integer(floor(current - skip)), as.integer(floor(progress.total - skip)), as.integer(floor(skip)), 
										as.integer(floor((current - skip) / (progress.total - skip) * 100)))
	message(format(now, "%H:%M:%S"), report, time.report, " ", comment)	
}


