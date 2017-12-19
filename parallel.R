# select parallel backend

opt.warn <- getOption("warn")
options(warn = -1)

library(foreach)
library(parallel)
library(iterators)

options(warn = opt.warn)

parallel.start <<- function(cores) { invisible(0) }
parallel.end <<- function() { invisible(0) }
parallel.backend <<- NA

if (sum(installed.packages()[,"Package"] == "doMC")) {
	parallel.backend <<- "doMC"
	message("loading ", parallel.backend)
	library(doMC)
	parallel.start <<- function(cores) { 
		registerDoMC(cores)
	}
	parallel.end <<- function() { 
		invisible(0)
	}
} else if (sum(installed.packages()[,"Package"] == "doSNOW")) {
	parallel.backend <<- "doSNOW"
	message("loading ", parallel.backend)
	library(doSNOW)
	parallel.start <<- function(cores) { 
		paralle.cl <<- makeCluster(cores, type="SOCK")
		registerDoSNOW( paralle.cl )
	}
	parallel.end <<- function() { 
		stopCluster(paralle.cl)
	}
} else if (sum(installed.packages()[,"Package"] == "doParallel")) {
	parallel.backend <<- "doParallel"
	message("loading ", parallel.backend)
	library(doParallel)
	parallel.start <<- function(cores) { 
		paralle.cl <<- makeCluster(cores)
		registerDoParallel( paralle.cl )
	}
	parallel.end <<- function() { 
		stopCluster(paralle.cl)
	}
}

