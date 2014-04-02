#######
#######
#######
#######
# Functions for creating/analyzing simulated data
#######
#######
#######
#######



######
# Function to switch subregions
# reps is either: vector of when to switch 
	# or number of monitors in each subregion
# i is monitor
# start is current counter
switchfun <- function(reps, i, start) {
	
	#if unequal
	if(length(reps) > 1) {
		#switch if changepoint
		startchange <- i %in% reps
		
	#otherwise switch every reps monitors
	}else {
		startchange <- i %% reps == 0
	}
	
	#update counter
	if(startchange) {
		start <- start + 1
	}
	start
}





#####
# Function to create data
# nr is number of observations
# vec is PCs from sample data
# sources is which sources in vec to keep
# cm is lognormal means
# sd is lognormal sds
createdat <- function(nr, vec, sources, cm, sd) {
	#keep which sources
	vec2 <- as.matrix(vec[, sources])
	
	
	#simulate source data
	scores <- matrix(nrow = nr, ncol = length(sources))
	for(i in 1 : length(sources)) {
		#get means/sd on logged scale
		m1 <- log(cm[i]^2/ sqrt(cm[i]^2 + sd[i]^2))
		s1 <- sqrt(log(sd[i]^2 / cm[i]^2 + 1))
		
		#generate sources
		scores[, i] <- rlnorm(nr, meanlog = m1, sdlog = s1)
	}
	
	#generate data
	dattemp <- as.matrix(scores) %*% t(vec2)
	
	#add lognormal errors
	n <- nrow(dattemp) * ncol(dattemp)
	errs <- matrix(rlnorm(n, sd = 0.01), nrow = nrow(dattemp))
	dattemp <- dattemp * errs
	
	
	list(dattemp, scores)
}





####
# function to generate data
# nmons is number of monitors
# reps is whether sources change by monitor
# ndays is number of observations
# PCs is varimax rotated PCs
# keeps is list of which sources are present at each monitor
# cms is vector of means corresponding to sources (vary by monitor)
# sds is vector of sds corresponding to sources
# unequal is number of monitors for each subregion if unequal
# days is number of days for each monitor if days
outerdat <- function(nmons, reps, ndays = 1000, PCs, keeps, 
	cms, sds, unequal = NULL, days = NULL) {
		
	datnew <- list()
	l <- 1
	
	#for each monitor
	for(i in 1 : nmons) {
		
		#determine number of observations
		if(!is.null(days)) {
			ndays = days[i]
		}
		
		#create data and save
		temp <- createdat(ndays, PCs, keeps[[l]], cms, sds)
		datnew[[i]] <- temp[[1]]
		
		#update subregion
		reps1 <- ifelse(!is.null(unequal), unequal, reps)
		l <- switchfun(reps1, i, l)
	}
	
	
	datnew
}




#######
# Function to run simulation for one set of data
# names is vector of source names 
#	(e.g. c("traffic", "fireworks", "soil")) corresponding to PCs
# nmons is number of monitors
# reps is number of monitors per subregion
# ndays is number of observations
# PCs is positive part of PC from sample data
# keeps is share info for creating data
# cms is vector of lognormal means for sources
# sds is vector of lognormal sds for sources
# unequal is vector of numbers to switch subregions
# days is vector of days for each monitor
# threli is threshold for number of sources
# thresang is threshold for share angle cutoff 
outerSIM <- function(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, unequal = NULL, days = NULL, threli = 1,
		thresang = pi/4) {
	
	#create dataset
	data <- outerdat(nmons, reps, ndays, PCs, keeps, 
		cms, sds, unequal, days)
		
	#perform SHARE
	dms <- domatchSIM(restrict.data = data, threli = threli, 
		usedata = T, thresang = thresang)
	share <- dms[[3]]
	reg <- dms[[1]][[1]]
	
		
	#perform mAPCA
	mapca <- as.matrix(spatial.apca(dat = data, lim = 50)[[6]][[2]]$load)

	
	#match share to PCs
	match1 <- as.vector(solve_LSAP(calc.dists(PCs, reg)))
	regnames1 <- names[match1]
	
	
	#match mAPCA to PCs
	match1 <- as.vector(solve_LSAP(calc.dists(PCs, mapca)))
	regnames2 <- names[match1]




	l <- 1
	matches <- matrix(nrow = length(share), ncol = 2)
	lens <- matrix(nrow = length(share), ncol = 2)

	#for each monitor
	for(i in 1 : nmons) {
		
		#get names
		s <- list()
		s[[1]] <- names[keeps[[l]]]
		s[[2]] <- regnames1[share[[i]]]
		s[[3]] <- regnames2
		
		
		#compare SHARE/mAPCA to PCs
		l1 <- sapply(s, length)
		for(j in 2 : 3) {
			s2 <- s
			#make same length
			if(l1[1] != l1[j]) {
				c1 <- which.min(l1[c(1, j)])
				c1 <- ifelse(c1 == 1, 1, j)
				seqs <- paste0("no", seq(1, 100))
				s2[[c1]] <- c(s2[[c1]], sample(seqs, abs(l1[j] - l1[1])))
			}
			
			#how many match
			matches[i, j - 1] <- length(which(s2[[1]] %in% s2[[j]]))
			#how many potential matches (max)
			lens[i, j - 1] <- sapply(s2, length)[1]

		}
		
		#update subregion
		reps1 <- ifelse(!is.null(unequal), unequal, reps)
		l <- switchfun(reps1, i, l)
	}
	
	#summarize
	out <- colSums(matches) / colSums(lens)

	
}





#####
# Function to run simulation across iterations
# nsims is number of iterations, remaining defined above
multsims <- function(nsims, names, nmons, 
	reps, ndays, PCs, keeps, cms, sds, 
	unequal = NULL, days = NULL, 
	threli = 1, thresang = pi/4) {
		
		
	outs <- matrix(nrow = nsims, ncol = 2)
	
	#for each simulation...
	for(i in 1 : nsims) {
		
		outs[i, ] <- outerSIM(names, nmons, 
			reps, ndays, PCs, keeps, 
			cms, sds, unequal, days, 
			threli, thresang)
	}
	colnames(outs) <- c("SHARE", "mAPCA")
	
	list(outs, colMeans(outs))
	
}