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
	cms, sds, unequal = NULL, days = NULL, sourceout = NULL) {
		
	datnew <- list()
	source <- list()
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
		source[[i]] <- temp[[2]]
		
		#update subregion
		reps1 <- ifelse(!is.null(unequal), unequal, reps)
		l <- switchfun(reps1, i, l)
	}
	
	
	out <- datnew
	
	if(!is.null(sourceout)) {
		out <- list(datnew, source)
	}
	out
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



#### Function to simulate hospitalization data
# nreps is number of simulations
# sources is matrix of source concentrations
# betas is vector of associations
# int is intercept
hospdat <- function(sources, betas, share, int = 5) {
	
	y <- list()
	#for each monitor
	for(i in 1 : length(sources)) {
		#get mean for poisson
		lincomb <- 5 + t(as.matrix(betas[share[[i]]])) %*% t(sources[[i]])
		mean <- exp(lincomb)
		
		#simulate y data
		y[[i]] <- rpois(nrow(sources[[i]]), mean)
	}
	
	y
}





#### Function to simulate mortality effects
outerSIMhosp <- function(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, etas, unequal = NULL, days = NULL, threli = 1,
		thresang = pi/4, int = 5) {
	
	#create dataset
	temp <- outerdat(nmons, reps, ndays, PCs, keeps, 
		cms, sds, unequal, days, sourceout = T)
	data <- temp[[1]]
	source <- temp[[2]]	
	
	#get all keeps
	shareT <- list()
	l <- 1
	for(i in 1 : length(data)) {
		shareT[[i]] <- keeps[[l]]
		
		reps1 <- ifelse(!is.null(unequal), unequal, reps)
		l <- switchfun(reps1, i, l)	
	}
	
	
	#get IQR
	iqrsALL <- lapply(source, function(x) apply(x, 2, IQR))
	iqrs <- vector()
	for(i in 1 : length(names)) {
		
		wh1 <- sapply(shareT, function(x) {
				ifelse(i %in% x, which(x == i), 0)})
		iqrs1 <- iqrsALL[which(wh1 > 0)]
		wh1 <- wh1[which(wh1 > 0)]	
		
		holds <-  vector()
		if(length(wh1))
		for(j in 1 : length(wh1)) {
			holds[j] <- iqrs1[[j]][wh1[j]]
		}	
				
		iqrs[i] <- median(holds)
		
	}
	names(iqrs) <- names
	
	
	
		
	#perform SHARE
	dms <- domatchSIM(restrict.data = data, threli = threli, 
		usedata = T, thresang = thresang)
	share <- dms[[3]]
	reg <- dms[[1]][[1]]
	
		
	#perform mAPCA
	mapca <- spatial.apca(dat = data, lim = 50)
	mapcasource <- mapca[[1]][[2]]
	mapca <- as.matrix(mapca[[6]][[2]]$load)
	
	
	#apply APCA
	apca <- list()
	for(i in 1 : length(data)) {
		nf1 <- length(share[[i]])
		temp <- abspca(data[[i]], nfactors = nf1)
		# wsds0 <- which(apply(data[[i]], 2, sd) == 0)
		# pc <- temp[[4]]
		# dat <- data[[i]]
		# if(length(wsds0) > 0) {
			# pc <- pc[-wsds0, ]
			# dat <- dat[, -wsds0]
		# }
		# vars <- diag(cov(dat %*% solve(cov(dat)) %*% pc))
		# vara <- vars / apply(temp[[1]], 2, var)
		# temp2 <- sweep(temp[[1]], 2, sqrt(vara), "*")
		# shifts <- (1 - vara) * apply(temp[[1]], 2, mean)
		# apca[[i]]  <- sweep(temp2, 2, shifts, "+")
		apca[[i]] <- temp[[1]]
	}


	#get hospiatlization data
	betas <- 1/iqrs * log(etas / 100 + 1)
	y <- hospdat(source, betas, shareT, int)

	
	#match share to PCs
	match1 <- as.vector(solve_LSAP(calc.dists(PCs, reg)))
	regnames1 <- names[match1]
	
	
	#match mAPCA to PCs
	match1 <- as.vector(solve_LSAP(calc.dists(PCs, mapca)))
	regnames2 <- names[match1]


	tlnmAPCA <- tlnout(regnames2, y, mapcasource, "mapca")
	tlnAPCA <- tlnout(regnames1, y, apca, "apca", share)
	tlntruth <- tlnout(names, y, source, "truth", shareT)

	#match results
	tlnAPCA <- tlnAPCA[match(rownames(tlntruth), rownames(tlnAPCA)), ]
	tlnmAPCA <- tlnmAPCA[match(rownames(tlntruth), rownames(tlnmAPCA)), ]
	
	#get IQR increase
	for(i in 1 : nrow(tlntruth)) {
		tlnAPCA[i, ] <- percinc(tlnAPCA[i, ], scale = iqrs[i])
		tlnmAPCA[i, ] <- percinc(tlnmAPCA[i, ], scale = iqrs[i])
		tlntruth[i, ] <- percinc(tlntruth[i, ], scale = iqrs[i])
	}
	out <- list(tlntruth, tlnAPCA, tlnmAPCA)
	names(out) <- c("truth", "APCA", "mAPCA")
	out
}




tlnout <- function(unsources, y, sourceconc, type, share = NULL) {
		
		
	#for each source
	glmout <- matrix(nrow = length(unsources), ncol = 2)
	colnames(glmout) <- c("est", "se")
	rownames(glmout) <- unsources
	
	for(i in 1 : length(unsources)) {
	
		#for each monitor
		glm <- matrix(ncol = 2, nrow = length(y))
		l <- 1
		for(j in 1 : length(y)) {
		
			#get source conc
			if(type %in% c("truth", "apca")) {
				share1 <- share[[j]]
				
				wh1 <- which(share1 == i)
				if(length(wh1) > 0) {
					sourceconc1 <- sourceconc[[j]][, wh1]
				}else{
					
					sourceconc1 <- NULL
				}
			}else{
				sourceconc1 <- sourceconc[[j]][, i]
				
			}
			
			if(!is.null(sourceconc1)) {
				glm[j, ] <- summary(glm(y[[j]] ~ sourceconc1, 
					family = "poisson"))$coef[-1, c(1, 2)]
			}
			

		
		}#end loop over monitors
		
	glm <- glm[complete.cases(glm),]
	if(nrow(glm) > 1) {
		
		glmout[i, ] <- tlnise(Y = glm[, 1], V = glm[, 2]^2)$gamma[1:2]
		
	}else if(nrow(glm) == 1){
		glmout[i, ] <- glm
	}
	
	
	}#end loop over sources
	
	glmout
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




reorderout <- function(out, nr, sources) {
	types <- c("Known", "SHARE", "mAPCA")
	info <- rep("0", 3)
	
	dat <- rep(0, 4)
	
	for(j in 1 : length(out)) {
		sim <- c("A", "B", "C", "D", "E")[j]
		for(i in 1 : 3) {
			d1 <- out[[j]][[i]]
			lb <- d1[, 1] - 1.96 * d1[, 2]
			ub <- d1[, 1] + 1.96 * d1[, 2]
			
			dat1 <- cbind(d1, lb, ub)
			info1 <- cbind(rep(types[i], nr), sources, rep(sim, nr))
			
			dat <- rbind(dat, dat1)
			info <- rbind(info, info1)
		}
	}
		
	dat <- dat[-1, ]
	info <- info[-1, ]
	
	out <- data.frame(dat, info)
	colnames(out) <- c("est", "se", "lb", "ub", "Type", "Source", "Sim")
	out
}