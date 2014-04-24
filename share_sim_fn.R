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




hospdat2 <- function(sources, betas, share, names, int = 5) {
	
	y <- list()
	#for each monitor
	for(i in 1 : length(sources)) {
		y[[i]] <- matrix(nrow = nrow(sources[[i]]), ncol = ncol(sources[[i]]))
		colnames(y[[i]]) <- names[share[[i]]]
		for(j in 1 : ncol(sources[[i]])) {
			#get mean for poisson
			lincomb <- 5 + betas[share[[i]][j]] * sources[[i]][, j]
			mean <- exp(lincomb)
			
			#simulate y data
			y[[i]][, j] <- rpois(nrow(sources[[i]]), mean)
			
		}
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
	sources <- dms[[4]]
	reg <- dms[[1]][[1]]
	match1 <- as.vector(solve_LSAP(calc.dists(reg, PCs)))

	
	
		
	#perform mAPCA
	mapca <- spatial.apca(dat = data, lim = 50)
	mapcasource <- mapca[[1]][[2]]
	mapca <- as.matrix(mapca[[6]][[2]]$load)
	
	
	#apply APCA
	apca <- list()
	betas <- list()
	vs <- list()
	for(i in 1 : length(data)) {
		nf1 <- length(share[[i]])
		temp <- abspca(data[[i]], nfactors = nf1)		
		apca[[i]] <- temp[[1]]
		betas[[i]] <- temp[[6]]
		vs[[i]] <- temp[[7]]
	}
	
	apca2 <- fixerror2(apca, data, betas, vs, share, sources)


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
	tlnAPCA2 <- tlnout(regnames1, y, apca2, "apca", share)
	tlntruth <- tlnout(names, y, source, "truth", shareT)
	

	#match results
	tlnAPCA <- tlnAPCA[match(rownames(tlntruth), rownames(tlnAPCA)), ]
	tlnAPCA2 <- tlnAPCA2[match(rownames(tlntruth), rownames(tlnAPCA2)), ]
	tlnmAPCA <- tlnmAPCA[match(rownames(tlntruth), rownames(tlnmAPCA)), ]
	
	#get IQR increase
	for(i in 1 : nrow(tlntruth)) {
		tlnAPCA[i, ] <- percinc(tlnAPCA[i, ], scale = iqrs[i])
		tlnAPCA2[i, ] <- percinc(tlnAPCA2[i, ], scale = iqrs[i])
		tlnmAPCA[i, ] <- percinc(tlnmAPCA[i, ], scale = iqrs[i])
		tlntruth[i, ] <- percinc(tlntruth[i, ], scale = iqrs[i])
	}
	out <- list(tlntruth, tlnAPCA, tlnAPCA2, tlnmAPCA)
	names(out) <- c("truth", "APCA", "APCA2", "mAPCA")
	out
}


fixerror <- function(apca, share, sources, sim = T) {
	
	apcaout <- list()
	
	#create list of sources
	for(i in 1 : length(sources)) {
		shares <- sapply(share, function(x) {	
			ifelse(i %in% x, which(x == i), 0)})
			
		temp <- 0	
		dates <- as.Date("1970-01-01")
		mons <- 0
		for(j in 1 : length(apca)) {
			if(shares[j] != 0) {
				temp2 <- apca[[j]][, shares[j]]
				if(sim == T) {
					dates1 <- as.Date(seq(1, nrow(apca[[j]])), 
					origin = "1970-01-01")
				}else{
					dates1 <- apca[[j]][, 1]
				}
				
				
				dates <- c(dates, dates1)
				temp <- c(temp, temp2)
				mons <- c(mons, rep(j, nrow(apca[[j]])))
			}
		}
		temp <- data.frame(dates, temp, mons)
		temp <- temp[-1, ]
		
		#fixed mixed model
		# lm1 <- lme(temp ~ 1,random = list(mons=~1, dates=~1), data = temp)
		# random = ~1 | mons + dates, data = temp)
		# lm1 <- lme(temp ~ 1, random = pdBlocked(list(pdSymm(~mons-1), pdSymm(~dates-1))), data = temp)
		lm1 <- lmer(temp ~ (1|dates) + (1|mons), data = temp)
		
		
		apcaout[[i]] <- temp
		
	}
	
	
	
}




#this is for each source
getsigma2 <- function(dat, betas, vs, shares) {
	
	#for each data
	all <- c(0, 0, 0, 0)
	outs <- vector(, length = length(shares))
	
	#for each monitor
	for(i in 1 : length(dat)) {
		# all <- c(0, 0, 0, 0)
		
		#if source at monitor
		if(shares[i] > 0) {
		
		
		#create iterations
		for(j in 1 : length(shares)) {
			
			#print(c(i, j))
			temp <- dat[[i]]
			
			#if source at monitor
			if(shares[j] > 0) {
				
				
				vorg <- vs[[i]][, shares[i]]
				vl <- vs[[j]][, shares[j]]
				if(sum(vorg * vl) < 0) {
					#print("rev")
					vl <- -vl
				}
				
				# names <- names(vl)
				names <- intersect(names(vl), names(vorg))
				vl <- vl[names]
				vorg <- vorg[names]
												
				if(length(names) != ncol(temp)) {
					temp <- temp[, names]					
				}				
				
				sd1 <-  apply(temp, 2, sd)
				wh0 <- which(sd1 == 0)
				if(length(wh0) > 0) {
					# browser()
					sd1 <- sd1[-wh0]
					temp <- temp[, -wh0]
					vl <- vl[-wh0]
				}
				
				sds <- diag(sd1)

				bl <- betas[[i]][shares[i]]
				# vlbl <- vl * bl
				vlbl <- vorg * bl
				
				#original source concentration
				f1o <-  temp %*% sds %*% solve(cor(temp)) %*% vlbl
				#new A
				a1 <- temp %*% sds %*% solve(cor(temp)) %*% vl

				#scaling param
				bet <- lm(f1o~ a1)$coef[-1]
				
				f1 <- a1 * bet

				if(length(which(is.na(f1))) != 0) { browser()}
				# f1 <- temp %*% vlbl
				others1 <- matrix(rep(c(i, j), length(f1)), 
					byrow= T, ncol = 2)
					

				time <- seq(1, length(f1))
				
				all <- rbind(all, cbind(f1, others1, time))
			}#end test share
			
			
		}#end loop over shares


		}#end test monitor
	}#end loop over monitor
	all <- all[-1, ]
	all <- data.frame(all)
	# all[, 3] <- paste0(all[, 2], all[, 4])
	
	colnames(all) <- c("conc", "mon", "iter", "time")
	
	all2 <- data.frame(all, paste0(all$mon, ":", all$time))
	colnames(all2) <- c(colnames(all), "montime")
	
	all2$mon <- factor(all2$mon)
	all2$time <- factor(all2$time)
	all2$iter <- factor(all2$iter)
	
	all2 <- all2[complete.cases(all2), ]
	

	outs <- vector()
	unmon <- as.numeric(as.character(unique(all2$mon)))
	
	apca <- list()
	for(k in 1 : length(dat)) {
		#print(k)
		if(k %in% unmon) {
			temp <- all2[which(all2$mon == k), ]
			apca[[k]] <- tapply(temp$conc, temp$time, mean, na.rm = T)
			# dat1 <- all2[which(all2$mon == k), ]
			# lm1  <- lmer(conc ~ (1 | iter) + 
				# (1 | time), data = dat1)
			# outs[k]	<- summary(lm1)$varcor$iter[[1]]
		}
	}
	# outs
	apca	
	
}

##### # code to get sigma, adjusted source est 
fixerror2 <- function(apca, dat, betas, vs, share, sources, sim = T) {
	
	apca1 <- list()
	apcaout <- list()
	#for each source
	for(i in 1 : length(sources)) {
		shares <- sapply(share, function(x) {	
			ifelse(i %in% x, which(x == i), 0)})
		
		#get vl/bl info
		
		#estimate sigma
		# sigma2 <- getsigma2(dat, betas, vs, shares)
		apca1[[i]] <- getsigma2(dat, betas, vs, shares)
		# print(sigma2)
		# sigma2 <- 0.1
		
# # 		#get empirical bayes estimate
		# for(j in 1 : length(apca)) {
			# if(i == 1) {
				# apcaout[[j]] <- apca[[j]]
			# }
			# if(shares[j] > 0) {
				# temp <- apca[[j]][, shares[j]]
				# mui <- mean(temp)
				
				# if(length(sigma2) > 1) {
					# sigma2U <- sigma2[j]
				# }else{
					# sigma2U <- sigma2
					# }
				# tau2i <- max(c(0, var(temp) - sigma2U))
				
				# B <- sigma2U / (sigma2U + tau2i)
				# # print(c(sigma2U, tau2i, B))
				# apcaout[[j]][, shares[j]] <- (1 - B) * temp + B * mui
			# }
		# }
		
		
	}
	
	
	#fix order 
	for(i in 1 : length(apca)) {
		share1 <- share[[i]]
		apcaout[[i]] <- matrix(nrow = nrow(apca[[i]]), ncol = length(share1))
		for(j in 1 : length(share1)) {
			apcaout[[i]][, j] <- apca1[[share1[j]]][[i]]
		}
	}
	
	
	apcaout
	
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
				temp <- try(glm(y[[j]] ~ sourceconc1, 
					family = "poisson"))
				if(class(temp)[1] != "try-error" & !(is.na(temp$coef[2]))) {
					#print(j)
					glm[j, ] <- summary(temp)$coef[-1, c(1, 2)]
				} #else{ browser()}
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
	n <- length(out[[1]]) - length(types)
	if(n > 0) {
		types <- c(types[1:2], paste0("SHARE", seq(1, n)), "mAPCA")
	}
	info <- rep("0", 3)
	
	dat <- rep(0, 4)
	
	for(j in 1 : length(out)) {
		sim <- c("A", "B", "C", "D", "E")[j]
		for(i in 1 : length(out[[j]])) {
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