#' Share simulations
#'
#' \code{multsims} Performs share simulation for multiple datasets
#'
#' This is a function to compare results from SHARE and mAPCA for 
#' multiple simulated datasets
#'
#' @param nsims number of iterations
#' @param names vector of source names (e.g. c("traffic", "fireworks", "soil")) corresponding to PCs
#' @param nmons number of monitors
#' @param reps number of monitors per subregion
#' @param ndays number of observations
#' @param PCs positive part of PC from sample data
#' @param keeps share info for creating data
#' @param cms vector of lognormal means for sources
#' @param sds vector of lognormal sds for sources
#' @param unequal vector of numbers to switch subregions
#' @param days vector of days for each monitor
#' @param cut cutoff for eigenvalues (see nmsource), default is 1.
#' @param thres threshold for share angle cutoff 
#' @param prnt Print simulation iteration (default = F)
#' @export
multsims <- function(nsims, names, nmons, 
                     reps, ndays, PCs, keeps, cms, sds, 
                     unequal = NULL, days = NULL, 
                     cut = 1, thres = pi/4, prnt = F, sderr = 0.01) {
    
    
    outs <- matrix(nrow = nsims, ncol = 2)
    extras <- matrix(nrow = nsims, ncol = 2)
    
    #for each simulation...
    for(i in 1 : nsims) {
        if(prnt) {print(i)}
        temp <- outerSIM(names, nmons, 
                              reps, ndays, PCs, keeps, 
                              cms, sds, unequal, days, 
                              cut, thres, sderr = sderr)
        outs[i, ] <- temp[["match"]]  
        extras[i, ] <- temp[["extra"]]                    
    }
    colnames(outs) <- c("SHARE", "mAPCA")
    colnames(extras) <- c("SHARE", "mAPCA")
    
    list(fulloutput = outs, summary = colMeans(outs), extra = extras)
    
}






#' \code{outerSIM} Performs share simulation for one dataset
#'
#' This is a function to compare results from SHARE and mAPCA for 
#' one simulated dataset
#'
#' @param names vector of source names (e.g. c("traffic", "fireworks", "soil")) corresponding to PCs
#' @param nmons number of monitors
#' @param reps number of monitors per subregion
#' @param ndays number of observations
#' @param PCs positive part of PC from sample data
#' @param keeps share info for creating data
#' @param cms vector of lognormal means for sources
#' @param sds vector of lognormal sds for sources
#' @param unequal vector of numbers to switch subregions
#' @param days vector of days for each monitor
#' @param cut cutoff for eigenvalues (see nmsource), default is 1.
#' @param thres threshold for share angle cutoff 
outerSIM <- function(names, nmons, reps, ndays, PCs, keeps, 
                     cms, sds, unequal = NULL, days = NULL, cut = 1,
                     thres = pi/4, sderr = 0.01) {
    
    #create dataset
    data <- outerdat(nmons, reps, ndays, PCs, keeps, 
                     cms, sds, unequal, days, sderr = sderr)

    #perform SHARE
    share1 <- share(data, cut = cut, thres = thres)
    share <- share1$share
    major.sig <- share1$major.sig
    
    
    #perform mAPCA
    mapca <- as.matrix(mAPCA(dat = data, lim = 50)[["apca"]][["vmax"]]$load)
    
    
    #match share to PCs
    match1 <- solveLSAP_nc(PCs, major.sig) 
    regnames1 <- names[match1]
    
    
    #match mAPCA to PCs
    match1 <- solveLSAP_nc(PCs, mapca) 
    regnames2 <- names[match1]
    
    
    l <- 1
    matches <- matrix(nrow = length(share), ncol = 2)
    lens <- matrix(nrow = length(share), ncol = 2)
    overid <- matrix(nrow = length(share), ncol = 2)
    
    #for each monitor
    for(i in 1 : nmons) {
        
        #get names for truth, SHARE, mAPCA
        s <- list()
        s[[1]] <- names[keeps[[l]]]
        s[[2]] <- regnames1[share[[i]]]
        s[[3]] <- regnames2
        
        
        #compare SHARE/mAPCA to PCs
        l1 <- sapply(s, length)
        
        #for SHARE and mAPCA
        for(j in 2 : 3) {
            s2 <- s
            #If length truth != length of SHARE/mAPCA
            if(l1[1] != l1[j]) {
                
                #find smaller length
                c1 <- which.min(l1[c(1, j)])
                c1 <- ifelse(c1 == 1, 1, j)
                #add source names to make up deficit
                seqs <- paste0("no", seq(1, 100))
                s2[[c1]] <- c(s2[[c1]], sample(seqs, abs(l1[j] - l1[1])))
            }
            
            #how many match
            matches[i, j - 1] <- length(which(s2[[1]] %in% s2[[j]]))
            #how many potential matches (max)
            lens[i, j - 1] <- sapply(s2, length)[1]
            
            overid[i, j - 1] <- length(which(is.na(s2[[j]])))
            
        }
        
        #update subregion
        reps1 <- ifelse(!is.null(unequal), unequal, reps)
        
        #Iterate through truth as necessary
        l <- switchfun(reps1, i, l)
    }
    
    #summarize
    out <- colSums(matches) / colSums(lens)
    
    list(match = out, extra = colMeans(overid))
    
}




#### # Function to get match 
solveLSAP_nc <- function(pcs, major) {
	ang1 <- angle(pcs, major)
	
	nr <- nrow(ang1)
    ns <-  nrow(ang1) - ncol(ang1)
    if(ns > 0) {
		adds <- matrix(rep(100, nr * ns), nrow = nr, ncol = ns)	
		ang1 <- cbind(ang1, adds)
	}
	
	as.vector(solve_LSAP(ang1))

}


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
createdat <- function(nr, vec, sources, cm, sd, sderr = 0.01) {
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
	errs <- matrix(rlnorm(n, sd = sderr), nrow = nrow(dattemp))
	dattemp <- dattemp * errs
	
    #add date
    dates <- as.Date(seq(1, nrow(dattemp)), origin = "1970-01-01")
    dattemp <- data.frame(dates, dattemp)
	
	list(data = dattemp, source = scores)
}





####
# function to generate simulated source data
####
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
	cms, sds, unequal = NULL, days = NULL, sourceout = NULL,
	sderr = 0.01) {
		
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
		temp <- createdat(ndays, PCs, keeps[[l]], cms, sds, sderr = sderr)
		datnew[[i]] <- temp$data
		source[[i]] <- temp$source
		
		#update subregion
		reps1 <- ifelse(!is.null(unequal), unequal, reps)
		l <- switchfun(reps1, i, l)
	}
	
	
	out <- datnew
	
	if(!is.null(sourceout)) {
		out <- list(data = datnew, source = source)
	}
	out
}








##################
# Functions for hospitalization simulation 

#### Function to simulate hospitalization data
# sources is matrix of source concentrations
# betas is vector of associations
# share is list of betas to use for each monitor
# int is intercept
hospdat <- function(sources, betas, share, int = 5) {
	
	y <- list()
	#for each monitor
	for(i in 1 : length(sources)) {
		#get mean for poisson
		lincomb <- int + t(as.matrix(betas[share[[i]]])) %*% t(sources[[i]])
		mean <- exp(lincomb)
		
		#simulate y data
		y[[i]] <- rpois(nrow(sources[[i]]), mean)
	}
	
	y
}





#' \code{outerSIMhosp} Performs share simulation for health effects for one dataset
#'
#' This is a function to compare results from SHARE and mAPCA for 
#' one simulated dataset
#'
#' @param names vector of source names (e.g. c("traffic", "fireworks", "soil")) corresponding to PCs
#' @param nmons number of monitors
#' @param reps number of monitors per subregion
#' @param ndays number of observations
#' @param PCs positive part of PC from sample data
#' @param keeps share info for creating data
#' @param cms vector of lognormal means for sources
#' @param sds vector of lognormal sds for sources
#' @param etas vector of percent increases for health effects
#' @param unequal vector of numbers to switch subregions
#' @param days vector of days for each monitor
#' @param cut cutoff for eigenvalues (see nmsource), default is 1.
#' @param thres threshold for share angle cutoff 
#' @param int baseline hosp rate for health effect regression
#' @param prnt Print simulation iteration (default = F)
#' @export
outerSIMhosp <- function(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, etas, unequal = NULL, days = NULL, cut = 1,
		thres = pi/4, int = 5, prnt = F, sderr = 0.01) {
	
	#create dataset
	temp <- outerdat(nmons, reps, ndays, PCs, keeps, 
		cms, sds, unequal, days, sourceout = T, sderr = sderr)
	data <- temp$data
	source <- temp$source	
	
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
	share1 <- share(data = data, cut = cut, thres = thres)
	share <- share1$share
	sources <- share1$Sources
	reg <- share1$major.sig

		
	#perform mAPCA
	mapca <- mAPCA(data = data, lim = 150)$apca
	mapcasource1 <- mapca$conc
	mapca <- as.matrix(mapca[["vmax"]]$load)
                       
	#get list of mapca results by monitor
    cn <- which(colnames(mapcasource1) == "mons")
    mons <- unique(mapcasource1$mons)
    mapcasource <- list()
    for(i in 1 : length(mons)) {
        mapcasource[[i]] <- mapcasource1[which(mapcasource1$mons == mons[i]), -cn]
    }

	#apply APCA
	apca <- list()
	for(i in 1 : length(data)) {
		nf1 <- length(share[[i]])
		apca[[i]] <- apca(data = data[[i]], nsources = nf1)$conc
	}


	#get hospiatlization data
	betas <- 1/iqrs * log(etas / 100 + 1)
	y <- hospdat(source, betas, shareT, int)

	#match share to PCs
	match1 <- as.vector(solve_LSAP(angle(PCs, reg)))
	regnames1 <- names[match1]

	
	#match mAPCA to PCs
	match1 <- as.vector(solve_LSAP(angle(PCs, mapca)))
	regnames2 <- names[match1]

    #get regional health effects
	tlnmAPCA <- tlnout(regnames2, y, mapcasource, "mapca")
	tlnAPCA <- tlnout(regnames1, y, apca, "apca", share)
	tlntruth <- tlnout(names, y, source, "truth", shareT)
	

	#match results
	tlnAPCA <- tlnAPCA[match(rownames(tlntruth), rownames(tlnAPCA)), ]
	tlnmAPCA <- tlnmAPCA[match(rownames(tlntruth), rownames(tlnmAPCA)), ]
	
	#get IQR increase
    tlnAPCApi <- matrix(nrow = nrow(tlnAPCA), ncol = 3)
	tlnmAPCApi <- matrix(nrow = nrow(tlnmAPCA), ncol = 3)
	tlntruthpi <- matrix(nrow = nrow(tlntruth), ncol = 3)
	colnames(tlnAPCApi) <- c("est", "lb", "ub")
	colnames(tlnmAPCApi) <- c("est", "lb", "ub")
	colnames(tlntruthpi) <- c("est", "lb", "ub")
	for(i in 1 : nrow(tlntruth)) {
		tlnAPCApi[i, ] <- percincCI(tlnAPCA[i, ], scale = iqrs[i])
		tlnmAPCApi[i, ] <- percincCI(tlnmAPCA[i, ], scale = iqrs[i])
		tlntruthpi[i, ] <- percincCI(tlntruth[i, ], scale = iqrs[i])
	}
    
    out1 <- list(truth = tlntruth, APCA = tlnAPCA, mAPCA = tlnmAPCA)
    out2 <- list(truth = tlntruthpi, APCA = tlnAPCApi, mAPCA = tlnmAPCApi)
	out <- list(regcoef = out1, percinc = out2, iqrs = iqrs)
	out
}



#### function to get CI for percinc
percincCI <- function(x, scale = 1) {
	lb <- x[1] - 1.96 * x[2]
	ub <- x[1] + 1.96 * x[2]
	cis <- c(x[1], lb, ub)
	percinc(cis, scale = scale)
}


#### function to find match for each row
#matchmat is matrix of 1/0 matches from matchfun
whichCS <- function(matchmat) {
    mins <- apply(matchmat, 1, function(x) suppressWarnings(min(which(x > 0))))
    mins
    
}


######
# Function to combine regression results
# unsources is vector of unique sources
# y is list of health outcome counts
# sourceconc is list of source concentrations
# type is type of source apportionment
# share is which sources for truth/apca
tlnout <- function(unsources, y, sourceconc, type, share = NULL, print = F) {
		
		
	#for each source
	glmout <- matrix(nrow = length(unsources), ncol = 2)
	colnames(glmout) <- c("est", "se")
	rownames(glmout) <- unsources
	
    #get results for each source
	for(i in 1 : length(unsources)) {
	
		glm <- matrix(ncol = 2, nrow = length(y))
		l <- 1
        
        #for each monitor
		for(j in 1 : length(y)) {
			
            #get appropriate source for apca/truth
			if(type %in% c("truth", "apca")) {
				share1 <- share[[j]]
				
				wh1 <- which(share1 == i)
				if(length(wh1) > 0) {
					sourceconc1 <- sourceconc[[j]][, wh1]
				}else{
					
					sourceconc1 <- NULL
				}
            #take original order for mAPCA
			}else{
				sourceconc1 <- sourceconc[[j]][, i]
				
			}
			
            #perform health effect regression
			if(!is.null(sourceconc1)) {
				temp <- try(glm(y[[j]] ~ sourceconc1, 
					family = "poisson"))
				if(class(temp)[1] != "try-error" & !(is.na(temp$coef[2]))) {

					glm[j, ] <- summary(temp)$coef[-1, c(1, 2)]
				} 
			}
			

		
		}#end loop over monitors
		
        #combine results across monitors
	    glm <- glm[complete.cases(glm),]
    	if(nrow(glm) > 1) {
		
	    	temp <- tlniseC(Y = glm[, 1], V = glm[, 2]^2, prnt = print, brief = 2)
            if(temp$converge == "no") {
                temp <- tlniseC(Y = glm[, 1], V = glm[, 2]^2, 
                    prnt = print, brief = 2, maxiter = 5000)
            }
		    glmout[i, ] <- temp$gamma[1:2]
	    }else if(nrow(glm) == 1){
	    	glmout[i, ] <- glm
	    }
	
	
	}#end loop over sources
	
	glmout
}






#' \code{reorderout} Get outcome
#'
#' @param out regression outcomes
#' @param sources sources
#' @param rn rownames
#' @export
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





#' \code{gethospsim} Get outcome
#'
#' @param outmult regression outcomes
#' @param iqrs iqrs
#' @export
gethospsim <- function(outmult, iqrs) {
    outres <- list()
    outresPI <- list()
    #for each method
    for(i in 1 : 3) {
        
        outres[[i]] <- matrix(nrow = nrow(outmult[[1]][[i]]), ncol = 2)
        rownames(outres[[i]]) <- rownames(outmult[[1]][[i]])
        
        outresPI[[i]] <- matrix(nrow = nrow(outmult[[1]][[i]]), ncol = 3)
        colnames(outresPI[[i]]) <- c("est", "lb", "ub")
        
        outsel <- sapply(outmult, function(x) x[[i]], simplify = F)
        for(j in 1 : nrow(outsel[[1]])) {
            outsel1 <- sapply(outsel, function(x) x[j, ])
            mn1 <- mean(outsel1[1, ], na.rm = T, trim = 0.1)
            within <- mean(outsel1[2, ]^2, na.rm = T, trim = 0.1)
            
            between1 <- (outsel1[1, ] - mn1)^2
            between1 <- mean(between1, na.rm = T, trim = 0.1)
            between <- ncol(outsel1)/(ncol(outsel1) - 1) * between1
            
            outres[[i]][j, 1] <- mn1
            sd1 <- sqrt(within + (1 + 1/ncol(outsel1)) * between)
            outres[[i]][j, 2]  <- sd1
            
            lb <- mn1 - 1.96 * sd1
            ub <- mn1 + 1.96 * sd1
            
            outresPI[[i]][j, 1] <- percinc(outres[[i]][j, 1], iqrs[j]) 
            outresPI[[i]][j, 2] <- percinc(lb, iqrs[j]) 
            outresPI[[i]][j, 3] <- percinc(ub, iqrs[j])  
        }

        

        
    }
    list(regcoef = outres, percinc = outresPI)
    
}




#' \code{msefun} Get MSE for results
#'
#' @param percinc results for percent increase
#' @param etas vector of percent increases for health effects
#' @param rn rownames
#' @export
msefun <- function(percinc, etas, rn)  {
    
    #get mse
    mses1 <- matrix(nrow = nrow(percinc[[1]][[1]]), ncol = 3)
    rownames(mses1) <- rn
    colnames(mses1) <- names(percinc[[1]])

    for(i in 1 : 3) {
    		#select estimates
        saps1 <- sapply(percinc, function(x) x[[i]][, 1])
        #difference from truth
        temp1 <- (sweep(saps1, 1, etas))^2
        mses1[, i] <- apply(temp1, 1, mean, na.rm = T, trim = 0.1)
    }
    mses1
    
}