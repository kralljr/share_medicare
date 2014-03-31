
med.path <- "/dexter/disk2/rpeng/medicare/MCAPS1999-2010"
#load library
library(tsModel, lib.loc = "/home/bst/student/jkrall")
library(splines)

#load APCA, mAPCA results
load("sourceconc_medicare.RData")


# regression formula
# # form1 <- "factor(agecat) + factor(dow) + ns(tmpd, df = 3) + 
	# ns(Lag(tmpd, 1, group = agecat),df = 3) + ns(date, df = 8 * "
form1 <- "factor(agecat) + factor(dow) + ns(tmpd, df = 6) + ns(tmp3, df = 6)  + ns(dptp, df = 3) + ns(dpt3, df = 3) + ns(date, df = 8 * "


#restrict data to what need
rest.data <- function(dat) {
	
	#get outcome data
	datout <- dat[, c("date", "agecat", "cardio", "resp", "denom")]
	#sum for same age/date
	datout <- aggregate(datout[, -c(1, 2)], by = list(datout$date, 
		datout$agecat), FUN = sum)
	colnames(datout) <- c("date", "agecat", "cardio", "resp", "denom")
	
	
	
	#get temperature data
	dattemp <- dat[, c("fips", "date", "tmpd", "dptp", 
		"dow")]
	dattemp <- dattemp[-which(duplicated(dattemp)), ]
	#get running means
	tmp3 <- runMean(dattemp$tmpd, lags = 1 : 3)
	dpt3 <-  runMean(dattemp$dptp, lags = 1 : 3)
	dattemp <- data.frame(dattemp, tmp3, dpt3)
	
	dat <- merge(dattemp, datout, by = "date", all.x = T, all.y = T)
}




get.hosp <- function(dat, sources, lag, formu = form1, outcome = "cardio") {
	
	#get running temp, format data	
	dat <- rest.data(dat)
	
	#merge with sources
	merged <- merge(dat, sources, all.x = T, 
		by.x = "date", by.y = "Date")

	#get lag info	
	#which columns are sources
	whSource <- which(substr(colnames(merged), 1, 4) == "sour")
	lagSource <- vector(, length = nrow(merged))
	
	temp2 <- rep(NA, nrow(merged))
	#for each source, #lag by agecat
	for(i in 1 : length(whSource)) {

		sour <- merged[, whSource[i]]	
		#for each source with at least 1 day of data
		if(length(which(!is.na(sour))) > 0) {
			temp <- Lag(sour, 
				k = lag, group = merged[, "agecat"])
			lagSource <- cbind(lagSource, temp)
		#else all NA
		}else{
			lagSource <- cbind(lagSource, temp2)
			}
	}
	lagSource <- lagSource[, -1]
	
	#add in lagged sources
	merged[, whSource] <- lagSource
	
	#number of years of data (should be 12)
	years <- length(unique(substr(merged$date, 1, 4)))
	
	#get formula
	formUSE1 <- paste0(formu, years)
	
	#set up estimates
	ests <- matrix(nrow = length(whSource), ncol = 2)
	#for each source
	for(l in 1 : length(whSource)) {
		# print(c("l", l))
		covar1 <- paste(colnames(merged)[whSource[l]], collapse = "+")
		formUSE <- paste0(outcome, " ~", covar1, "+", formUSE1, ")")
		
		#run model
		glm1 <- try(glm(formula = eval(formUSE), 
			data = merged, family = "quasipoisson",
			offset = log(denom)), silent = T)
			
		#save results
		if(class(glm1) != "try-error") {
			out1 <- summary(glm1)$coef
			whS <- which(substr(rownames(out1), 1, 4) == "sour")
			out <- out1[whS, c(1, 2)]
		}else{
			out <- c(NA, NA)
		}
		ests[l, ] <- out

	}

	rownames(ests) <- colnames(merged)[whSource]
	colnames(ests) <- c("est", "se")
	ests
}





set.seed(8996)
mortsmAPCA <- list()
mortsAPCA <- list()
lag <- 1
i <- 1

# for(lag in 1 : 3) {
	mortsmAPCA[[lag]] <- list()
	mortsAPCA[[lag]] <- list()
	
	#for each county
	# for(i in 1 : length(mAPCA1[[1]])) {
		print(i)
		
		fips <- names(mAPCA1[[1]])[i]
		dat <- try(readRDS(file.path(med.path, paste0(fips, ".rds"))))
		
		if(class(dat) == "try-error"){
			print(fips)
			stop("Error: data doesn't exist")
		}
		mortsmAPCA[[lag]][[i]] <- get.hosp(dat = dat, 
			sources = mAPCA1[[1]][[i]], lag = lag - 1)
		mortsAPCA[[lag]][[i]] <- get.hosp(dat = dat, 
			sources = APCA1[[1]][[i]], lag = lag - 1)
	# }
# }

# save(mortsmAPCA, mortsAPCA, file = "mort_sources_medicare.RData")
