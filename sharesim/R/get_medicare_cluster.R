


#' \code{rest.data} Clean medicare data
#' 
#' @param dat medicare data for one fips county
#' @export
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

