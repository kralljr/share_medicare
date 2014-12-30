#####  File to get associations between hospitalizations and sources in east


rm(list = ls())
args <- commandArgs(TRUE)
lagi <- as.numeric(args[[1]])
season <- args[[2]]
print(c(lagi, season))



#load my libraries
library(share)
library(handles)
library(sharesim)
library(dplyr)
library(splines)




#######
# get constituent data
data(unmonlist)
load("speciation_medicare.RData")
data.rr <- lapply(datall, function(x) {
    x[, -c(2)]
})
names(data.rr) <- monsKEEP[, 3]
pm25 <- lapply(datall, function(x) x[, 2])


iswarm <- function(dat) {
	months <- as.numeric(substr(dat[, 1], 6, 7))
	ifelse(between(months, 4, 9), 1, 0)
}
warm <- lapply(data.rr, iswarm)
names1 <- names(data.rr)


seas <- ifelse(season == "w", 1, 0)

data.rr <- lapply(data.rr, function(x) {
	warm <- iswarm(x)
	x[warm == seas, ]
	
})


for(i in 1 : length(pm25)) {
	pm25[[i]] <- pm25[[i]][warm[[i]] == seas]
}

##########




# fix for less than 50
lens.keep <- which(sapply(pm25, length) >= 50)
data.rr <- data.rr[lens.keep]
pm25 <- pm25[lens.keep]



#get new list of monitors
names1 <- names(data.rr)
names(names1) <- substr(names1, 1, 5)
un.names <- unique(substr(names1, 1, 5))
unmonlist <- list()
for(i in 1 : length(un.names)) {
	unmonlist[[i]] <- names1[which(names(names1) == un.names[i])]
}
names(unmonlist) <- un.names








#######
# get medicare data
med.path <- "/dexter/disk2/rpeng/medicare/MCAPS1999-2010"
fips <- names(unmonlist)
healthdat <- list()
for(i in 1 : length(fips)) {
    dat <- try(readRDS(file.path(med.path, paste0(fips[i], ".rds"))))
    #get running temp, format data    
    healthdat[[i]] <- rest.data(dat)
    
}
names(healthdat) <- fips

########








########
# Perform SHARE and mAPCA
form1 <- "factor(agecat) + factor(dow) + ns(tmpd, df = 6) + ns(tmp3, df = 6)  + ns(dptp, df = 3) + ns(dpt3, df = 3) + ns(date, df = 8 * "
gv <- "agecat"


share <- sharehealth(data.rr, healthdata = healthdat, 
    tots = pm25, list = unmonlist, 
    formula = form1, lag = lagi, groupvar = gv)
    
iqrs <- share$summary[, "IQR"]
    
mapca <- sharehealth(data.rr, healthdata = healthdat, 
    tots = pm25, 
    list = unmonlist, method = "mapca", 
        formula = form1, iqrs = iqrs, lag = lagi, groupvar = gv)



save(share, mapca, file = paste0("east_hosp_lag", lagi, 
	"_seas", season, ".RData"))

#########











