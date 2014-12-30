#####  File to get associations between hospitalizations and sources in east


rm(list = ls())
args <- commandArgs(TRUE)
lagi <- as.numeric(args[[1]])
print(lagi)



#load my libraries
library(share)
library(handles)
library(sharesim)
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
##########








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



seeds <- c(2746, 2723, 6158)
set.seed(seeds[lagi + 1])

share <- sharehealth(data.rr, healthdata = healthdat, 
    tots = pm25, list = unmonlist, 
    formula = form1, lag = lagi, groupvar = gv)
    
    
mapca <- sharehealth(data.rr, healthdata = healthdat, 
    tots = pm25, 
    list = unmonlist, method = "mapca", 
        formula = form1, lag = lagi, groupvar = gv)



save(share, mapca, file = paste0("east_hosp_lag", lagi, ".RData"))

#########











