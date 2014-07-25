#####  File to get associations between hospitalizations and sources in east

# set working directories
dircode <- ("~/Dropbox/SpatialFA/rcode")
home.dir <- "~/Dropbox/SpatialFA"



#get data
load("/Users/jennakrall/Dropbox/PM25cons_mort/FinalR_815/Healtheffectsdata_17aug11.RData")
# load("/Users/jennakrall/Dropbox/PM25cons_mort/FinalR_815/Prednaivedata_17aug11.RData")
# spec <- readRDS("~/Dropbox/SpatialFA/data/speciation_monitors.rds")
counties <- readRDS("/Users/jennakrall/Dropbox/PM25cons_mort/cities/nmmaps_counties.rds")

load(file.path(home.dir, "data/pvdeast_medicare.RData"))
load(file.path(home.dir, "data/mort_sources_medicare.RData"))
load(file.path(home.dir, "data/sourceconc_medicare.RData"))

#datall, monsKEEP, unmons
load(file.path(home.dir, "data/speciation_medicare.RData"))

load("/Users/jennakrall/Dropbox/SpatialFA/data/link_east_county.RData")


#load other libraries
library(splancs)
library(gpclib)
library(splines)
library(tsModel)
library(tlnise)
library(chron)
library(timeDate)
library(ggplot2)


#load my libraries
library(share)
library(handles)
library(sharesim)








########

########
# Perform SHARE and mAPCA
sharehealth <- function(consdata, healthdata = NULL, list = NULL, 
                        formula = NULL, iqrs = NULL,
                        lag = NULL, groupvar = NULL, print = F, 
                        cut = 1, thres = pi/4,
                        tots = NULL, method = "SHARE") {
    
share <- list()
mapca <- list()
for(i in 1 : 3) {
    share[[i]] <- sharehealth(data.rr, tots = pm25, list = unmonlist, 
        formula = form1, lag = i, groupvar = gv)
    iqrs <- share$summary[, "IQR"]
    mapca[[i]] <- sharehealth(data.rr, tots = pm25, list = unmonlist, method = "mapca", 
        formula = form1, iqrs = iqrs, lag = i, groupvar = gv)
}
names(share) <- paste0("lag", seq(0, 2))
names(mapca) <- paste0("lag", seq(0, 2))













