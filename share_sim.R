#working directories
dir1 <- "/Users/jennakrall/Dropbox/SpatialFA/data"
dircode <- "/Users/jennakrall/Dropbox/SpatialFA/rcode"
home.dir <- "/Users/jennakrall/Dropbox/SpatialFA"



#load data
load(file.path(dir1, "speciation_medicare.RData"))
load("/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/cache/mnvar.RData")


#get code
source("/Users/jennakrall/Dropbox/SpatialFA/rcode/functions/mortdat_31oct13.R")
source("/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/src/absolutepca_5mar12.R")
source(file.path(dircode, "load_file_spatialfa_17may13.R"))
source(file.path(dircode, "medicare/share_sim_fn.R"))



#use Allegheny County, PA data 
#(city with large pollution, many days of data)
dat <- datall[[61]][, -c(1, 2)]









#######
#######
#######
#######
# Get info about data using PCA
#######
#######
#######
#######

stddat <- stdize(dat)
pr <- prcomp(stddat, retx = T)
nc <- length(which(pr$sd > 1))
vec <- as.matrix(varimax(pr$rot[, 1 : nc])$load[1: 24, ])


#make positive
cs <- colSums(as.matrix(vec))
if(length(which(cs < 0)) > 0) {
	vec[, which(cs < 0)] <- -vec[, which(cs < 0)]
}
for(i in 1 : ncol(vec)) {
	wh0 <- which(vec[, i] < 0)
	if(length(wh0) > 0) {
		vec[wh0, i] <- rep(0, length(wh0))
	}
}


#name sources based on vec
names <- c("traffic", "fireworks", "soil", 
	"sec sulf", "salt", "metals", "P/V")
	
cms <- rep(mnvar[1, ], 2)
sds <- rep(sqrt(mnvar[2, ]), 2)
	
	
#set up 5 subregions
keeps <- list()
#with all
keeps[[1]] <- seq(1, nc)
#minus 2
keeps[[2]] <- seq(1, nc)[-c(1, 2)]
#minus 3 
keeps[[3]] <- seq(1, nc)[-c(1, 2, 3)]
#minus 3
keeps[[4]] <- seq(1, nc)[-c(1, 2, 4)]
#minus 3
keeps[[5]] <- seq(1, nc)[-c(1, 2, 5)]












	
	
	
	
	





















#######
#######
#######
#######
# Simulate data
#######
#######
#######
#######
ns <- 10
nd <- 1000

seeds1 <- c(3474, 4866, 3451, 4672, 9165, 2165)
simout <- list()

set.seed(seeds1[1])
simout[[1]] <- multsims(ns, names, 25, 5, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds)



set.seed(seeds1[2])
simout[[2]] <- multsims(ns, names, 100, 20, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds)
	
	
set.seed(seeds1[3])
simout[[3]] <- multsims(ns, names, 5, 1, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds)	
	
	
	
set.seed(seeds1[4])
simout[[4]] <- multsims(ns, names, 25, 25, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds)	
	

#unequal
stops <- c(12, 15, 16, 20)
#12, 3, 1, 4, 5
set.seed(seeds1[5])
simout[[5]]<- multsims(ns, names, 25, 0, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds, unequal = stops)	
	


#days
days <- rep(c(200, 500, 1000, 1000, 5000), each = 5)
set.seed(seeds1[6])
simout[[6]] <- multsims(ns, names, 25, 5, ndays = nd, 
	PCs = vec, keeps = keeps, 
	cms = cms, sds = sds, days = days)	
	


library(xtable)
tabs <- t(sapply(simout, function(x) x[[2]]))
rownames(tabs) <- c("25", "100", "5", "same", "unequal subregions", "different days")
xtable(tabs)






save(dat, vec, names, keeps, simout, 
	file = file.path(home.dir, "data", "share_sim_res.RData"))






