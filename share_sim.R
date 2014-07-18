# File to run simulation study for SHARE
# 


########
#set  directories
basedir <- "/Users/jennakrall/Dropbox"
basedir <- "C:/Users/jrkrall/Dropbox"
home.dir <- file.path(basedir, "SpatialFA")
dir1 <- file.path(home.dir, "data")
dircode <- file.path(home.dir, "rcode")


#github
gh <- "https://github.com/kralljr/share_medicare/blob/master/"
source(file.path(gh, "share_sim_fn.R"))

########
#load packages
#library(devtools)
#install_github("handles", "kralljr")
#library(handles)
#install_github("share", "kralljr")
#library(share)

#for now, source file
source(file.path(dircode, "loadfiles_handles_share.R"))
########



########
#load data
load(file.path(dir1, "data_share_sim.RData"))



########
# load functions
source(file.path(dircode, "medicare", "share_sim_fn.R"))











#######
#######
#######
#######
# Simulate data
#######
#######
#######
#######
ns <- 5
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






save(simout, file = file.path(cwd, "share_sim_res.RData"))






