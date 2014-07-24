# File to run simulation study for SHARE
# No health, just source identification


########
#load packages
library(devtools)
install_github("handles", "kralljr")
library(handles)
install_github("share", "kralljr")
library(share)
install_github("share_medicare", "kralljr", subdir = "sharesim")
library(sharesim)

########



########
#load data
data(data_share_sim)












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






