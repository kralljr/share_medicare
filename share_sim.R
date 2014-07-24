# File to run simulation study for SHARE
# No health, just source identification

rm(list = ls())
args <- commandArgs(TRUE)
seed <- as.numeric(args[[1]])
print(seed)

########
#load packages (make sure recent versions are installed)
library(devtools)
library(handles)
library(share)
library(sharesim)

########






#######
#######
#######
#######
# Simulate data
#######
#######
#######
#######

#universal values
ns <- 100
nd <- 1000
seeds1 <- c(3474, 4866, 3451, 4672, 9165, 2165)
mons <- c(25, 100, 5, 25, 25, 25)
regs <- c(5, 20, 1, 25, 0, 5)

#set params for simulation
keeps1 <- keeps
uneqs = NULL
day1 <- NULL
    
if(seed == 5) {
    uneqs <- c(12, 15, 16, 20)
}else if(seed == 6) {
    day1 <- rep(c(200, 500, 1000, 1000, 5000), each = 5)
}

    
    
mon <- mons[seed]
reg <- regs[seed]
set.seed(seeds1[seed])



simout <- multsims(nsims = ns, names = names, 
    nmons = mon, reps = reg, ndays = nd, 
    PCs = vec, keeps = keeps1, 
    cms = cms, sds = sds, 
    unequal = uneqs, days = day1, cut = 1, thres = pi/4)	


save(simout, file = paste0("simout", seed, ".RData"))






