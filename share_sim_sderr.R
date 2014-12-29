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



########
# try new combination
load("data_share_sim_revisedcombo.RData")




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
seeds1 <- c(4658, 4260, 8656, 8745, 5075)
mon <- 25
reg <- 5

#set params for simulation
keeps1 <- keeps
uneqs <- NULL
day1 <- NULL
    
    
set.seed(seeds1[seed])

sderrs <- c(0.001, 0.1, 0.5, 1, 10)
sderr1 <- sderrs[seed]

simout <- multsims(nsims = ns, names = names, 
    nmons = mon, reps = reg, ndays = nd, 
    PCs = vec, keeps = keeps1, 
    cms = cms, sds = sds, 
    unequal = uneqs, days = day1, cut = 1, thres = pi/4, prnt = T,
    sderr = sderr1)	


save(simout, file = paste0("simout_sderr_", seed, ".RData"))






