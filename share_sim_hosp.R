### File to simulate mortality effects
rm(list = ls())
args <- commandArgs(TRUE)
sim <- as.numeric(args[[1]])
print(sim)






#load others
library(share)
library(handles)
library(sharesim)






#load data
data(keeps, cms, sds, vec, names)

# new combo of sources
load("~/SHARE/sharesim/data_share_sim_revisedcombo.RData")


#set up info
nmons <- 25
ndays <- 500
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)
seeds <- c(9763, 398)
reps1 <- c(25, 5)


###############
reps <- reps1[sim]
outmult <- list()
set.seed(seeds[sim])
for(i in 1 : 100) {
    print(i)
    outmult[[i]] <- outerSIMhosp(names, nmons, reps, ndays, vec, keeps, 
                                  cms, sds, etas)
}




regcoef <- sapply(outmult, function(x) x[[1]], simplify = F)
percinc1 <- sapply(outmult, function(x) x[[2]], simplify = F)
iqrs <- apply(sapply(outmult, function(x) x[[3]]), 1, median)
outrest <- gethospsim(regcoef, iqrs)[[2]]
mse <- msefun(percinc1, etas, rownames(outmult[[1]][[1]][[1]]))
#first iteration, percinc
out <- outmult[[1]][[2]]








# save output
sims <- c("A", "B")
save(out, outmult, mse, outrest, iqrs, percinc1, regcoef,
     file = paste0(sims[sim], "simhosp_multres.RData"))


