# File to get results for share simulation for hospitalizations
setwd("~/Dropbox/SpatialFA/rcode/medicare/rdata")

library(sharesim)
library(share)
library(Hmisc)
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)



for(i in 1 : 3) {

load(paste0("Asimhosp_multres_sderr", i, ".RData"))
outAmult <- outmult

mseA <- mse
outAres <- outrest

load(paste0("Bsimhosp_multres_sderr", i, ".RData"))
outBmult <- outmult

mseB <- mse
outBres <- outrest










### 
# Table of average across 100 simulated datasets

#format confidence interval and make into dataframe
rn <- rownames(outBmult[[1]][[1]][[1]])
ovall <- sapply(outAres, function(x) {
    x <- round(x, 2)
    x <- cbind(x[, 1], paste0( "(", x[, 2],", ", x[, 3], ")"))
    x
    }, simplify = F)[2 :3 ]
ovall <- data.frame(ovall[[1]], ovall[[2]])
rownames(ovall) <- rn
ovallA <- ovall


#same for sim B
ovall <- sapply(outBres, function(x) {
    x <- round(x, 2)
    x <- cbind(x[, 1], paste0( "(", x[, 2],", ", x[, 3], ")"))
    x
    }, simplify = F)[2 :3 ]
ovall <- data.frame(ovall[[1]], ovall[[2]])
rownames(ovall) <- rn
ovallB <- ovall


#get MSES
shareO <- data.frame(ovallA[, 1 : 2], ovallB[, 1 : 2])
mAPCAO <- data.frame(ovallA[, 3 : 4], ovallB[, 3 : 4])
colnames(shareO) <- rep(c("A", "B"), each = 2)
colnames(mAPCAO) <- rep(c("A", "B"), each = 2)
ovallO <- rbind(shareO, mAPCAO)


#format results
mseAPCA <- round(c(mseA[, "APCA"], mseB[, "APCA"]), 2)
names(mseAPCA) <- NULL
msemAPCA <- round(c(mseA[, "mAPCA"], mseB[, "mAPCA"]), 2)
names(msemAPCA) <- NULL
rownames(ovall) <- NULL
ovallO <- data.frame(rep(c("SHARE", "mAPCA"), each = nrow(ovallO)/2),
                     rep(rn, 2),
                     rep(etas, 2), 
                     ovallO[, 1 : 2], mseAPCA,
                     ovallO[, 3 : 4], msemAPCA)
colnames(ovallO) <- c("Method", "Source", "eta", "estA","ciA", "mseA", "estB","ciB", "mseB")
ovallO <- ovallO[-which(ovallO$Source == "P/V"), ]


ovallO$Source <- capitalize(as.character(ovallO$Source))

library(xtable)
print(xtable(ovallO), include.rownames = F)

}

