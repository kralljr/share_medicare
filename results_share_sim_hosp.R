# File to get results for share simulation for hospitalizations
setwd("~/Dropbox/SpatialFA/rcode/medicare/")

library(sharesim)
library(Hmisc)
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)

load("Asimhosp_multres.RData")
outA <- out
outAmult <- outmult

mseA <- msefun(percinc1, etas, rownames(outmult[[1]][[1]][[1]]))
outAres <- gethospsim(regcoef, iqrs)[[2]]

load("Bsimhosp_multres.RData")
outB <- out
outBmult <- outmult

mseB <- msefun(percinc1, etas, rownames(outmult[[1]][[1]][[1]]))
outBres <- gethospsim(regcoef, iqrs)[[2]]











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
shareO <- data.frame(ovall1[, 1 : 2], ovall[, 1 : 2])
mAPCAO <- data.frame(ovall1[, 3 : 4], ovall[, 3 : 4])
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









