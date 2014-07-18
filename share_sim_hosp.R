### File to simulate mortality effects

######
# set working directories
basedir <- "/Users/jennakrall/Dropbox"
basedir <- "C:/Users/jrkrall/Dropbox"


home.dir <- file.path(basedir, "SpatialFA")
dir1 <- file.path(home.dir, "data")
dircode <- file.path(home.dir, "rcode")
dirplot <- file.path(home.dir, "plots")
medicare <- file.path(dircode, "medicare")



#load data
#datall from speciation data
load(file.path(dir1, "speciation_medicare.RData"))
#vec, keeps from share_sim.R
load(file.path(basedir, "SpatialFA/data/share_sim_res.RData"))
#data(mnvar)



load(file.path(basedir, "MDL_sourceapp/MDL_project_oct2012/cache/mnvar.RData"))




#load libraries
library(tlnise)
library(ggplot2)
library(nlme)
library(Hmisc)
library(clue)
library(devtools)


#source code
source(file.path(medicare, "share_sim_fn.R"))



#Install packages and relevant functions
#library(devtools)
#install_github("handles", "kralljr")
#library(handles)
#install_github("share", "kralljr")
#library(share)

#for now, source file
source(file.path(dircode, "loadfiles_handles_share.R"))




#set up info
names <- names
nmons <- 25
ndays <- 500
PCs <- vec
keeps <- keeps
cms <- rep(mnvar[1, ], 2)
sds <- rep(sqrt(mnvar[2, ]), 2)
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)




###############
# Case A: sources same across monitors
reps <- nmons
outAmult <- list()
set.seed(9763)
for(i in 1 : 100) {
    print(c("A", i))
    outAmult[[i]] <- outerSIMhosp(names, nmons, reps, ndays, PCs, keeps, 
                                  cms, sds, etas)
}
regcoef <- sapply(outAmult, function(x) x[[1]], simplify = F)
percinc1 <- sapply(outAmult, function(x) x[[2]], simplify = F)
iqrs <- apply(sapply(outAmult, function(x) x[[3]]), 1, median)
outAres <- gethospsim(regcoef, iqrs)[[2]]
mseA <- msefun(percinc1, etas, rownames(outAmult[[1]][[1]][[1]]))
outA <- outAmult[[1]][[2]]



###############
# Case B: sources differ across monitors
reps <- 5
outBmult <- list()
set.seed(398)
for(i in 1 : 100) {
    print(c("B", i))
    outBmult[[i]] <- outerSIMhosp(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, etas)
}


regcoef <- sapply(outBmult, function(x) x[[1]], simplify = F)
percinc1 <- sapply(outBmult, function(x) x[[2]], simplify = F)
iqrs <- sapply(outBmult, function(x) x[[3]])
outBres <- gethospsim(regcoef, iqrs)[[2]]
mseB <- msefun(percinc1, etas, rownames(outBmult[[1]][[1]][[1]]))
outB <- outBmult[[1]][[2]]







# save output
save(outAmult, outBmult, file = file.path(dircode, "simhosp_multres.RData"))



#get output as one matrix
out2 <- suppressWarnings(reorderout(list(outA, outB), 
	nr = nrow(outA[[1]]), 	sources = rownames(outAmult[[1]][[1]][[1]])))
out2$Source <- factor(out1$Source, levels = names)
out2$Type <- factor(out1$Type, levels = c("Known", "SHARE", "mAPCA"))

names(etas) <- names
sim1 <- matrix(rep(etas, each = 3), byrow = T, ncol = 3)
sim1 <- data.frame(sim1[, 1], rep(0, 7), sim1[, -1], rep("Truth", 7), 
	names)
simA <- data.frame(sim1, rep("A", 7))
colnames(simA) <- colnames(out1)
out2 <- rbind(out1, simA)
simB <- data.frame(sim1, rep("B", 7))
colnames(simB) <- colnames(out1)
out2 <- rbind(out2, simB)

out2 <- out2[-which(out2$Source == "P/V"), ]
sour1 <- capitalize(as.character(unique(out2$Source)))
sour1[which(sour1 == "Sec sulf")] <- "Sec. sulfate"
out2$Source <- factor(out2$Source, levels = names[-7], 
	labels = sour1)

out2 <- out2[-which(out2$Type == "Known"), ]
out2$Type <- factor(out2$Type, levels = c("Truth", "SHARE", "mAPCA"))

out2 <- out2[order(out2[, 7], out2[, 5]), ]









#####
#  table of results for one simulated dataset

out2A <- out2[which(out2$Sim == "A"), 
	c(1, 3, 4, 5, 6)]
out2A <- reshape(out2A, direction = "wide", 
	idvar = "Source", 
	timevar = "Type")
out2A <- out2A[, -c(3, 4)]
out2A[, -1] <- round(out2A[,-1], 2)
out <- data.frame(out2A[, c(1 : 3)], 
	paste0("(", out2A[, 4], ", ", out2A[, 5], ")"), 
	out2A[, 6], 
	paste0("(", out2A[, 7], ", ", out2A[, 8], ")"))

colnames(out) <- c("Source", "Eta","SHARE est", 
	"SHARE ci", "mAPCA est", "mAPCA ci")


out2B <- out2[which(out2$Sim == "B"), 
	c(1, 3, 4, 5, 6)]
out2B <- reshape(out2B, direction = "wide", 
	idvar = "Source", 
	timevar = "Type")
out2B <- out2B[, -c(3, 4)]
out2B[, -1] <- round(out2B[,-1], 2)
outB <- data.frame(out2B[, c(1 : 3)], 
	paste0("(", out2B[, 4], ", ", out2B[, 5], ")"), 
	out2B[, 6], 
	paste0("(", out2B[, 7], ", ", out2B[, 8], ")"))
colnames(outB) <- c("Source", "Eta", "SHARE est", 
	"SHARE ci", "mAPCA est", "mAPCA ci")
# Replace truth with simulated value

shareS <- data.frame(out[, 1 : 4], outB[, 3 : 4])
mAPCAS <-  data.frame(out[, c(1, 2, 5, 6)], outB[, 5 : 6])
colnames(shareS) <- c("Source", "eta", "A.est", "A.ci", "B.est", "B.ci")
colnames(mAPCAS) <- colnames(shareS)
outS <- rbind(shareS, mAPCAS)
outS <- data.frame(rep(c("SHARE", "mAPCA"), each = nrow(shareS)), outS)
colnames(outS)[1] <- "Method"









### 
# Table of average across 100 simulated datasets

ovall <- sapply(outBres, function(x) {
    x <- round(x, 2)
    paste0(x[, 1], " (", x[, 2], ")")})[, -1]
rownames(ovall) <- rownames(outBmult[[1]][[1]][[1]])
ovall1 <- ovall
ovall <- sapply(outAres, function(x) {
    x <- round(x, 2)
    paste0(x[, 1], " (", x[, 2], ")")})[, -1]
rownames(ovall) <- rownames(outAmult[[1]][[1]][[1]])

shareO <- data.frame(ovall1[, 1], ovall[, 1])
mAPCAO <- data.frame(ovall1[, 2], ovall[, 2])
colnames(shareO) <- c("A", "B")
colnames(mAPCAO) <- c("A", "B")
ovallO <- rbind(shareO, mAPCAO)
mseAPCA <- round(c(mseB[, "APCA"], mseA[, "APCA"]), 2)
names(mseAPCA) <- NULL
msemAPCA <- round(c(mseB[, "mAPCA"], mseA[, "mAPCA"]), 2)
names(msemAPCA) <- NULL
rownames(ovall) <- NULL
ovallO <- data.frame(rep(c("SHARE", "mAPCA"), each = nrow(ovallO)/2),
    rep(rownames(shareO), 2),
    rep(etas, 2), ovall[, 1], mseAPCA,
    ovallO[, 2], msemAPCA)
colnames(ovallO) <- c("Method", "Source", "eta", "estA", "mseA", "estB", "mseB")



library(xtable)
print(xtable(ovallO), include.rownames = F)
print(xtable(outS), include.rownames = F)























######
# plot, old

library(RColorBrewer)
ub1 <- 7
lb1 <- -6
ub2 <- ifelse(out2$ub >ub1, Inf, out2$ub)
lb2 <- ifelse(out2$lb < lb1, -Inf, out2$lb)




pd <- position_dodge(.4)
size1 <- 22
sizep <- 1.7


col1 <- brewer.pal(5, "Dark2")

col1 <- brewer.pal(8, "Set1")[1:2]
cols <- c("indianred3", col1[2])
col1 <- c(1, cols)
g1 <- ggplot(out2, aes(x = Source, y = est, colour = Type), 
             ylim = c(-0.01, 0.01)) + 
    geom_hline(aes(yintercept = 0), colour = "grey80", 
               linetype = "dashed") +
    geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = Type), 
                    width = 0.1, position = pd, size = sizep) +
    scale_color_manual(name = "",
                       values = col1) + 
    # labels =c("Known", "SHARE", "mAPCA")
    
    ylab(expression(atop("% inc. in hospitalizations per", 
                         "IQR increase in source"))) +
    xlab("") +
    # ggtitle("Mortality effects by\nsource apportionment method") +
    scale_y_continuous(limits=c(lb1, ub1))      +         
    theme_bw() +        
    theme(axis.text.y=element_text(size=size1)) +
    theme(legend.text=element_text(size=size1), legend.position = "bottom")  +
    theme(axis.title=element_text(size=size1)) 

g1 <- g1 +theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1, size = size1))
g1 <- g1	+ theme(strip.text.y = element_text(size = size1)) +
    facet_wrap(~Sim, ncol = 2) +theme(strip.text.x = element_text(size = size1))
g1 <- g1 + theme(legend.justification=c(0,0), legend.position="right")



#pdf(file.path(dirplot, "hosp_sim.pdf"), height = 6, width = 15)
g1
#graphics.off()