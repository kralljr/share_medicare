# File to get results for share simulation for hospitalizations
setwd("~/Dropbox/SpatialFA/rcode/medicare/rdata")

library(sharesim)
library(share)
library(Hmisc)
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)

load("Asimhosp_multres.RData")
outA <- out
outAmult <- outmult

mseA <- mse
outAres <- outrest

load("Bsimhosp_multres.RData")
outB <- out
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




#Make plot for presentations
outAres <- outAres[2 : 3]
names(outAres) <- c("SHARE", "mAPCA")
outBres <- outBres[2 : 3]
names(outBres) <- c("SHARE", "mAPCA")
rownames(outAres[[1]]) <- capitalize(rn)
rownames(outAres[[2]]) <- capitalize(rn)
rownames(outBres[[1]]) <- capitalize(rn)
rownames(outBres[[2]]) <- capitalize(rn)

library(reshape2)
mA <- melt(outAres)
mA <- data.frame(mA, rep("A", nrow(mA)))
colnames(mA) <- c("source", "variable", "value",
	"method", "sim")
cA <- reshape(mA, 
	 idvar = c("source", "sim", "method"), 
	 timevar = "variable",
	 direction = "wide")
colnames(cA)[4 : 6] <- c("est", "lb", "ub")

mB <- melt(outBres)
mB <- data.frame(mB, rep("B", nrow(mB)))
colnames(mB) <- c("source", "variable", "value",
	"method", "sim")
cB <- reshape(mB, 
	 idvar = c("source", "sim", "method"), 
	 timevar = "variable",
	 direction = "wide")
colnames(cB)[4 : 6] <- c("est", "lb", "ub")


merged <- merge(cA, cB, all = T)







etas1 <- data.frame(matrix(rep(etas, 2)), matrix(rep(capitalize(names), 2)))
etas1 <- data.frame(etas1, etas1[, 1], etas1[, 1])
etas1 <- data.frame(etas1, rep(c("A", "B"), each = length(etas)))
etas1 <- data.frame(etas1, rep("Truth", nrow(etas1)))
colnames(etas1) <- c("est", "source", "lb", "ub", "sim", "method")

merged1 <- merge(merged, etas1, all = T)


merged1 <- merged1[-which(merged1$source == "P/V"), ]
sour1 <- capitalize(as.character(unique(merged1$source)))
sour1[which(sour1 == "Sec sulf")] <- "Sec. sulfate"
merged1$source <- factor(merged1$source, levels = unique(merged1$source)[-7], 
	labels = sour1)
 




ub1 <- 6
lb1 <- -6
ub2 <- ifelse(merged1$ub >ub1, Inf, merged1$ub)
lb2 <- ifelse(merged1$lb < lb1, -Inf, merged1$lb)



merged1 <- data.frame(merged1, ub2, lb2)

merged1$method <- factor(merged1$method, levels = c("Truth", "SHARE", "mAPCA"))


library(ggplot2)
library(RColorBrewer)



pd <- position_dodge(.4)
size1 <- 22
sizep <- 1.7
col1 <- brewer.pal(5, "Dark2")

col1 <- brewer.pal(8, "Set1")[1:2]
cols <- c("indianred3", col1[2])
col1 <- c(1, cols)
g1 <- ggplot(merged1, aes(x = source, y = est, colour = method), 
	ylim = c(-0.01, 0.01)) + 
	geom_hline(aes(yintercept = 0), colour = "grey80", 
		linetype = "dashed") +
    geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = method), 
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
	facet_wrap(~sim, ncol = 2) +theme(strip.text.x = element_text(size = size1))
g1 <- g1 + theme(legend.justification=c(0,0), legend.position="right")



pdf(file.path("~/Dropbox/SpatialFA/plots", "hosp_sim.pdf"), height = 6, width = 15)
g1
graphics.off()


