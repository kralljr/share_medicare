### File to simulate mortality effects


# set working directories
dircode <- ("/Users/jennakrall/Dropbox/SpatialFA/rcode")
dirplot <- ("/Users/jennakrall/Dropbox/SpatialFA/plots")
medicare <- file.path(dircode, "medicare")


load("/Users/jennakrall/Dropbox/SpatialFA/data/share_sim_res.RData")
load("/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/cache/mnvar.RData")





#load libraries
library(tlnise)
library(ggplot2)
library(lme4)
library(Hmisc)

#load code
source(file.path(dircode, "functions/mort_fn_east.R"))
source(file.path(dircode, "load_file_spatialfa_17may13.R"))
source(file.path(medicare, "share_sim_fn.R"))








#set up info
names <- names
nmons <- 25
reps <- 5
ndays <- 500
PCs <- vec
keeps <- keeps
cms <- rep(mnvar[1, ], 2)
sds <- rep(sqrt(mnvar[2, ]), 2)
nreps <- 10
etas <- c(3, 1, 0.75, 0.5, 1, 1, 2)




set.seed(398)
outA <- outerSIMhosp(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, etas)


reps <- nmons
set.seed(9763)
outB <- outerSIMhosp(names, nmons, reps, ndays, PCs, keeps, 
		cms, sds, etas)

# outA <- outB

out2 <- reorderout(list(outB, outA), 
	nr = nrow(outA[[1]]), 	sources = rownames(outA[[1]]))
out1 <- out2

# out1$Type <- factor(out1$Type, levels = c("Known", "SHARE", "mAPCA"))
out1$Source <- factor(out1$Source, levels = names)
# out1$Sim <- factor(out1$Sim, levels = c("B", "A"), )



# out1 <- out1[-which(out1$Type == "SHARE1"), ]
out1$Type <- factor(out1$Type, levels = c("Known", "SHARE", "mAPCA"))



ub1 <- 6
lb1 <- -2.5
ub2 <- ifelse(out1$ub >ub1, Inf, out1$ub)
lb2 <- ifelse(out1$lb < lb1, -Inf, out1$lb)


out2 <- out1
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



pdf(file.path(dirplot, "hosp_sim.pdf"), height = 6, width = 15)
g1
graphics.off()








#####
# try table

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
colnames(out) <- c("Source", "Eta", "A: SHARE est", 
	"A: SHARE ci", "A: mAPCA est", "A: mAPCA ci")



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
colnames(outB) <- c("Source", "Eta", "B: SHARE est", 
	"B: SHARE ci", "B: mAPCA est", "B: mAPCA ci")


outall <- data.frame(out, outB[, -c(1, 2)])

library(xtable)
print(xtable(out), include.rownames = F)
print(xtable(outB), include.rownames = F)





######
# try usual plot

getdat <- function(type, sim) {
	out <- out1[which(out1$Type == type & out1$Sim == sim), ]
	out
}


getpoints <- function(try, cols, pchs, shift = 0) {
	seqs <-(seq(1, nrow(try))) + shift
	points(seqs, try$est, col = cols,
		pch = pchs)
	segments(x0 = seqs, y0 = try$lb, 
		y1 = try$ub, col = cols)
	# points(try$est, seqs, col = cols,
		# pch = pchs)
	# segments(y0 = seqs, x0 = try$lb, 
		# x1 = try$ub, col = cols)	
}


segtruth <- function(try, etas, cols) {
	seqs <- (seq(1, nrow(try)))
	# segments(x0 = seqs - 0.5, x1 = seqs + 0.5, 
		# y0 = etas, col = cols, lwd = 2)
	# segments(y0 = seqs - 0.5, y1 = seqs + 0.5, 
		# x0 = etas, col = cols, lwd = 2)	
	# points(etas, seqs + 0.25, col = cols, pch = "+")
	points(seqs - 0.25, etas, col = cols, pch = "+")	
}


plotall <- function(sim, etas, ylims = c(-2, 5), side1 = 2, side2 = 1) {
	
		
	cols <- brewer.pal(4, "Dark2")
	col1 <- brewer.pal(8, "Set1")[1:2]
	cols <- c("indianred3", col1[2])
	pchs <- c(16, 15)
	pchs <- c(16, 16)
	
	try1 <- getdat("SHARE", sim)
	try1 <- try1[-which(try1$Source == "P/V"), ]
	# try1$Source <- factor(try1$Source, levels = c("traffic", 
		# "fireworks", "metals", "salt", "soil", "sec sulf"), 
		# labels = c("Traffic", "Fireworks", "Metals", "Salt",
		# "Soil", "Secondary Sulfate"))
	ords <- order(try1$Source)
	try1 <- try1[ords, ]	
	etas <- etas[-7]
	plot(1, 1, xlim = c(0.5, nrow(try1) + 0.5), ylim = ylims, 
		axes = F, type = "n", xlab = "",
		ylab = "", main = "")
	# plot(1, 1, ylim = c(0.5, nrow(try1) + 0.5), xlim = ylims, 
		# axes = F, type = "n", xlab = "",
		# ylab = "", main = "")	
	mtext(paste0(sim, "."), side = 3, line = -2)
	if(sim == "A") {
		mtext(expression(atop("% increase in CVD hospitalizations per", 	
		"IQR increase in source concentration")), side = side1, line = 2,
		outer = T)
		# mtext("% increase in CVD hospitalizations per IQR increase in source concentration", side = side1, line = 2,
		# outer = T)

		# axis(side2, at = seq(1, nrow(try1)),
		 # labels = sour1, las = 1)
		axis(side1)
	}
	
			sour1 <- as.character((try1$Source))
		sour1 <- capitalize(sour1)
		sour1[which(sour1 == "Sec sulf")] <- "Secondary sulfate"
		
	box()
	axis(side2, at = seq(1, nrow(try1)),
		 labels = sour1, las = 2)
	# axis(side1)
	abline(h = 0, col = "grey20", lty = 2)	 
	# abline(v = 0, col = "grey20", lty = 2)	 
	
	segtruth(try1, etas, 1)	 
	getpoints(try1, cols[1], pchs[1])
	
	try2 <- getdat("mAPCA", sim)
	try2 <- try2[-which(try2$Source == "P/V"), ]
	# try2$Source <- factor(try2$Source, levels = c("traffic", 
		# "fireworks", "metals", "salt", "soil", "sec sulf"), 
		# labels = c("Traffic", "Fireworks", "Metals", "Salt",
		# "Soil", "Secondary Sulfate"))
	try2 <- try2[ords, ]
	getpoints(try2, cols[2], pchs[2], shift = 0.25)
	
	if(sim == "A") {
		legend("bottomleft", legend = c("Truth", "SHARE", "mAPCA"),
			col = c(1, cols), pch = c(3, 16, 16))
	}
	}



pdf("sim_hosp_baser_exty.pdf", height = 7, width = 13)
par(mfrow = c(1, 2), oma = c(0, 5.5, 0, 2), mar = c(8, 0, 3, 0))
# par(mfrow = c(1, 2), oma = c(4, 0, 0, 2), mar = c(3, 0, 2, 0))
plotall("A", etas, ylims = c(-6, 7))
# par(mar = c(7, 3, 3, 4.5))
plotall("B", etas, ylims = c(-6, 7))
#, side1 = 1, side2 = 2
graphics.off()



pdf("sim_hosp_baser.pdf")
par(mfrow = c(1, 2), oma = c(0, 5.5, 0, 2), mar = c(7, 0, 3, 0))
plotall("A", etas, ylims = c(-2, 5))
# par(mar = c(7, 3, 3, 4.5))
plotall("B", etas, ylims = c(-2, 5))
graphics.off()