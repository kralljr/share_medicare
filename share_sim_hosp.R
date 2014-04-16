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

out1 <- reorderout(list(outB, outA), 
	nr = nrow(outA[[1]]), 	sources = rownames(outA[[1]]))


out1$Type <- factor(out1$Type, levels = c("Known", "SHARE", "mAPCA"))
out1$Source <- factor(out1$Source, levels = names)
# out1$Sim <- factor(out1$Sim, levels = c("B", "A"), )

pd <- position_dodge(.4)
size1 <- 18
sizep <- 1.2
ub1 <- 5
lb1 <- -2.5
ub2 <- ifelse(out1$ub >ub1, Inf, out1$ub)
lb2 <- ifelse(out1$lb < lb1, -Inf, out1$lb)


col1 <- brewer.pal(5, "Dark2")
g1 <- ggplot(out1, aes(x = Source, y = est, colour = Type), 
	ylim = c(-0.01, 0.01)) + 
	geom_hline(aes(yintercept = 0), colour = "grey80", 
		linetype = "dashed") +
    geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = Type), 
  	  width = 0.1, position = pd, size = sizep) +
    scale_color_manual( labels =c("Known", "SHARE", "mAPCA"), name = "",
                values = col1) + 

    ylab(expression(atop("% increase in CVD hospitalizations per", 
		 "IQR increase in source concentration"))) +
    xlab("") +
    # ggtitle("Mortality effects by\nsource apportionment method") +
    scale_y_continuous(limits=c(lb1, ub1))      +         
    theme_bw() +        
    theme(axis.text.y=element_text(size=size1)) +
  	theme(legend.text=element_text(size=size1), legend.position = "bottom")  +
  	theme(axis.title=element_text(size=size1)) 

g1 <- g1 +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size = size1))
g1 <- g1	+ theme(strip.text.y = element_text(size = size1)) +
	facet_wrap(~Sim, ncol = 2) +theme(strip.text.x = element_text(size = size1))
g1 <- g1 + theme(legend.justification=c(0,0), legend.position=c(0,0))



# pdf(file.path(dirplot, "hosp_sim.pdf"), height = 7, width = 10)
g1
# graphics.off()

