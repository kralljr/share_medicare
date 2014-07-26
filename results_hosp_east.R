# Health effects results for each

setwd("~/Dropbox/SpatialFA/rcode/medicare/rdata")
plot.dir <- "~/Dropbox/SpatialFA/paper_spatialfa/figs/"


#load libraries
library(ggplot2)
library(plyr)
library(RColorBrewer)

name <- "east_hosp_lag"
mapca1 <- list()
share1 <- list()
for(i in 1 : 3) {
    load(paste(name, i, ".RData"))
    mapca1[[i]] <- mapca
    share1[[i]] <- share
}


names(mapca1) <- paste0("lag", seq(0, 2))
names(share1) <- paste0("lag", seq(0, 2))



#MATCH SHARE AND MAPCA (lag doesn't matter)
sharem <- share1[[1]]$major.sig
mapcam <- mapca1[[1]]$major.sig
matches <- whichCS(matchfun(list(sharem), mapcam, thres = 70 * pi/180)$match[[1]])




names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", 
            "Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")

sumfun <- function(x, type = "share") {
    x <- x$iqrinc
    if(type == "share") {
        x <- data.frame(x[[1]], x[[2]])
    }else if(type != "share") {
        x <- suppressWarnings(x[matches, ])
    }
    x <- data.frame(names1, x)
    colnames(x) <- c("source", "est", "se")
    x
}
ishare <- ldply(sapply(share1, sumfun, simplify = F), data.frame)
imapca <- ldply(sapply(mapca1, sumfun, type = "mapca", simplify = F), data.frame)
colnames(ishare) <- c("lag","source",  "est", "se")
colnames(imapca) <- c("lag", "source", "est", "se")

res <- list(ishare, imapca)
names(res) <- c("share", "mapca")
res <- ldply(res, data.frame)
colnames(res)[1] <- "method"


lb <- res$est - 1.96 * res$se
ub <- res$est + 1.96 * res$se

res <- data.frame(res, lb, ub)









########
########
########
########
########
########
########
########
### Get result



res <- res[complete.cases(res), ]

lb1 <- -1
ub1 <- 2

res$lag <- factor(res$lag, levels = c("lag0", "lag1", "lag2"), 
                   labels = c("Lag 0", "Lag 1", "Lag 2"))
res$method <- factor(res$method, levels = c("share", "mapca"), 
                    labels = c("SHARE", "mAPCA"))
ubsind <- ifelse(res$ub > ub1, 1, 0)
lbsind <- ifelse(res$lb < lb1, 1, 0)

inds <- ifelse(ubsind == 0 & lbsind == 0, 0, ifelse(ubsind == 1 & lbsind == 0, 
    1, ifelse(ubsind == 0 & lbsind == 1, 2, 3)))

ub2 <- ifelse(res$ub >ub1, Inf, res$ub)
lb2 <- ifelse(res$lb < lb1, -Inf, res$lb)

res <- data.frame(res, inds, ub2, lb2)

limits <- aes(ymax = ub, ymin= lb)


# s1 <- c("Metals", "Traffic", "Residual oil", "Soil", "Salt", "Sec. Sulfate")
# res <- res[which(res$source %in% s1), ]
# res$source <- factor(res$source, levels = s1)












########
########
########
########
########
########
########
########
### PLOT RESULTS




size1 <- 18
sizep <- 1.2
pd <- position_dodge(.4)


col1 <- brewer.pal(5, "Blues")[3:5]
col1 <- brewer.pal(5, "Dark2")[c(1, 4)]
col1 <- brewer.pal(8, "Set1")[1:2]
col1 <- c("indianred3", col1[2])
gplot1 <- function(lag, dat = tln3) {
    dat <- dat[which(dat$lag %in% paste("Lag", lag)), ]
    
    g1 <- ggplot(dat, aes(x = source, y = est, colour = method), 
                 ylim = c(-0.01, 0.01)) + 
        geom_hline(aes(yintercept = 0), colour = "grey80", 
                   linetype = "dashed") +
        geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = method), 
                        width = 0.1, position = pd, size = sizep) +
        scale_color_manual( labels=c("SHARE", "mAPCA"), name = "",
                            values = col1) + 	
        ylab(expression(atop("% increase in CVD hospitalizations per", 
                             "IQR increase in source concentration"))) +
        xlab("") +
        # ggtitle("Mortality effects by\nsource apportionment method") +
        scale_y_continuous(limits=c(lb1, ub1))      +         
        theme_bw() +        
        theme(axis.text.y=element_text(size=size1)) +

        theme(axis.title=element_text(size=size1)) 
    theme(legend.justification=c(1,0), 
          legend.position=c(1,0)) # Position legend in bottom right
    g1 <- g1 + theme(axis.text.x=element_text(size = size1, hjust = 1,
                                              vjust = 1, angle = 45))
    g1 <- g1 + theme(strip.text.x = element_text(size = size1))
    if(length(lag) > 1) {
        g1 <- g1 + facet_wrap(~lag, ncol = 1) 
    }
    g1 <- g1	+ theme(strip.text.y = element_text(size = size1))
    
    g1
}




setwd(plot.dir)
# pdf("hosp_east_lag012.pdf", height = 10, width = 7)
gplot1(c(0,1, 2), dat = res)
# graphics.off()


# ###########	
###########	
###########	
###########	
###########	













