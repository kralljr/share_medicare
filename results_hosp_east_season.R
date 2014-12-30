# Health effects results for each

setwd("~/Dropbox/SpatialFA/rcode/medicare/rdata")
plot.dir <- "~/Dropbox/SpatialFA/paper_spatialfa/figs/"


#load libraries
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(share)
library(sharesim)
library(Hmisc)

seas <- c("c", "w")
name <- "east_hosp_lag"
mapca1 <- list()
share1 <- list()
iqrs <- list()
for(j in 1 : 2) {
	mapca1[[j]] <- list()
	share1[[j]] <- list()
	for(i in 1 : 3) {

    		load(paste0(name, i - 1, "_seas", seas[j], ".RData"))
    		mapca1[[j]][[i]] <- mapca
    		share1[[j]][[i]] <- share
    	}
    	names(mapca1[[j]]) <- paste0("lag", seq(0, 2))
    	names(share1[[j]]) <- paste0("lag", seq(0, 2))
    	
    	iqrs[[j]] <- share$summary[, "IQR"]
    	
}
names(mapca1) <- paste0("season", seas)
names(share1) <- paste0("season", seas)




#MATCH SHARE AND MAPCA (lag doesn't matter)
sharem.cold <- share1[["seasonc"]][[1]]$major.sig
mapcam.cold <- mapca1[["seasonc"]][[1]]$major.sig
matches.cold <- whichCS(matchfun(list(sharem.cold), mapcam.cold, thres = 70 * pi/180)$match[[1]])

names1 <- rownames(sharem.cold)
apply(sharem.cold, 2, function(x) names1[order(abs(x), decreasing = T)[1 : 5]])
names.cold <- c("Traffic", "Sec. Sulfate","Soil", 
"Residual oil", "Salt", "As/Se/Br", "P", "Metals")

#MATCH SHARE AND MAPCA (lag doesn't matter)
sharem.warm <- share1[["seasonw"]][[1]]$major.sig
mapcam.warm <- mapca1[["seasonw"]][[1]]$major.sig
matches.warm <- whichCS(matchfun(list(sharem.warm), mapcam.warm, thres = 70 * pi/180)$match[[1]])

names1 <- rownames(sharem.warm)
apply(sharem.warm, 2, function(x) names1[order(abs(x), decreasing = T)[1 : 5]])
names.warm <- c("Sec. Sulfate", "Soil", "Fireworks", "Metals", "Salt", "P/V", "Residual oil", "As/Se/Br")


names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", 
            "Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")



sumfun <- function(x, iqrs, matches, names, type = "share") {
    #x <- x$iqrinc
    x <- x$regcoef
    lb <- x[, 1] - 1.96 * x[, 2]
    ub <- x[, 1] + 1.96 * x[, 2]
    x <- data.frame(x[, 1], lb, ub)
    
    if(type == "share") {
        #x <- data.frame(x[[1]], x[[2]])
    }else if(type != "share") {
        x <- suppressWarnings(x[matches, ])
        names <- names[matches]
        names[is.na(names)] <- ""
    }
    
    for(j in 1 : nrow(x)) {
        x[j, ] <- percinc(x[j, ], scale = iqrs[j])
    }
    
    rownames(x) <- NULL
    x <- data.frame(names, x)
    colnames(x) <- c("source", "est", "lb", "ub")
    x
}


#get IQRs, make merge better

ishare.cold <- ldply(sapply(share1[[1]], sumfun, iqrs = iqrs[[1]], 
	matches = matches.cold, names = names.cold, simplify = F), data.frame)
imapca.cold <- ldply(sapply(mapca1[[1]], sumfun, iqrs = iqrs[[1]],
	matches = matches.cold, names = names.cold, type = "mapca", simplify = F), data.frame)
colnames(ishare.cold)[1] <- "lag"
colnames(imapca.cold)[1] <- "lag"

ishare.warm <- ldply(sapply(share1[[2]], sumfun, iqrs = iqrs[[2]], 
	matches = matches.warm, names = names.warm, simplify = F), data.frame)
imapca.warm <- ldply(sapply(mapca1[[2]], sumfun, iqrs = iqrs[[2]],
	matches = matches.warm, names = names.warm, type = "mapca", simplify = F), data.frame)
colnames(ishare.warm)[1] <- "lag"
colnames(imapca.warm)[1] <- "lag"

res <- list(ishare.cold, imapca.cold)
names(res) <- c("share", "mapca")
res <- ldply(res, data.frame)
colnames(res)[1] <- "method"
res.cold <- res

res <- list(ishare.warm, imapca.warm)
names(res) <- c("share", "mapca")
res <- ldply(res, data.frame)
colnames(res)[1] <- "method"

res <- list(res.cold, res)
names(res) <- c("cold", "warm")
res <- ldply(res, data.frame)
colnames(res)[1] <- "season"









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

lb1 <- -2
ub1 <- 4

res$lag <- factor(res$lag, levels = c("lag0", "lag1", "lag2"), 
                   labels = c("Lag 0", "Lag 1", "Lag 2"))
res$method <- factor(res$method, levels = c("share", "mapca"), 
                    labels = c("SHARE", "mAPCA"))
res$season<- factor(res$season, levels = c("cold", "warm"), labels = c("Cold", "Warm"))    
ubsind <- ifelse(res$ub > ub1, 1, 0)
lbsind <- ifelse(res$lb < lb1, 1, 0)

inds <- ifelse(ubsind == 0 & lbsind == 0, 0, ifelse(ubsind == 1 & lbsind == 0, 
    1, ifelse(ubsind == 0 & lbsind == 1, 2, 3)))

ub2 <- ifelse(res$ub >ub1, Inf, res$ub)
lb2 <- ifelse(res$lb < lb1, -Inf, res$lb)

res <- data.frame(res, inds, ub2, lb2)

limits <- aes(ymax = ub, ymin= lb)

#order by value of share for lag 0



s1 <- c("Metals", "Traffic", "Residual oil", "Soil", "Salt", "Sec. Sulfate")

est <- share1[[1]][[1]]$iqrinc
rownames(est) <- capitalize(names.cold)
est <- est[rownames(est) %in% s1, , drop = F]
ord1 <- rownames(est)[order(est[, 1], decreasing = T)]

res <- res[which(res$source %in% s1), ]

res$source <- factor(res$source, levels = ord1)

res1 <- res

load("hospeast_year.RData")
res <- data.frame(rep("Year", nrow(res)), res)
colnames(res)[1] <- "season"

resall <- merge(res, res1, all = T)







########
########
########
########
########
########
########
########
### PLOT RESULTS


resall <- resall[-which(resall$source == "Metals" & resall$method == "SHARE" & resall$season == "Cold"), ]



size1 <- 18
sizep <- 1.2
pd <- position_dodge(.6)


#color
col1 <- c("indianred3", col1[2])
#bw
col1 <- c("black", "grey70", "grey40")
gplot1 <- function(lag, dat = tln3, nc = 1) {
    dat <- dat[which(dat$lag %in% paste("Lag", lag)), ]
    
    g1 <- ggplot(dat, aes(x = method, y = est, 
    		colour = season, shape = season), 
        	ylim = c(-0.01, 0.01)) + 
        geom_hline(aes(yintercept = 0), colour = "grey80", 
                   linetype = "dashed") +
        geom_pointrange(aes(ymin = lb2, 
        	ymax = ub2, colour = season), 
            width = 0.1, position = pd, size = sizep) +
        theme_bw() + 
        scale_color_manual( labels=c("Year", "Cold", "Warm"), 
        	name = "", values = col1) + 
        scale_shape( labels=c("Year", "Cold", "Warm"), 
        	name = "") + 	
        ylab(expression(atop("% increase in CVD hospitalizations per", 
           	"IQR increase in source concentration"))) +
        xlab("") +
        # ggtitle("Mortality effects by\nsource apportionment method") +
        #scale_y_continuous(limits=c(lb1, ub1))      +         
        theme_bw() +  
        theme(axis.text.y=element_text(size=size1)) +

        theme(axis.title=element_text(size=size1)) 
    theme(legend.justification=c(1,0), 
          legend.position=c(1,0)) # Position legend in bottom right
    g1 <- g1 + theme(axis.text.x=element_text(size = size1, hjust = 1,
                                              vjust = 1, angle = 45))
    g1 <- g1 + theme(strip.text.x = element_text(size = size1))
    #if(length(lag) > 1) {
        g1 <- g1 + facet_wrap(~source, scales = "free_y", ncol = nc) 
    #}
    g1 <- g1	+ theme(strip.text.y = element_text(size = size1))
    
    g1
}




setwd(plot.dir)
pdf("hosp_east_lag0_season.pdf", height = 7, width = 7)
gplot1(c(0),  dat = resall, nc = 2)
graphics.off()




# check sig
ressig <- resall[which(resall$lb > 0 | resall$ub < 0), ]
ressig[order(ressig$source), ]

# ###########	
###########	
###########	
###########	
###########	













