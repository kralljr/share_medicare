

########
########
########
########
########
########
########
########
### PLOT RESULTS
names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", 
            "Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")





########
########
########
########
########
########
########
########
### PLOT RESULTS





### first, get data
#check which sources, name
round(pvd.east[[1]][[1]], 2)
names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", 
            "Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")


Sources <- Sources[-which(is.infinite(Sources))]
keeps <- seq(1, length(Sources))
bounds <- c(-10, 10)
namesTLN <- c("APCA", "mAPCAsw")
pchs <- c(16, 17, 15)
cols <- c(1, "grey60", "grey20")


skips1 <- c(0, 0.1, 0.2)
start1 <- 1
tlnout <- matrix(nrow = 1, ncol = 3)
tlnoutc <- matrix(nrow = 1, ncol = 3)
#plot individual mortalities
for(lagi in 1 : 3) {
    
    
    xstart <- seq(1, length(keeps)) + skips1[lagi]
    
    skips <- c(-.2, 0.2)
    is1 <- length(namesTLN)
    for(i in 1 : is1) {
        print(c(lagi, i))
        xs <- xstart + skips[i]
        
        tlnest <- get(paste0("tln", namesTLN[i]))
        
        tlns <- sapply(tlnest[[1]], function(x) x[[1]][lagi, 1 : 2], simplify = F)
        lens <- sapply(tlns, length)
        wlen <- which(lens == 0)
        if(length(wlen) > 0) {
            for(m in 1 : length(wlen)) {
                tlnest[[wlen[m]]] <- list()
                tlnest[[wlen[m]]][[1]] <- matrix(rep(NA, 6), ncol = 2)
            }
        }
        # #first, get only main results
        tlns <- sapply(tlnest[[1]], function(x) {
            if(!is.null(x)){
                x[[1]][lagi, 1 : 2]
            }else{ c(NA, NA)}})
        for(j in 1 : ncol(tlns)) {
            tlns[, j] <- 100 * (exp(tlns[, j] *  iqrs[j]) - 1)
        }
        
        tlns <- tlns[, keeps]
        
        lb <- tlns[1,] - 1.96 * tlns[2,]
        ub <- tlns[1,] + 1.96 * tlns[2,]
        
        stop <- start1 + length(keeps) - 1
        tlnout <- rbind(tlnout, cbind(tlns[1, ], lb, ub))
        tlnoutc <- rbind(tlnoutc, cbind(paste0("lag", 
                                               rep(lagi, length(keeps))),
                                        names1[keeps], 
                                        rep(namesTLN[i], length(keeps))))
        
        start <- stop + 1
    }
    
}
tlnoutall <- data.frame(tlnout, tlnoutc)
tlnoutall <- tlnoutall[-1, ]
colnames(tlnoutall) <- c("est", "lb", "ub", "lag", "source", "type")





# save(list = ls(), file = '~/Dropbox/SpatialFA/data/all_hosp_results.RData')










###########	
###########	
###########	
###########	
###########	


pd <- position_dodge(.4)
tln2 <- tlnoutall[complete.cases(tlnoutall), ]
lb1 <- -1
ub1 <- 2
tln2$lag <- factor(tln2$lag, levels = c("lag1", "lag2", "lag3"), 
                   labels = c("Lag 0", "Lag 1", "Lag 2"))
tln2$type <- factor(tln2$type, levels = c("APCA", "mAPCAsw"), 
                    labels = c("SHARE", "mAPCA"))
ubsind <- ifelse(tln2$ub > ub1, 1, 0)
lbsind <- ifelse(tln2$lb < lb1, 1, 0)

inds <- ifelse(ubsind == 0 & lbsind == 0, 0, ifelse(ubsind == 1 & lbsind == 0, 
                                                    1, ifelse(ubsind == 0 & lbsind == 1, 2, 3)))

ub2 <- ifelse(tln2$ub >ub1, Inf, tln2$ub)
lb2 <- ifelse(tln2$lb < lb1, -Inf, tln2$lb)

tln2 <- data.frame(tln2, inds, ub2, lb2)

limits <- aes(ymax = ub, ymin= lb)


tln3 <- tln2
names2 <- names1
# names2 <- names1[c(1, 2, 3, 5, 7, 4, 8, 6)]
# tln3 <- tln3[-which(tln3$source %in% names2[6 : 8]), ]
# names2 <- names2[1 : 5]
tln3$source <- factor(tln3$source, levels = names2)

s1 <- c("Metals", "Traffic", "Residual oil", "Soil", "Salt", "Sec. Sulfate")
# tln3 <- tln3[-which(tln3$source %in% s1), ]
tln3 <- tln3[which(tln3$source %in% s1), ]
tln3$source <- factor(tln3$source, levels = s1)
tlnhold <- tln3[which(tln3$type == "SHARE" & tln3$lag == "Lag 0"), ]
sources  <- as.character(tlnhold[order(tlnhold$est, 	
                                       decreasing = T), "source"])
tln3$source <- factor(tln3$source, levels = sources)


# tln3 <- tln2

size1 <- 18
sizep <- 1.2

col1 <- brewer.pal(5, "Blues")[3:5]
col1 <- brewer.pal(5, "Dark2")[c(1, 4)]
col1 <- brewer.pal(8, "Set1")[1:2]
col1 <- c("indianred3", col1[2])
gplot1 <- function(lag, dat = tln3) {
    dat <- dat[which(dat$lag %in% paste("Lag", lag)), ]
    
    g1 <- ggplot(dat, aes(x = source, y = est, colour = type), 
                 ylim = c(-0.01, 0.01)) + 
        geom_hline(aes(yintercept = 0), colour = "grey80", 
                   linetype = "dashed") +
        geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = type), 
                        width = 0.1, position = pd, size = sizep) +
        # scale_colour_hue(drop = FALSE, name="",
        # # breaks=c("lag0", "lag1", "lag2"),
        # labels=c("Lag 0", "Lag 1", "Lag 2"),
        # l=20, c = 150, h = 150) + 
        # scale_colour_brewer(palette="Blues", labels=c("Lag 0", "Lag 1", "Lag 2"), name = "") +   
        # scale_color_manual( labels=c("Lag 0", "Lag 1", "Lag 2"), name = "",
        # values = col1) + 	
        scale_color_manual( labels=c("SHARE", "mAPCA"), name = "",
                            values = col1) + 	
        ylab(expression(atop("% increase in CVD hospitalizations per", 
                             "IQR increase in source concentration"))) +
        xlab("") +
        # ggtitle("Mortality effects by\nsource apportionment method") +
        scale_y_continuous(limits=c(lb1, ub1))      +         
        theme_bw() +        
        theme(axis.text.y=element_text(size=size1)) +
        # theme(legend.text=element_text(size=size1), legend.position = c(0.9, 0.85),
        # legend.direction = "vertical")  +
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

setwd("/Users/jennakrall/Dropbox/SpatialFA/plots")
pdf("hosp_east_lag0.pdf", height = 7, width = 13)
gplot1(0)
graphics.off()


setwd("/Users/jennakrall/Dropbox/SpatialFA/plots")
pdf("hosp_east_lag012.pdf", height = 10, width = 7)
gplot1(c(0,1, 2))
graphics.off()


setwd("/Users/jennakrall/Dropbox/SpatialFA/plots")
pdf("hosp_east_6.pdf", height = 10, width = 10)	
g1
graphics.off() 	
# ###########	
###########	
###########	
###########	
###########	











# ###########	
###########	
###########	
###########	
###########	
#table 1
tab1 <- data.frame(xtab, iqrs)
colnames(tab1) <- c("Monitors", "Counties", "Major constituents", "IQR")
xtable(tab1)	





