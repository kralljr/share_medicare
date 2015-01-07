# File to estimate major sources in NE US


#Set directory
home.dir <- "~/Dropbox/SpatialFA/"
dir1 <- file.path(home.dir, "data")
paperdir <- file.path(home.dir, "paper_spatialfa/figs")

#Load my packages
library(share)
library(handles)
library(sharesim)

#load other packages
library(devtools)
library(RColorBrewer)
library(xtable)
library(maps)
library(dplyr)



########
# get DATA
load(file.path(dir1, "speciation_medicare.RData"))
data.rr <- lapply(datall, function(x) {
    x[, -c(2)]
})
pm25 <- lapply(datall, function(x) x[, 2])
names(data.rr) <- monsKEEP[, 3]




# get data specifically for each season

# indicator function for warm (ito 2010)
iswarm <- function(dat) {
	months <- as.numeric(substr(dat[, 1], 6, 7))
	ifelse(between(months, 4, 9), 1, 0)
}
warm <- lapply(data.rr, iswarm)

cold.dat <- list()
warm.dat <- list()
cold.pm <- list()
warm.pm <- list()
for(i in 1 : length(data.rr)) {
	print(i)
	warm1 <- warm[[i]]
	cold.dat[[i]] <- data.rr[[i]][warm1 == 0, ]
	warm.dat[[i]] <- data.rr[[i]][warm1 == 1, ]
	
	cold.pm[[i]] <- pm25[[i]][warm1 == 0]
	warm.pm[[i]] <- pm25[[i]][warm1 == 1]
}




# restrict to monitors with more than 50 days
cold.great <- which(sapply(cold.pm, length) >= 50)
cold.pm <- cold.pm[cold.great]
cold.dat <- cold.dat[cold.great]
names(cold.dat) <- names(data.rr)[cold.great]
names(cold.pm) <- names(data.rr)[cold.great]


warm.great <- which(sapply(warm.pm, length) >= 50)
warm.pm <- warm.pm[warm.great]
warm.dat <- warm.dat[warm.great]
names(warm.dat) <- names(data.rr)[warm.great]
names(warm.pm) <- names(data.rr)[warm.great]


#get new list of monitors
names.cold <- names(data.rr)[cold.great]
names(names.cold) <- substr(names.cold, 1, 5)
un.cold <- unique(substr(names.cold, 1, 5))
cold.mons <- list()
for(i in 1 : length(un.cold)) {
	cold.mons[[i]] <- names.cold[which(names(names.cold) == un.cold[i])]
}
names(cold.mons) <- un.cold

names.warm <- names(data.rr)[warm.great]
names(names.warm) <- substr(names.warm, 1, 5)
un.warm <- unique(substr(names.warm, 1, 5))
warm.mons <- list()
for(i in 1 : length(un.warm)) {
	warm.mons[[i]] <- names.warm[which(names(names.warm) == un.warm[i])]
}
names(warm.mons) <- un.warm





########
# Perform SHARE and mAPCA

set.seed(10)

share.warm <- sharehealth(warm.dat, tots = warm.pm, list = warm.mons)
share.cold <- sharehealth(cold.dat, tots = cold.pm, list = cold.mons)

set.seed(10)

mapca.warm <- sharehealth(warm.dat, tots = warm.pm, list = warm.mons, method = "mapca")
mapca.cold <- sharehealth(cold.dat, tots = cold.pm, list = cold.mons, method = "mapca")




#######
# Match results between mAPCA and SHARE

#COLD

nc <- ncol(cold.dat[[1]][, -1])
mapcaload <- mapca.cold$major.sig
#use high threshold so match all sources
match1 <- matchfun(list(mapcaload), share.cold$major.sig, thres = 70 * pi/180)$match[[1]]
rownames(match1) <- paste0("mapca", seq(1, ncol(mapcaload)))
colnames(match1) <- paste0("share", seq(1, ncol(share.cold$major.sig)))
match1 <- whichCS(t(match1))


#reorder mapca
mapca.cold$major.sig <- suppressWarnings(mapca.cold$major.sig[, match1])
mapcaconc.cold <- sapply(mapca.cold$sources, function(x) {
    dat <- as.matrix(x[, -1])
    dat <- data.frame(x[, 1], suppressWarnings(dat[, match1]))
    colnames(dat) <- c("date", paste0("source", seq(1, ncol(dat) - 1)))
    dat
    }, simplify = F)
names(mapcaconc.cold) <- names(cold.dat)



#WARM

nc <- ncol(warm.dat[[1]][, -1])
mapcaload <- mapca.warm$major.sig
#use high threshold so match all sources
match1 <- matchfun(list(mapcaload), share.warm$major.sig, thres = 70 * pi/180)$match[[1]]
rownames(match1) <- paste0("mapca", seq(1, ncol(mapcaload)))
colnames(match1) <- paste0("share", seq(1, ncol(share.warm$major.sig)))
match1 <- whichCS(t(match1))


#reorder mapca
mapca.warm$major.sig <- suppressWarnings(mapca.warm$major.sig[, match1])
mapcaconc.warm <- sapply(mapca.warm$sources, function(x) {
    dat <- as.matrix(x[, -1])
    dat <- data.frame(x[, 1], suppressWarnings(dat[, match1]))
    colnames(dat) <- c("date", paste0("source", seq(1, ncol(dat) - 1)))
    dat
    }, simplify = F)
names(mapcaconc.warm) <- names(warm.dat)













########
# Results table
#Number of days
days <- sapply(cold.dat, nrow)
c(min(days), max(days), median(days))

days <- sapply(warm.dat, nrow)
c(min(days), max(days), median(days))

simpleCap <- function(x) {
    s1 <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s1, 1, 1)), substring(s1, 2),
          sep = "", collapse = " ")
}

names1 <- rownames(share.cold$major.sig)
apply(round(share.cold$major.sig, 3), 2, function(x) names1[order(abs(x), decreasing = T)[1 : 5]])

names.cold <- c("Traffic", "Sec. Sulfate", "Soil", "Residual oil", "Salt", "As/Se/Br", "P", "Metals")


names1 <- rownames(share.warm$major.sig)
apply(round(share.warm$major.sig, 3), 2, function(x) names1[order(abs(x), decreasing = T)[1 : 5]])

names.warm <- c("Sec. Sulfate", "Soil", "Fireworks", "Metals", "Salt", "P/V", "Residual oil", "As/Se/Br")














#### Table for cold
summ <- data.frame(names.cold, share.cold$summary[, c("monitor", "counties", "cons", "Median", "IQR")])
summ$IQR <- round(summ$IQR, 2)
temp <- sapply(as.character(summ$cons), function(x) {
    x <- strsplit(x, ", ")[[1]]
    for(i in 1 : length(x)) {
        x1 <- x[i]
        x1 <- simpleCap(as.character(x1))
        x1 <- ifelse(x1 == "Elemental_carbon", "EC", x1)
        x1 <- ifelse(x1 == "Sodium_ion", "Sodium ion", x1)
        x1 <- ifelse(x1 == "Ammonium_ion", "Ammonium", x1)
        x[i] <- x1
    }
    paste(x, collapse = ", ")
    }, simplify = T)
summ$cons <- temp

colnames(summ) <- c("Sources", "Monitors", "Counties", "Major constituents", "Median", "IQR")

summ$Monitors <- as.integer(summ$Monitors)
summ$Counties <- as.integer(summ$Counties)

# colnames(summ)[-1] <- paste0(colnames(summ)[-1], ".cold")
summ[, 1] <- as.character(summ[, 1])
summ[which(summ[, 1] == "P"), 1] <- "P/V"
summ.cold <- summ








##### TABLE FOR WARM


summ <- data.frame(names.warm, share.warm$summary[, c("monitor", "counties", "cons", "Median", "IQR")])
summ$IQR <- round(summ$IQR, 2)
temp <- sapply(as.character(summ$cons), function(x) {
    x <- strsplit(x, ", ")[[1]]
    for(i in 1 : length(x)) {
        x1 <- x[i]
        x1 <- simpleCap(as.character(x1))
        x1 <- ifelse(x1 == "Elemental_carbon", "EC", x1)
        x1 <- ifelse(x1 == "Sodium_ion", "Sodium ion", x1)
        x1 <- ifelse(x1 == "Ammonium_ion", "Ammonium", x1)
        x[i] <- x1
    }
    paste(x, collapse = ", ")
    }, simplify = T)
summ$cons <- temp

colnames(summ) <- c("Sources", "Monitors", "Counties", "Major constituents", "Median", "IQR")

summ$Monitors <- as.integer(summ$Monitors)
summ$Counties <- as.integer(summ$Counties)

# colnames(summ)[-1] <- paste0(colnames(summ)[-1], ".warm")

summ <- summ[match(names1, summ[, 1]), ]
summ.cold <- summ.cold[match(names1, summ.cold[, 1]), ]

print(xtable(summ.cold), include.rownames = F)
print(xtable(summ), include.rownames = F)

summ <- merge(summ, summ.cold, all = T)








########
# Results map



#set map defaults

#set size, bw
cex1 <- 2.2
cex2 <- 1.75
lwd1 <- 1
cols <- c(1, "grey50")

#set size, color
# lwd1 <- 2
# cex1 <- 3
# cex2 <- 2
# cols <- brewer.pal(8, "Dark2")
#cols <- c("darkolivegreen3", "orange")

m1 <- map("state", fill = TRUE, col = "grey60",
          ylim = c(36, 45), xlim = c(-90, -69),plot = F)
          
          
names1 <- c("Metals", "Soil", "Sec. Sulfate", "Fireworks", "Salt", "P/V", "Residual oil", "As/Se/Br", "Traffic")



pdf(file.path(paperdir, "map_east_sources_warm.pdf"), height = 7, width = 11)
par(mfrow = c(3, 3), mar = c(4, 2.5, 1.4, 0),  oma = c(3,2,1,0))
reps <- seq(1, 9)
reps2 <- 0
for(j in c(1, 4, 7)) {
    reps3 <- rep(reps[j : (j + 2)], 2)
    reps2 <- c(reps2, reps3)
}
reps2 <- reps2[-1]
Sources <- sort(unique(unlist(share.warm$share)))
Sources <- Sources[!is.infinite(Sources)]
#layout(matrix(c(reps2, rep(10, 3)), 7, 3, byrow = TRUE))
l <- 1
for(k in names1) {
	
	if(k %in% names.warm) {
	
	i <- which(names.warm == k)
	
    plot(m1$x, m1$y, type = "l", xlim = c(-90, -69),
         ylim = c(36, 45), xlab = "", ylab = "", 
         cex.axis = 2, axes = F, col = "grey70")
    keeps <- sapply(share.warm$share, function(x) {
        ifelse(i %in% x, 1, 0)
    })
    wh1 <- which(keeps == 1)
    mtext(names.warm[i], cex = 2)
    points(monsKEEP[-wh1, c(2, 1)],  col = cols[1], pch = "+", cex = cex1)
    points(monsKEEP[wh1, c(2, 1)],  col = cols[2], pch = 16, cex = cex2)

    }else{
    	
    	plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
    	}
    	
    	    if(l == 8) {
        
        par(xpd = T)
        legend(-87, 35.5, col = c(cols[2], cols[1]), 
               legend = c("Source found", "Source not found"), 
               pch = c(16, 3), cex = 1.5, border = NULL, bty = "n")   
        par(xpd = F)

    }
    l <- l + 1
}

graphics.off()







pdf(file.path(paperdir, "map_east_sources_cold.pdf"), height = 7, width = 11)
par(mfrow = c(3, 3), mar = c(4, 2.5, 1.4, 0),  oma = c(3,2,1,0))
reps <- seq(1, 9)
reps2 <- 0
for(j in c(1, 4, 7)) {
    reps3 <- rep(reps[j : (j + 2)], 2)
    reps2 <- c(reps2, reps3)
}
reps2 <- reps2[-1]
Sources <- sort(unique(unlist(share.cold$share)))
Sources <- Sources[!is.infinite(Sources)]
#layout(matrix(c(reps2, rep(10, 3)), 7, 3, byrow = TRUE))
l <- 1
names.cold[7] <- "P/V"
for(k in names1) {
	
	if(k %in% names.cold) {
	
	i <- which(names.cold == k)
	
	
    plot(m1$x, m1$y, type = "l", xlim = c(-90, -69),
         ylim = c(36, 45), xlab = "", ylab = "", 
         cex.axis = 2, axes = F, col = "grey70")
    keeps <- sapply(share.cold$share, function(x) {
        ifelse(i %in% x, 1, 0)
    })
    wh1 <- which(keeps == 1)
    mtext(names.cold[i], cex = 2)
    points(monsKEEP[-wh1, c(2, 1)],  col = cols[1], pch = "+", cex = cex1)
    points(monsKEEP[wh1, c(2, 1)],  col = cols[2], pch = 16, cex = cex2)

        }else{
    	
    	plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
    	}
    	
    	    if(l == 8) {
        
        par(xpd = T)
        legend(-87, 35.5, col = c(cols[2], cols[1]), 
               legend = c("Source found", "Source not found"), 
               pch = c(16, 3), cex = 1.5, border = NULL, bty = "n")   
        par(xpd = F)

    }
    l <- l + 1
}

graphics.off()
	