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




########
# get DATA
load(file.path(dir1, "speciation_medicare.RData"))
data.rr <- lapply(datall, function(x) {
    x[, -c(2)]
})
pm25 <- lapply(datall, function(x) x[, 2])
names(data.rr) <- monsKEEP[, 3]






########
# Perform SHARE and mAPCA

share <- sharehealth(data.rr, tots = pm25, list = unmonlist)
mapca <- sharehealth(data.rr, tots = pm25, list = unmonlist, method = "mapca")




#######
# Match results between mAPCA and SHARE


nc <- ncol(data.rr[[1]][, -1])
mapcaload <- mapca$major.sig
#use high threshold so match all sources
match1 <- matchfun(list(mapcaload), share$major.sig, thres = 70 * pi/180)$match[[1]]
rownames(match1) <- paste0("mapca", seq(1, ncol(mapcaload)))
colnames(match1) <- paste0("share", seq(1, ncol(share$major.sig)))
match1 <- whichCS(t(match1))


#reorder mapca
mapca$major.sig <- suppressWarnings(mapca$major.sig[, match1])
mapcaconc <- sapply(mapca$sources, function(x) {
    dat <- as.matrix(x[, -1])
    dat <- data.frame(x[, 1], suppressWarnings(dat[, match1]))
    colnames(dat) <- c("date", paste0("source", seq(1, ncol(dat) - 1)))
    dat
    }, simplify = F)
names(mapcaconc) <- names(data.rr)













########
# Results table
#Number of days
days <- sapply(data.rr, nrow)
c(min(days), max(days), median(days))

simpleCap <- function(x) {
    s1 <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s1, 1, 1)), substring(s1, 2),
          sep = "", collapse = " ")
}

names <- c("Metals", "Soil", "Sec. Sulfate", "Fireworks", 
           "Salt", "P/V", "Residual oil", "As/Se/Br", "Traffic")


summ <- data.frame(names, share$summary[, c("monitor", "counties", "cons", "IQR")])
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
colnames(summ) <- c("Sources", "Monitors", "Counties", "Major constituents", "IQR")

print(xtable(summ), include.rownames = F)






########
# Results map



#set map defaults

#set size, bw
cex1 <- 2.2
cex2 <- 1.75
lwd1 <- 1
cols <- c(1, "grey50")

#set size, color
lwd1 <- 2
cex1 <- 3
cex2 <- 2
cols <- brewer.pal(8, "Dark2")
#cols <- c("darkolivegreen3", "orange")

m1 <- map("state", fill = TRUE, col = "grey60",
          ylim = c(36, 45), xlim = c(-90, -69),plot = F)

pdf(file.path(paperdir, "map_east_sources.pdf"), height = 7, width = 11)
par(mfrow = c(3, 3), mar = c(4, 2.5, 1.4, 0),  oma = c(3,2,1,0))
reps <- seq(1, 9)
reps2 <- 0
for(j in c(1, 4, 7)) {
    reps3 <- rep(reps[j : (j + 2)], 2)
    reps2 <- c(reps2, reps3)
}
reps2 <- reps2[-1]
Sources <- sort(unique(unlist(share$share)))
Sources <- Sources[!is.infinite(Sources)]
#layout(matrix(c(reps2, rep(10, 3)), 7, 3, byrow = TRUE))
for(i in Sources) {
    plot(m1$x, m1$y, type = "l", xlim = c(-90, -69),
         ylim = c(36, 45), xlab = "", ylab = "", 
         cex.axis = 2, axes = F)
    keeps <- sapply(share$share, function(x) {
        ifelse(i %in% x, 1, 0)
    })
    wh1 <- which(keeps == 1)
    mtext(names[i], cex = 2)
    points(monsKEEP[-wh1, c(2, 1)],  col = cols[1], pch = "+", cex = cex1)
    points(monsKEEP[wh1, c(2, 1)],  col = cols[2], pch = 16, cex = cex2)
    if(i == 8) {
        
        par(xpd = T)
        legend(-87, 35.5, col = c(cols[2], cols[1]), 
               legend = c("Source found", "Source not found"), 
               pch = c(16, 3), cex = 1.5, border = NULL, bty = "n")   
        par(xpd = F)

    }
}

graphics.off()
	