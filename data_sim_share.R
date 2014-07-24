#######
#######
#######
#######
# Create simulation data
#######
#######
#######
#######


basedir <- "C:/Users/jrkrall/Dropbox"
home.dir <- file.path(basedir, "SpatialFA")
dir1 <- file.path(home.dir, "data")
load(file.path(dir1, "speciation_medicare.RData"))


#use Allegheny County, PA data 
#(city with large pollution, many days of data)
dat <- datall[[61]][, -c(1, 2)]


# Get principal components from data
stddat <- stdize1(dat)$dat
pr <- prcomp(stddat, retx = T)
nc <- length(which(pr$sd > 1))
vec <- as.matrix(varimax(pr$rot[, 1 : nc])$load[1: 24, ])


#Get positive part of PCs


cs <- colSums(as.matrix(vec))
#Take opposite sign if sum is negative
if(length(which(cs < 0)) > 0) {
    vec[, which(cs < 0)] <- -vec[, which(cs < 0)]
}
#Get positive part
for(i in 1 : ncol(vec)) {
    wh0 <- which(vec[, i] < 0)
    if(length(wh0) > 0) {
        vec[wh0, i] <- rep(0, length(wh0))
    }
}
#Normalize so columns sum to 1
vec <- sweep(vec, 2, colSums(vec), "/")


#Name sources
names <- c("traffic", "fireworks", "soil", 
           "sec sulf", "salt", "metals", "P/V")

cms <- rep(mnvar[1, ], 2)
sds <- rep(sqrt(mnvar[2, ]), 2)


#set up 5 subregions
keeps <- list()
#with all
keeps[[1]] <- seq(1, nc)
#minus 2
keeps[[2]] <- seq(1, nc)[-c(1, 2)]
#minus 3 
keeps[[3]] <- seq(1, nc)[-c(1, 2, 3)]
#minus 3
keeps[[4]] <- seq(1, nc)[-c(1, 2, 4)]
#minus 3
keeps[[5]] <- seq(1, nc)[-c(1, 2, 5)]


save(dat, names, vec, keeps, cms, sds, 
    file = file.path(dir1, "data_share_sim.RData"))