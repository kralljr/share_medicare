#File to match monitors to counties

home.dir <- "~/Dropbox/SpatialFA"
save.dir <- file.path(home.dir, "rcode/medicare/sharesim/data")


load(file.path(home.dir, "data/speciation_medicare.RData"))


# get fips code from each monitor
fips.mons <- sapply(strsplit(monsKEEP[, 3], "\\."), 
                    function(x) x[1])
monsKEEP2 <- monsKEEP[, 3]
names(monsKEEP2) <- fips.mons
unmon <- unique(fips.mons)
unmonlist <- list()
for(i in 1 : length(unmon)) {
    unmonlist[[i]] <- monsKEEP2[which(names(monsKEEP2) == unmon[i])]
}
names(unmonlist) <- unmon

save(unmonlist, file = file.path(save.dir, "unmonlist.rda"))
