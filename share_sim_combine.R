simout1 <- list()
for(i in 1 : 6) {
    load(paste0("simout", i, ".RData"))
    simout1[[i]] <- simout
}

library(xtable)
tabs <- t(sapply(simout1, function(x) x[[2]]))
rownames(tabs) <- c("25", "100", "5", "same", 
    "unequal subregions", "different days")
xtable(tabs)
