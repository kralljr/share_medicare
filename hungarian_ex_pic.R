### File to create hungarian method example
library(clue)
library(xtable)
setwd("~/Dropbox/SpatialFA/paper_spatialfa/")

### Run east_share.R to get SHARE results

ms <- t(share$major.sig)
#ms <- t(apca(data.rr[[2]])$vmax[[1]][1 : 24, ])
mon1 <- data.rr[[2]]
apca1 <- t(apca(mon1)$vmax[[1]][1 : 24, ])

# Get max direction
mins <- sign(apply(ms, 1, sum))
ms <- sweep(ms, 1, mins, "*")
# Get max direction
mins <- sign(apply(apca1, 1, sum))
apca1 <- sweep(apca1, 1, mins, "*")

#standardize each
norm1 <- function(x) sqrt(sum(x^2))
ms <- (sweep(ms, 1, apply(ms, 1, norm1), "/"))
apca1 <-  (sweep(apca1, 1, apply(apca1, 1, norm1), "/"))

# get psi_i
psii <- abs(apca1 %*% t(ms))
# get angles
ang1 <- round(acos(psii) * 180 / pi, 1)


ang2 <- ang1[c(2, 3, 4, 5), -c(1, 3, 5, 9)]
colnames(ang2) <- paste0("M", seq(1, ncol(ang2)))
rownames(ang2) <- paste0("i", seq(1, nrow(ang2)))

ls <- as.vector(solve_LSAP(ang2))
ls <- cbind(seq(1, nrow(ang2)), ls)
ls[-which(ang2[ls] > 45), ]


hatpsi <- diag(x = 0, nrow = nrow(ang2), ncol = ncol(ang2))
hatpsi[ls] <- 1
colnames(hatpsi) <- paste0("M", seq(1, ncol(ang2)))
rownames(hatpsi) <- paste0("i", seq(1, nrow(ang2)))

write.csv(round(psii[c(2, 3, 4, 5), -c(1, 3, 5, 9)], 2), file = "psii.csv")
write.csv(hatpsi, file = "hatpsi.csv")
write.csv(ang2, file = "hun_angles.csv")
