#File to get monitors for medicare analysis


#specify directories
dir1 <- "/Users/jennakrall/Dropbox/SpatialFA/data"
dircode <- "/Users/jennakrall/Dropbox/SpatialFA/rcode"
home.dir <- "/Users/jennakrall/Dropbox/SpatialFA"

#get data
# monitors <- readRDS("/Users/jennakrall/Dropbox/PM25cons_mort/healthest_epipaper_28may11/monitor-subset_all.rds")
specmons <- readRDS(file.path(dir1, "speciation_monitors.rds"))
monitors <- readRDS(file.path(dir1, "monitor_locations.rds"))
load(file.path(dir1, "all_fips_medicare.rds"))




#source functions (monsall is created here)
source(file.path(dircode, "load_file_spatialfa_17may13.R"))


#load libraries
library(splancs)
library(gpclib)
library(maps)



########
# Define region of interest
map("state")
axis(1)
axis(2)
points(unique(monsall[, c(2, 3)]))


x1 <- -89
y1 <- 36.5
abline(v = x1, col = "red")
abline(h = y1, col = "red")






###################
# Obtain location information for monitors
monitors2 <- names(specmons)
mon2 <- unique(as.character(monitors[, 3]))
try1 <- monitors2[-which(monitors2 %in% mon2)]
# "32003.0540"???

#really huge
try2 <- mon2[-which(mon2 %in% monitors2)]

#get all monitors in NE
monsKEEP <- monitors[which(monitors[, 2] > x1 & 
	monitors[, 1] > y1), ]
monsKEEP <- unique(monsKEEP)
colnames(monsKEEP)[3] <- "monname"


#restrict to speciation data
monsKEEP <- monsKEEP[which(monsKEEP[, 3] %in% monitors2), ]


#double check by plotting
map("state")
points(monsKEEP[, 2: 1], col = "blue", pch = 16)
abline(v = x1, col = "red")
abline(h = y1, col = "red")





###################
# Match monitors to medicare
monsKEEP <- monsKEEP[which(substr(monsKEEP$monname, 1, 5) %in% all.fips), ]










#################
# Find correct days
#for each monitors

#all columns are the same, get cols to keep
cn <- colnames(specmons[[1]])
cn1 <- strsplit(cn, "\\.")
cntype <- sapply(cn1, function(x) x[2])
cnname <- sapply(cn1, function(x) x[1])
#limit to columns of interest
whCN <- which(is.na(cntype) & cnname %in%  cons)
whPM <- which(is.na(cntype) & cnname == "PM25_SPEC")
whCN <- c(1, whPM, whCN)

keep <- vector(, length = nrow(monsKEEP))
datall <- list()
unmons <- list()
for(i in 1 : nrow(monsKEEP)) {
	
	#get speciation data
	monname <- as.character(monsKEEP[i, 3])
	dat <- specmons[[monname]]
	
	#limit to columns of interest
	dat <- get.cols2(dat, monname, whCN)
	
	#get result formatted
	temp <- dat[[1]]
	montemp <- dat[[2]]
	
	#get complete cases
	cc <- complete.cases(temp)
	datall[[i]] <- temp[cc, ]
	unmons[[i]] <- montemp[cc]
	
	#ensure enough data
	if(class(datall[[i]]) == "data.frame"){	
		if(nrow(datall[[i]]) > 50) {
			keep[i] <- 1
		}
	}
	
}
monsKEEP <- monsKEEP[which(keep == 1), ]
datall <- datall[which(keep == 1)]
unmons <- unmons[which(keep == 1)]



save(monsKEEP, datall, unmons, file = "speciation_medicare.RData")







#########
# Plot map
cols <- brewer.pal(8, "Dark2")


setwd(file.path(home.dir, "plots", "region_maps"))

pdf("map_east_medicare.pdf", height = 7)
map("state", col = "grey90",
	ylim = c(36, 45), xlim = c(-90, -69), fill = T, 
	mar = c(0, 5, 0, 0), oma = c(0,0,0,0))
map.axes()	


points(monsKEEP[, c(2, 1)], col = "grey10", lwd = 2)
graphics.off()



