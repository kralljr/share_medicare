#File to get monitors for medicare analysis


#specify directories
dir1 <- "/Users/jennakrall/Dropbox/SpatialFA/data"
dircode <- "/Users/jennakrall/Dropbox/SpatialFA/rcode"


#get data
monitors <- readRDS("/Users/jennakrall/Dropbox/PM25cons_mort/healthest_epipaper_28may11/monitor-subset_all.rds")
counties <- readRDS("/Users/jennakrall/Dropbox/PM25cons_mort/cities/nmmaps_counties.rds")
specmons <- readRDS(file.path(dir1, "speciation_monitors.rds"))
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
monitors <- names(specmons)
mon2 <- unique(as.character(monsall[, 4]))
try1 <- monitors[-which(monitors %in% mon2)]
mon2[-which(mon2 %in% monitors)]

#get all monitors in NE
monsKEEP <- monsall[which(monsall[, 2] > x1 & 
	monsall[, 3] > y1), ]
monsKEEP <- monsKEEP[, c(2 : 4)]
monsKEEP <- unique(monsKEEP)
colnames(monsKEEP)[3] <- "monname"






###################
# Match monitors to medicare
nrow(monsKEEP)
monsKEEP <- monsKEEP[which(substr(monsKEEP$monname, 1, 5) %in% all.fips), ]










#################
# Find correct days
#for each monitors

keep <- vector(, length = nrow(monsKEEP))
datall <- list()
unmons <- list()
for(i in 1 : nrow(monsKEEP)) {
	
	#get speciation data
	monname <- as.character(monsKEEP[i, 3])
	dat <- specmons[[monname]]
	
	#limit to columns of interest
	dat <- get.cols(dat, monname)
	
	#get result formatted
	nc <- ncol(dat[[4]])
	temp <- data.frame(dat[[4]][, nc], 
		dat[[3]], dat[[4]][, -nc])
	colnames(datall)[1 : 2] <- c("Date", "PM2.5")
	montemp <- dat[[2]]
	
	#get complete cases
	cc <- complete.cases(temp)
	datall[[i]] <- temp[cc, ]
	unmons[[i]] <- montemp[cc]
	
	#ensure enough data
	if(class(dat) == "data.frame"){	
		if(nrow(datall[[i]]) > 50) {
			keep[i] <- 1
		}
	}
	
}
monsKEEP <- monsKEEP[which(keep == 1), ]
datall <- datall[which(keep == 1), ]
unmons <- unmons[which(keep == 1)]



save(monsKEEP, datall, unmons, file = "speciation_medicare.RData")








