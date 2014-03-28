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





###################
# Save data
#substring for FIPS code
n1 <- substr(monsKEEP[, 4], 1, 5)
EAST1 <- data.frame(n1, monsKEEP)
colnames(EAST1) <- c("FIPS", "monname", "LATITUDE", "LONGITUDE", "monname1")
save(EAST1, file = file.path(dir1, "link_east_county.RData"))







################
# CHECK USING BOUNDARIES OF COUNTY
################

#######
# get associated fips with each monitor using boundaries

#load fips data
data(county.fips)
data(countyMapEnv)
namesall <- map("county", plot = F)$names




#for boundary associated with each county, get name associated
# assign fips code to each row
monsKEEP1 <- monsKEEP[, c(2, 3)]
EAST <- matrix(nrow = nrow(monsKEEP), ncol = 2)
for(i in 1 : length(namesall)) {
	# print(i)
	
	#get county boundaries for each county
	countdat <- matrix(c(map("county", namesall[i], 
		fill = TRUE, exact = T,
		plot = FALSE)$x,
		map("county", namesall[i], exact = T, fill = TRUE,
		plot = FALSE)$y), byrow = FALSE, ncol = 2)
	colnames(countdat) <- c("x", "y")	
		
	#if is.na values (separate parts of county)	
	if(length(which(is.na(countdat))) > 0) {
		
		#which is.na
		whna <- which(is.na(countdat[, 1]))
		io <- matrix(nrow = nrow(monsKEEP), 
			ncol = (length(whna) + 1))
		start <- 1
		stop <- whna[1] - 1
		
		#for each part of county, determine which monitors in
		for(k in 1 : (length(whna) + 1)) {
			
			#which monitors in the county?
			io[, k] <- 1 * inout(monsKEEP1, 
				countdat[start : stop, ], bound = T)
			start <- whna[k] + 1
			
			#update stop
			if(length(whna) > k) {
				stop <- whna[k + 1] - 1
			}else{
				stop <- nrow(countdat)
			}
		}
		io <- rowSums(io)
		
		
	#else determine which monitors in	
	}else{
		io <- 1 * (inout(monsKEEP1, countdat, bound = T))
	}
	
	
	#if at least one monitor in
	if(sum(io) > 0) {
		
		#which monitor in
		wh1 <- which(io == 1)
		if(!sum(is.na(EAST[wh1, 1]))) {
			print("error: already assigned")
			browser()
		}

		for(j in 1 : length(wh1)) {
			EAST[wh1[j], ] <- as.matrix(county.fips[which(
				county.fips[, 2] == namesall[i]), ])
		}
	}
	
}






#################
#Fix NA values
wna <- which(is.na(EAST[, 1]))




#fix 2, 3: in NEW YORK CITY/ in WATER
par(mfrow = c(1, 2))
m <- map("county", "new york", xlim = c(-74.5, -73.5), 
	ylim = c(40, 41), namesonly = T) 
points(monsKEEP[wna[2], c(2, 3)], col = "red")	
points(monsKEEP[wna[3], c(2, 3)], col = "blue")	
map("county", "new york", xlim = c(-74.5, -73.5), 
	ylim = c(40, 41))
for(i in 1 : length(m)) {
	map("county", m[i], add = T, col = i) 
}	
plot(seq(1, length(m)), type = "n")
for(i in 1 : length(m)) {
	abline(v = i, col = i)
}
	

EAST[69, ] <- as.matrix(county.fips[which(county.fips[, 2] == "new york,bronx"),])
#bronx matches FIPS from monitor code
EAST[65, ] <- as.matrix(county.fips[which(county.fips[, 2] == "new york,bronx"),])



#fix 1: CONNECTICUT
m <- map("county", "connecticut", namesonly = T)
par(mfrow = c(1, 2))
map("county", "connecticut")
points(monsKEEP[wna[1], c(2, 3)], col = "red")	
for(i in 1 : length(m)) {
	map("county", m[i], add = T, col = i)
}
plot(seq(1, length(m)), type = "n")
for(i in 1 : length(m)) {
	abline(v = i, col = i)
}
EAST[1, ] <- as.matrix(county.fips[which(county.fips[, 2] == m[1]), ])

par(mfrow = c(1, 1))
map("county", m[1], xlim = c(-73.4, -73.3), 
	ylim = c(41.05, 41.14))
points(monsKEEP[wna[1], c(2, 3)], col = "red")	







####
#check that fips code assigned for monitor makes sense
fips <- as.numeric(as.character(substr(monsKEEP[, 4], 1, 5)))
inouts <- vector(, length = nrow(monsKEEP))
for(i in 1 : length(fips)) {
	print(i)
	
	#cities within counties
	if(!(i %in% c(124, 125, 126))) {
	mon <- monsKEEP[i, c(2, 3)]
	
	#get name of fips
	n1 <- county.fips[which(county.fips[, 1] == fips[i]), 2]
	
	#find boundary of county
	countdat <- matrix(c(map("county", n1, 
		fill = TRUE, exact = T,
		plot = FALSE)$x,
		map("county", n1, exact = T, fill = TRUE,
		plot = FALSE)$y), byrow = FALSE, ncol = 2)
	colnames(countdat) <- c("x", "y")
	
	# determine whether monitor is in county
	inouts[i] <- inout(mon, countdat, bound = T)	
	}
}

whF <- which(inouts == F)
monsKEEP[whF[4:6],]


##########
#keep an eye on these... cities within larger counties
m <- map("county", "virginia", namesonly = T)
map("county", "virginia,washington")
points(monsKEEP[124, c(2, 3)])

map("county", "virginia,chesterfield")
points(monsKEEP[125, c(2, 3)])

map("county", "virginia,roanoke")
points(monsKEEP[126, c(2, 3)])

###########
# on border
map("county", c("west virginia"))
map("county", "ohio", add = T, col = "green")
map("county", "kentucky", add = T, col = "blue")
points(monsKEEP[whF[4:6],c( 2, 3)], col = "red")



EAST1 <- cbind(EAST, monsKEEP)
m1 <- as.numeric(as.character(EAST1[, 1]))
m2 <- as.numeric(substr(EAST1$temp, 1, 5))
all.equal(m1, m2)
EAST1[which(m1 != m2), c(1, 2)]
m1[which(m1 != m2)]
m2[which(m1 != m2)]







