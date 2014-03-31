#File to get SHARE information from Medicare data, apply mAPCA


dir1 <- "/Users/jennakrall/Dropbox/SpatialFA/data"
dircode <- "/Users/jennakrall/Dropbox/SpatialFA/rcode"


#source code
source(file.path(dircode, "load_file_spatialfa_17may13.R"))


load(file.path(dir1, "speciation_medicare.RData"))




########
# perform SHARE
data.rr <- lapply(datall, function(x) {
	x[, -c(1, 2)]
})
set.seed(70878)
pvd.east <- domatchSIM(restrict.data = data.rr, threli = 1, 
	usedata = T, thresang = pi/4) #thresang = 70 * pi/180)
share <- pvd.east[[3]]
Sources <- pvd.east[[4]]

# save(pvd.east, share, Sources, file = "pvdeast_medicare.RData")






#how many sites/source
nsources <- length(Sources[-which(is.infinite(Sources))])
vec <- vector()
for(i in 1 : nsources) {

		temp <- sapply(share, function(x) i %in% x)
		vec[i] <- length(which(temp == T))
		
}
#how many sources for each monitors
lens <- vector()
for(i in  1: length(share)) {
	lens[i] <- length(share[[i]])
}








########
# Get names of sources
reg <- round(pvd.east[[1]][[1]], 2)
reg <- pvd.east[[1]][[1]]
n1 <- vector()
n2 <- vector()
for(i in 1 : nsources) {
	regs <- sort(reg[, i], decreasing = T)
	nr <- names(regs)[which(regs  > 0.4)]
	n1[i] <- paste(nr, collapse = "/")
	regs <- rev(regs)
	nr <- names(regs)[which(regs  < -0.4)]
	n2[i] <- paste(nr, collapse = "/")
}

cbind(paste(n1, n2), vec)



plot(pvd.east[[1]][[6]][[2]])
round(pvd.east[[1]][[1]], 3)

                                                # vec 
 # [1,] " lead/zinc/manganese"                    "42"
 # [2,] "silicon/aluminum/titanium/calcium/iron " "79"
 # [3,] "ammonium_ion/sulfate/OC/selenium "       "74"
 # [4,] " strontium/potassium/copper"             "70"
 # [5,] " chlorine/sodium_ion/nitrate/bromine"    "49"
 # [6,] " phosphorus/vanadium"                    "37"
 # [7,] " nickel/iron/vanadium/manganese"         "37"
 # [8,] " arsenic/selenium/bromine"               "11"
 # [9,] "elemental_carbon/OC/iron "               "29"
names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", "Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")



















#restrict to lat/long
monsKEEP1 <- monsKEEP[, c(2, 1)]


############
# Create plot of source locations
setwd(file.path(home.dir, "plots"))
pdf("map_east_sources_medicare.pdf", height = 7, width = 11)

#set margins
par(mfrow = c(3, 3), mar = c(3.4, 2.5, 2, 0), 
	oma = c(2,2,0,0))

#Get rid of infinite sources
whinf <- which(is.infinite(Sources))
if(length(whinf) > 0) {
	Sources <- Sources[-whinf]
}

#set size
cex1 <- 2.2
cex2 <- 1.75


# save map
m1 <- map("state", fill = TRUE, col = "grey60",
	ylim = c(36, 45), xlim = c(-90, -69),plot = F)

# set colors
cols <- brewer.pal(8, "Dark2")


#for each source
if(length(Sources) > 0) {
	for(i in 1 : length(Sources)) {
		if(!is.infinite(Sources[i])) {
	
		#plot map
		plot(m1$x, m1$y, type = "l", xlim = c(-90, -69),
			ylim = c(36, 45), xlab = "", ylab = "", 
			cex.axis = 2, axes = F)
		
		#add name
		mtext(names1[i], cex = 2)
		
		#which monitors have source
		monsS <- sapply(share, whichS, num = i)
		monsS <- which(monsS == 1)
		print(length(monsS))
		
		#add points for monitors without source
		if(length(monsS) != nrow(monsKEEP1)) {
			points(monsKEEP1[-monsS,], col = cols[1],
			 pch = "+", cex = cex1)
		}
		
		#add points for monitors with source
		points(monsKEEP1[monsS, ], col = cols[2], pch = 16, cex = cex2)
		
		
		#add legend for 8
		if(i == 8) {
			par(xpd = T)
			legend(-87, 35.5, col = c(cols[2], cols[1]), 
				legend = c("Source present", "Source absent"), 
				pch = c(16, 3), cex = 1.1, border = NULL)
			par(xpd = F)	
		}
		}
	}

}


graphics.off()













#######
# Plot for thresholds
#create sequence of thresholds
allthres <- seq(pi/8, pi/2, length = 10)


pdf("map_east_sourcesorder_medicare_allthres.pdf", height = 7)
j <- 1
print(j)
par(mfrow = c(3, 3))

	
set.seed(9268)
#get large threshold to start
pvd.east <- domatchSIM(restrict.data = data.rr, threli = 1, 
	usedata = T, thresang = 1.5)
share <- pvd.east[[3]]
Sources <- pvd.east[[4]]
whinf <- which(is.infinite(Sources))
	
if(length(whinf) > 0) {
	Sources <- Sources[-whinf]
}

cols <- brewer.pal(8, "Dark2")
cols2 <- brewer.pal(8, "Set1")
cols <- c(cols, cols2, cols)
if(length(Sources) > 0) {
	
	#for each source
	for(i in 1 : length(Sources)) {
		if(!is.infinite(Sources[i])) {
			
		#plot source name
		plot(1, 1, type = "n", axes = F, ylab = "", xlab = "")
		mtext(paste0("source = ", names1[i]))		
		
		#for each threshold
		for(k in 1 : length(allthres)) {
			
			
	
			#map threshold
			map("state", fill = TRUE, col = "grey90",
				ylim = c(37, 45), xlim = c(-90, -69),
				mar = c(0, 2, 0, 2), oma = c(0,0,0,0))
			map.axes()
	
			#add threshold to map
			threang <- round(allthres[k] * 180 / pi, 2)
			mtext(paste0("thresh=", threang))
		
			# SHARE for threshold
			set.seed(9268)
			pvd.east <- domatchSIM(restrict.data = data.rr, 
				threli = 1, 
				usedata = T, thresang = allthres[k])
			share <- pvd.east[[3]]
			Sources <- pvd.east[[4]]	
					
					
			#get monitors with source		
			monsS <- sapply(share, whichS, num = i)
			monsS <- which(monsS == 1)
			print(length(monsS))
			
			
			#plot monitors
			if(length(monsS) != nrow(monsKEEP1)) {
				points(monsKEEP1[-monsS, ], 
					col = "grey40", pch = 1, cex = .8)
			}
			points(monsKEEP1[monsS, ], col = 1, 
				pch = 16, cex = 1)
			}
		}
	}

}

graphics.off()






