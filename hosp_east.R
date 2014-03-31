#####  File to get associations between hospitalizations and sources

# set working directories
dircode <- ("/Users/jennakrall/Dropbox/SpatialFA/rcode")
home.dir <- "/Users/jennakrall/Dropbox/SpatialFA"


#set working directory
setwd(home.dir)


#get data
load("/Users/jennakrall/Dropbox/PM25cons_mort/FinalR_815/Healtheffectsdata_17aug11.RData")
# load("/Users/jennakrall/Dropbox/PM25cons_mort/FinalR_815/Prednaivedata_17aug11.RData")
# spec <- readRDS("~/Dropbox/SpatialFA/data/speciation_monitors.rds")
counties <- readRDS("/Users/jennakrall/Dropbox/PM25cons_mort/cities/nmmaps_counties.rds")
load(file.path(home.dir, "data/speciation_medicare.RData"))
load(file.path(home.dir, "data/pvdeast_medicare.RData"))
load(file.path(home.dir, "data/mort_sources_medicare.RData"))
load("/Users/jennakrall/Dropbox/SpatialFA/data/sourceconc_medicare.RData")

#load libraries
library(splancs)
library(gpclib)
library(splines)
library(tsModel)
library(tlnise)
library(chron)
library(timeDate)
library(ggplot2)


#load code
source(file.path(dircode, "functions/mort_fn_east.R"))
source(file.path(dircode, "load_file_spatialfa_17may13.R"))







#######
#first create list of same county monitors
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
#######
	


########
########
########
########
########
########
########
########
# get source concentrations
#######


set.seed(949069)
#####


#get data minus date, PM2.5
data.rr <- lapply(datall, function(x) {
	x[, -c(1, 2)]
})
pm25 <- lapply(datall, function(x) x[, 2])
dates <- lapply(datall, function(x) x[, 1])

#First, do mAPCA
mAPCAall <- spatial.apca(dat = data.rr, 
	monitors = monsKEEP[, 3], lim = 50, tots = pm25, 
	dates = dates)	
mAPCA <- mAPCAall[[1]][[2]]

	
	
#get APCA/SHARE results, format APCA

#set up
APCA <- list()
rms <- 0
for(i in 1 : length(datall)) {
	dat <- data.rr[[i]]
	dates <- datall[[i]][, 1]
	
	#if at least 50 days of data
	if(nrow(dat) >= 50) {
		print(i)
		
		#names for apca based on share
		cn <- c("Date", paste0("source", i, ".", 
			share[[i]]))
		#number of sources
		nf1 <- length(share[[i]])


		#get APCA results
		APCA[[i]] <- getsources(dat, nf1, type = "apca", 
			tots = pm25[[i]])
		#format output
		APCA[[i]] <- data.frame(dates, APCA[[i]])
		colnames(APCA[[i]]) <- cn
		
		#format mAPCA output
		mAPCA[[i]] <- data.frame(dates, mAPCA[[i]])
		colnames(mAPCA[[i]]) <- c("Date", 
			paste0("source", i, ".", seq(1, ncol(mAPCA[[i]]) - 1)))
		
	}else{
		#which to remove
		rms <- c(rms, i)
		}
}
rms <- rms[-1]
print(rms)


#remove any missing
if(length(rms) > 0) {
	APCA <- APCA[-rms]
	mAPCA <- mAPCA[-rms]
}
#errors okay-- replacing Inf matches with NA










######
#create share info for mAPCA
#no date
share1 <- seq(1, ncol(mAPCA[[1]]) - 1)
sharesmAPCA <- list()
for(i in 1 : nrow(monsKEEP)) {sharesmAPCA[[i]] <- share1}












######
#average source conc for those monitors in same city,
#get final source conc




	
#combine source info for each monitor
mAPCA1 <- getallsource(mAPCA, sharesmAPCA, unmonlist, monsKEEP2)
APCA1 <- getallsource(APCA, share, unmonlist, monsKEEP2 )



#compare number of sources between APCA, mAPCA
SourcesmAPCA <- getS(mAPCA1)
SourcesAPCA <- getS(APCA1)
length(SourcesmAPCA)
length(SourcesAPCA)
Sources <- SourcesAPCA


names(mAPCA1[[1]]) <- names(unmonlist)
names(APCA1[[1]]) <- names(unmonlist)
save(APCA1, mAPCA1, SourcesmAPCA, SourcesAPCA, 
	file = "sourceconc_medicare.RData")















########
########
########
########
########
########
########
########
#get mortality risk for each city

#Medicare data on enigma, run on there

















#########
# get TLNise results for all counties
names(mAPCA1[[2]]) <- names(unmonlist)
names(APCA1[[2]]) <- names(unmonlist)
tlnmAPCA <- tlncomb(mortsmAPCA, SourcesmAPCA, mAPCA1[[2]],
	 "tln_mapca_est.txt", names = names(unmonlist))	
tlnAPCA <- tlncomb(mortsAPCA, Sources, APCA1[[2]], "tln_apca_est.txt",
	names = names(unmonlist))	

sapply(tlnAPCA[[2]], length)

readTLNconv("tln_mapca_est.txt")
readTLNconv("tln_apca_est.txt")

#all converged












########
#
#which cities have which sources
SourcesCit <- list()
tln <- tlnAPCA[[1]]
for(lag in 1 : 3) {
for(i in 1 : length(Sources)) {
	print(i)
	if(lag == 1) {
		SourcesCit[[i]] <- rownames(tln[[i]][[2]][[lag]])
	}else{
		print(all.equal(SourcesCit[[i]], rownames(tln[[i]][[2]][[lag]])))
		}
}
}

SourcesCitAPCA <- SourcesCit








######
# switch mAPCA results
mS <- mAPCAall[[6]][[2]]$load[1:24, ]
oS <- pvd.east[[1]][[1]]
switch <- whichCS(find.local.reg(list(oS), mS, thres = 1)[[1]][[1]])
whInf <- which(is.infinite(switch))
if(length(whInf) > 0) {
	switch[whInf] <- NA
}


tlnmAPCAsw <- sapply(tlnmAPCA, function(x) x[switch], simplify = F)



















########
########
########
########
########
########
########
########
### PLOT RESULTS





### first, get data
#check which sources, name
round(pvd.east[[1]][[1]], 2)
names1 <- c("Metals","Soil", "Sec. Sulfate","Fireworks", 
	"Salt", "P/V",  "Residual oil",  "As/Se/Br", "Traffic")



keeps <- seq(1, length(Sources))
bounds <- c(-10, 10)
namesTLN <- c("APCA", "mAPCAsw")
pchs <- c(16, 17, 15)
cols <- c(1, "grey60", "grey20")


skips1 <- c(0, 0.1, 0.2)
start1 <- 1
tlnout <- matrix(nrow = 1, ncol = 3)
tlnoutc <- matrix(nrow = 1, ncol = 3)
#plot individual mortalities
for(lagi in 1 : 3) {
	
	
	xstart <- seq(1, length(keeps)) + skips1[lagi]
	
	skips <- c(-.2, 0.2)
	is1 <- length(namesTLN)
	for(i in 1 : is1) {
		print(c(lagi, i))
		xs <- xstart + skips[i]
		
		tlnest <- get(paste0("tln", namesTLN[i]))
		
		tlns <- sapply(tlnest[[1]], function(x) x[[1]][lagi, 1 : 2], simplify = F)
		lens <- sapply(tlns, length)
		wlen <- which(lens == 0)
		if(length(wlen) > 0) {
			for(m in 1 : length(wlen)) {
				tlnest[[wlen[m]]] <- list()
				tlnest[[wlen[m]]][[1]] <- matrix(rep(NA, 6), ncol = 2)
			}
		}
		tlns <- sapply(tlnest[[1]], function(x) {
			if(!is.null(x)){
				percinc(x[[1]][lagi, 1 : 2])
			}else{ c(NA, NA)}})
		tlns <- tlns[, keeps]
	
		lb <- tlns[1,] - 1.96 * tlns[2,]
		ub <- tlns[1,] + 1.96 * tlns[2,]
		
		stop <- start1 + length(keeps) - 1
		tlnout <- rbind(tlnout, cbind(tlns[1, ], lb, ub))
		tlnoutc <- rbind(tlnoutc, cbind(paste0("lag", 
			rep(lagi, length(keeps))),
			names1[keeps], 
			rep(namesTLN[i], length(keeps))))

		start <- stop + 1
	}

}
tlnoutall <- data.frame(tlnout, tlnoutc)
tlnoutall <- tlnoutall[-1, ]
colnames(tlnoutall) <- c("est", "lb", "ub", "lag", "source", "type")

	
	


	
	
	


	
	
	



###########	
###########	
###########	
###########	
###########	


pd <- position_dodge(.4)
tln2 <- tlnoutall[complete.cases(tlnoutall), ]
lb1 <- -10
ub1 <- 11
tln2$type <- factor(tln2$type, levels = c("APCA", "mAPCAsw"), 
	labels = c("SHARE", "mAPCA"))
ubsind <- ifelse(tln2$ub > ub1, 1, 0)
lbsind <- ifelse(tln2$lb < lb1, 1, 0)

inds <- ifelse(ubsind == 0 & lbsind == 0, 0, ifelse(ubsind == 1 & lbsind == 0, 
	1, ifelse(ubsind == 0 & lbsind == 1, 2, 3)))
	
ub2 <- ifelse(tln2$ub >ub1, Inf, tln2$ub)
lb2 <- ifelse(tln2$lb < lb1, -Inf, tln2$lb)
	
tln2 <- data.frame(tln2, inds, ub2, lb2)

limits <- aes(ymax = ub, ymin= lb)


tln3 <- tln2
names2 <- names1
# names2 <- names1[c(1, 2, 3, 5, 7, 4, 8, 6)]
# tln3 <- tln3[-which(tln3$source %in% names2[6 : 8]), ]
# names2 <- names2[1 : 5]
tln3$source <- factor(tln3$source, levels = names2)

s1 <- c("P/V", "As/Se/Br")
s1 <- c("Salt", "Traffic")
# tln3 <- tln3[-which(tln3$source %in% s1), ]
tln3 <- tln3[which(tln3$source %in% s1), ]
tln3$source <- factor(tln3$source, levels = s1)


tln3 <- tln2

size1 <- 18
sizep <- 1.2

col1 <- brewer.pal(5, "Blues")[3:5]
g1 <- ggplot(tln3, aes(x = type, y = est, colour = lag), 
	ylim = c(-0.01, 0.01)) + 
	geom_hline(aes(yintercept = 0), colour = "grey80", 
		linetype = "dashed") +
    geom_pointrange(aes(ymin = lb2, ymax = ub2, colour = lag), 
  	  width = 0.1, position = pd, size = sizep) +
    # scale_colour_hue(drop = FALSE, name="",
                     # # breaks=c("lag0", "lag1", "lag2"),
                     # labels=c("Lag 0", "Lag 1", "Lag 2"),
                     # l=20, c = 150, h = 150) + 
     # scale_colour_brewer(palette="Blues", labels=c("Lag 0", "Lag 1", "Lag 2"), name = "") +   
      scale_color_manual( labels=c("Lag 0", "Lag 1", "Lag 2"), name = "",
                values = col1) + 

    ylab(expression(atop("% increase in mortality per", 
		"10-" * mu * 
		 "g/m"^"3"* " increase in source concentration"))) +
    xlab("") +
    # ggtitle("Mortality effects by\nsource apportionment method") +
    scale_y_continuous(limits=c(lb1, ub1))      +         
    theme_bw() +        
    theme(axis.text.y=element_text(size=size1)) +
  	theme(legend.text=element_text(size=size1), legend.position = "bottom")  +
  	theme(axis.title=element_text(size=size1)) 
    # theme(legend.justification=c(1,0), 
    	# legend.position=c(1,0)) # Position legend in bottom right
g1 <- g1 + theme(axis.text.x=element_text(size = size1))
g1 <- g1 + facet_wrap(~source, ncol = 2) +theme(strip.text.x = element_text(size = size1))
g1 <- g1	+ theme(strip.text.y = element_text(size = size1))

g1







setwd("/Users/jennakrall/Dropbox/SpatialFA/plots")
pdf("hosp_east.pdf", height = 12, width = 10)	
g1
graphics.off() 	
# ###########	
###########	
###########	
###########	
###########	
	
	
	
	


	
	
	
	
	
	
	














###########	
###########	
###########	
###########	
###########		
######
# empirical bayes

#### function to get county info from fips
data(county.fips)
data(state.fips)
#add state
split <- strsplit(as.character(county.fips[, 2]), ",")
states <- sapply(split, function(x) x[1])
county <- sapply(split, function(x) x[2])


#switch to abbreviaton
abbrev <- unique(state.fips[, c("abb", "polyname")])
abbrev[, 2] <- sapply(strsplit(as.character(abbrev[, 2]), ":"), 
	function(x) x[1])
abbrevsall <- as.character(abbrev[match(states, abbrev[, 2]), 1])

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}


county <- sapply(county, simpleCap)
namesC <- paste0(county, ", ", abbrevsall)
matchname <- data.frame(county.fips, namesC)


getcounty <- function(fips) {
	fips <- as.numeric(fips)
	name <- as.character(matchname[which(matchname[, 1] == fips), 3])
	
	if(length(name) == 0) {
		print(fips)
	}
	
	name
	
}


####
#empirical bayes estimates
types <- c("tlnAPCA", "tlnmAPCAsw")
lags <- c("Lag 0", "Lag 1", "Lag 2")
sources <- names1
others <- paste(rep(0, 5))
type <- c("SHARE", "mAPCA")
sources[4] <- "Fireworks" 
dat <- rep(0, 3)
for(i in 1 : 9) {
	for(j in 1 : 2) {
		for(l in 1 : 3) {
			print(c(i, j, l))
			hold <- get(types[j])
			ests <- hold[[1]][[i]]
			
			if(!is.null(ests)) {
				A <- ests[[3]][l]
				citest <- ests[[2]][[l]]
				totest <- ests[[1]][l, 1]
				
				# temp <- eb(citest[, 1], citest[, 2]^2, totest, A)
				temp <- hold[[1]][[i]][[4]][,, l]
				emp <- percinc(temp[, 1])
				# sd <- percinc(sqrt(temp[[2]]))
				sd <- percinc(temp[, 2])
				
				ests <- c(emp, percinc(citest[, 1]))
				ses <- c(sd, percinc(citest[, 2]))
				lb <- ests - 1.96 * ses
				ub <- ests + 1.96 * ses
				
				dat1 <- cbind(ests, lb, ub)
				cities1 <- sapply(rownames(citest), getcounty)
				n <- length(cities1)
				others1 <- cbind(rep(cities1, 2), rep(sources[i], n * 2), 
					rep(lags[l], n * 2), rep(type[j], n * 2),
					rep(c("EB", "MLE"), each = n))
				
				dat <- rbind(dat, dat1)
				others <- rbind(others, others1)
			}
		}
	}
	
}

dat <- dat[-1, ]
others <- others[-1, ]




datU <- data.frame(dat, others)
colnames(datU) <- c("est", "lb", "ub", "city", "source", "lag", "SA", "type")

citfull <- unique(matchname[which(matchname$fips %in% as.numeric(names(unmonlist))), 3])


# #add in extra for missing
# dathold <- datU[which(datU$type == "EB" & 
	# datU$SA == "SHARE" & datU$lag == "Lag 1" & datU$source == "Salt"), ]
# datcit <- as.character(dathold[order(dathold$est, decreasing = T), "city"]	)
# citno <- cities1[-which(cities1 %in% datcit)]
# datcit <- c(datcit, citno)
# citfull <- citfull[match(datcit, citfull)]

# datU$city <- factor(datU$city, levels = citfull)
datU$source <- factor(datU$source, levels = names1)	
datU$SA <- factor(datU$SA, levels = c("SHARE", "mAPCA"))

s1 <- sources
s1 <- s1[-which(s1 %in% c("P/V", "As/Se/Br", "Fireworks"))]
datU2 <- datU
datU2$source <- factor(datU$source, levels = s1)


# 

# cits <- unique(as.character(datU2[which(abs(datU2$est) > 50), "city"]))
# datU2 <- datU2[-which(datU2$city %in% cits), ]


# wh1 <- (which(abs(datU2$est) > ub2))
# if(length(wh1) > 0) {
	# datU2[wh1, "est"] <- Inf * sign(datU[wh1, "est"])
	# }

ggplotColours <- function(n=6, h=c(0, 360) +15, l = 65){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l )
}
col2 <- ggplotColours(n = 3, l = 40)[-1]
size1 <- 14
pd <- position_dodge(.4)
g1fun <- function(dat) {
	lb2 <- -25
	ub2 <- 25
	wh1 <- which(abs(dat$lb) > ub2)	
	if(length(wh1 ) > 0) {
		dat[wh1, "lb"] <- -Inf
	}
	wh1 <- which(abs(dat$ub) > ub2)	
	if(length(wh1 ) > 0) {
		dat[wh1, "ub"] <-  Inf
	}

	g1 <- ggplot(dat, aes(x = city, y = est, colour = SA), 
		ylim = c(-0.01, 0.01)) + 
	geom_hline(aes(yintercept = 0), colour = "grey80", 
		linetype = "dashed") +
    geom_pointrange(aes(ymin = lb, ymax = ub, colour = SA), 
  	  width = 0.1, position = pd) +

          scale_color_manual( labels=c("SHARE", "mAPCA"), name = "",
                values = col2) +            

        ylab(expression(atop("% increase in mortality per", 
		"10-" * mu * 
		 "g/m"^"3"* " increase in source concentration"))) +
    xlab("") +
    scale_y_continuous(limits=c(lb2, ub2)) +     # Set y range
                       # breaks=0:20*4) +                       # Set tick every 4
    theme_bw() +
        theme(axis.text.y=element_text(size=size1)) +
  	theme(legend.text=element_text(size=size1), legend.position = "bottom")  +
  	theme(axis.title=element_text(size=size1)) 

	g1 <- g1 + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size = size1))


	g1 <- g1 + facet_wrap(~source, ncol = 1)  
	g1
}




g1fun(datU3)


setwd("/Users/jennakrall/Dropbox/SpatialFA/plots")

i <- 1
lag2 <- paste0("lag", i-1)
pdf(paste0("hosp_empbayes_east_", lag2, ".pdf"), height = 15, width = 9)
lag <- paste("Lag", i- 1)

datU3 <- datU2[which(datU2$type == "EB" & datU2$lag == lag & 
	datU2$source %in% s1 ), ]
g1fun(datU3)	

graphics.off()





