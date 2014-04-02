### File to simulate mortality effects

load("/Users/jennakrall/Dropbox/SpatialFA/data/share_sim_res.RData")
load("/Users/jennakrall/Dropbox/MDL_sourceapp/MDL_project_oct2012/cache/mnvar.RData")


#means and sds	
cms <- rep(mnvar[1, ], 2)
sds <- rep(sqrt(mnvar[2, ]), 2)