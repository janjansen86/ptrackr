###################################################################################################################################
## ptrackR
## J.Jansen
## 2016
###################################################################################################################################

## data() load these

# only one slice, in this case an average over the summer currents
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/ToyROMS_data.Rdata")
list2env(ToyROMS_data ,.GlobalEnv)

#multiple consecutive current-slices in a row to track particles in changing currents
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/toyROMS58.Rdata")
list2env(toyROMS58 ,.GlobalEnv)
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/toyROMS59.Rdata")
list2env(toyROMS59 ,.GlobalEnv)
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/toyROMS60.Rdata")
list2env(toyROMS60 ,.GlobalEnv)
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/toyROMS61.Rdata")
list2env(toyROMS61 ,.GlobalEnv)

# > > this would be data(toyroms)


# load surface chlorophyll raster
load("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/data/surface_chl.Rdata")
jChlA <- crop(surface_chl,extent(136.6,138.9,-66.31,-65.75))

# > > this would be data(surface_chl)

## transform values into particle distribution
source("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/R/create_points_pattern.R")
cpp <- create_points_pattern(jChlA,10)   # (a_chl_raster,a_multiplicator=1,a_mixed_layer_depth=0)
                                         # requires: library(geostatsp)
pts <- cpp$pts

## We need an id for each particle to follow individual tracks
id_vec <- seq_len(nrow(pts))

## setup lookup table
source("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/R/setup_knn.R")
sknn <- setup_knn()             # (lon_roms=lon_u, lat_roms=lat_u, depth_roms=hh)
kdtree <- sknn$kdtree
kdxy <- sknn$kdxy


## either static ROMS,
curr_vector <- rep(1,400)   #for 1 slice over 100days
all_i_u <- list(i_u)        #for multiple slices choose multiple i_u
all_i_v <- list(i_v)        #for multiple slices choose multiple i_v
all_i_w <- list(i_w)
## or multiple ROMS
curr_vector <- rep(1:4,400) #for 4 slices in 100days
all_i_u <- list(i_u1,i_u2,i_u3,i_u4)
all_i_v <- list(i_v1,i_v2,i_v3,i_v4)
all_i_w <- list(i_w1,i_w2,i_w3,i_w4)


source("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/R/R_package_ptrack_trackit_3D.R")
source("C:/Users/jjansen/Desktop/PhD/ptrackR_package/ptrackr/R/R_package_ptrack_loopit_V1.R")

###########################################################################################################################
## time of each individual run in trackit-function, this should be left 0.5 days for lots of particles
days <- 0.5

## choose a number of different sinking speeds (in m/day) to be tested
all_speeds <- 200 #seq(50,200,50)
#with the max(depth) of the ROMS-area and the sinking speed, calculate the time the model should run
all_days <- ceiling(max(h)/all_speeds)

## this is actually only important to save output from different runs
library(stringr)
all_speeds2 <- str_pad(all_speeds, 3, pad = "0")        #generate a 3digit-vector of speed (e.g. 050m instead of 50m)
all_stuff <- rbind(all_speeds2,as.character(all_days))

pts_seeded <- pts
pts_settled <- list()
pts_floating <- list()

## tracking can be either run with scoring particle-locations at each time-step (slow, but neccessary for some applications)
#trajectories = TRUE
## or without scoring particle-locations at each time-step (faster)
trajectories = FALSE

## for-loop to run multiple sinking speed
for(inumber in 1:ncol(all_stuff)){
  message(paste0(inumber,".run"))
  xspeed <- as.numeric(all_stuff[1,inumber])

  ## re-assign original points because pts gets overwritten in the loop
  pts <- pts_seeded

  ## copy to check stopping conditions of the loop
  check <- pts

  ## run the tracking for a given sinking speed
  run <- loopit(xspeed,trajectories=trajectories)

  ## use the output
  xpend <- run$pend
  pts_settled[[inumber]] <- xpend
  pfloat <- pts[run$stopindex==0,]
  pts_floating[[inumber]] <- pfloat

  if(trajectories == TRUE){
    idx_list <- run$idx_list
    idx_list_2D <- run$idx_list_2D
    id_list <- run$id_list
  }

  ## for saving output of each run individually
  #xspeed2 <- all_stuff[1,inumber]
  #xdays <- all_stuff[2,inumber]
  #save(xpend, file = paste0(res,"_",xspeed2,"_",xdays,"_",mld,"_pend.Rdata"))
  #save(pfloat, file = paste0(res,"_",xspeed2,"_",xdays,"_",mld,"_pfloat.Rdata"))
  #save(idx_list, file = paste0(res,"_",trackspeed,"_",n_days,"_",mld,"_idx.Rdata"))
  #save(id_list, file = paste0(res,"_",trackspeed,"_",n_days,"_",mld,"_id.Rdata"))
}

###########################################################################################################################

###################################################################
## not part of the package, just for checking the output
## look at the results
plot(pts_settled[[1]][,1:2])
points(pts_seeded[,1:2],col="red",cex=0.5)


## plot the trajectories (this doesn't work automated yet)
points(lon_u[idx_list_2D[[1]]],lat_u[idx_list_2D[[1]]],col="green")

library(plyr)
mat_list <- list()
for(islices in 1:8){
  mat_list[[islices]] <- matrix(unlist(idx_list_2D[[islices]]),ncol=24)
}
mat_list[[9]] <- matrix(unlist(idx_list_2D[[9]]),nrow=nrow(idx_list_2D[[9]][[1]]))
#mat_list[[9]] <- cbind(mat_list[[9]],NA)
testmatrix <- do.call(rbind, mat_list)
testid <- unlist(id_list)
## no split the matrix by the vector, resulting ids will not be in order!
flux_list <- split(testmatrix,testid)
points(lon_u[flux_list[[1]]],lat_u[flux_list[[1]]],col="green")
points(lon_u[flux_list[[15]]],lat_u[flux_list[[15]]],col="green")
points(lon_u[flux_list[[5]]],lat_u[flux_list[[5]]],col="green")
