#' Loopit 2D/3D
#' 
#' WORKING FINE FOR BOTH 2D AND 3D
#'
#' Function to run the functions loopit_trackit_2D/loopit_trackit_3D to follow particles through different consecutive ROMS-sclices. Looping can also increase performance when using very large number of particles by looping through shorter time steps.
#' Loops are set to run in half day intervals. If no runtime is defined, the function will loop depending on the depth of the deepest cell and the sinking speed to allow each particle to possibly sink to the seafloor (2*max(h)/speed)
#'
#' @param speed (w_sink) sinking rate m/days
#' @param runtime (time) total number fo days to run the model
#' @param trajectories TRUE/FALSE statement to define whether to store particle trajectories (model runs much faster without storing trajectories). Default is FALSE.
#' 
#' @return list(pts=pts, pend=pend, stopindex=obj$stopindex, ptrack=obj$ptrack, lon_list=lon_list, idx_list=idx_list, idx_list_2D=idx_list_2D, id_list=id_list)
#' @export
#' @examples 
#' data(surface_chl)
#' data(toyROMS)
#' 
#' ########## Case 1:
#' pts_seeded <- create_points_pattern(surface_chl, multi=100)
#' run <- loopit_2D3D(pts_seeded = pts_seeded, romsobject = toyROMS, speed = 100, runtime = 50, domain = "3D", trajectories = TRUE)
#' 
#' ## testing the output (for domain="3D")
#' library(rasterVis)
#' library(rgdal)
#' library(rgl)
#' ra <- raster(nrow=50,ncol=50,ext=extent(surface_chl))
#' r_roms <- rasterize(x = cbind(as.vector(toyROMS$lon_u), as.vector(toyROMS$lat_u)), y= ra, field = as.vector(-toyROMS$h))
#' pr <- projectRaster(r_roms, crs = "+proj=laea +lon_0=137 +lat_0=-66")  #get the right projection (through the centre)
#' 
#' plot3D(pr, adjust = FALSE, zfac = 50)                    # plot bathymetry with 50x exaggerated depth
#' pointsxy <- project(as.matrix(run$pend[,1:2]), projection(pr))  #projection on Tracking-points
#' points3d(pointsxy[,1],pointsxy[,2],run$pend[,3]*50)#,xlim=xlim,ylim=ylim)
#' 
#' 
#' ########## Case 2:
#' pts_seeded <- create_points_pattern(surface_chl, multi=100)
#' run <- loopit_2D3D(pts_seeded = pts_seeded, romsobject = toyROMS, speed = 100, runtime = 50)
#' 
#' plot(pts_seeded)
#' points(run$pend, col="red", cex=0.6)
#' points(run$pts , col="blue", cex=0.6)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ########## Case 2
#' ## work with trajectories to get a flux (added presences/absences) of particles
#' run <- loopit_2D3D(pts_seeded = pts_seeded, romsobject = toyROMS, speed = 100, runtime = 50, trajectories = TRUE)
#'  
#' ## THIS WAS WORKING BEFORE, only run once the functions are fixed:
#' ## this should be abother function to handle the output
#' mat_list <- list()
#' for(islices in 1:length(run$idx_list_2D)){
#'   mat_list[[islices]] <- matrix(unlist(run$idx_list_2D[[islices]]),ncol=24)
#' }
#' testmatrix <- do.call(rbind, mat_list)
#' testid <- unlist(run$id_list)
#' flux_list <- split(testmatrix,testid)
#' for(k in 1:nrow(pts_seeded)){
#'   ## cells visited by a particle ("presence-only")
#'   flux_list[[k]] <- unique(flux_list[[k]])
#'   ## drop first and last value (input and setting cell)
#'   flux_list[[k]] <- flux_list[[k]][-c(1,length(flux_list[[k]]))]
#' } 
#' flux <- as.vector(unlist(flux_list))
#' 


loopit_2D3D <- function(pts_seeded, romsobject, speed, runtime = 10, domain = "2D", 
                        looping_time = 0.25, roms_slices = 1, trajectories = FALSE){
  params <- NULL
  if(domain == "2D"){
     buildparams(speed)  # loopit_trackit_2D only needs testFunct
#     params <- list(
#     ## from Jenkins & Bombosch (1995)
#     p0 =1030,             #kg/m^3 seawater density
#     p1 =1100,             #kg/m^3 Diatom density (so far a quick-look-up-average density from Ierland & Peperzak (1984))
#     cosO =1,              #its 1 for 90degrees
#     g =9.81,              #accelaration due to gravity
#     K =0.0025,            #drag coefficient
#     E =1,                 #aspect ration of settling flocks (spherical = 1 ??)
#     r =0.00016,           #particle-radius
#     Wd =speed/24/3600,
#     Ucsq =-(0.05*(p0-p1)*g*2*(1.5*E)^(1/3)*r)/(p0*K),
#     testFunct =function(U_div,dens) 1800*-(p1*(dens)*Wd*cos(90)*(U_div)*(U_div))/p0
#     )
#     
    }
  romsparams <- list()
  romsparams$h <- romsobject$h
  
# h <<- romsobject$h
  all_i_u <- romsobject$i_u
  all_i_v <- romsobject$i_v
  all_i_w <- romsobject$i_w
  pts <- pts_seeded
  loop_length <- looping_time*24*2
  
  ## create lists to store all particles that settled at the end of each tracking-loop
  lon_list <- list()
  lat_list <- list()
  if(domain == "3D") 
    depth_list <- list()
  
  ## create lists to store the positions of each particle in each time-step
  if(missing(trajectories)){
    trajectories <- FALSE
  } else if(trajectories == TRUE){
    idx_list <- list()
    idx_list_2D <- list()
    id_list <- list()
    id_vec <- seq_len(nrow(pts_seeded))
  }
  if(missing(runtime)){
    runtime <- ceiling(max(h)/speed)                  ## no runtime defined
  } else runtime <- runtime                           ## runtime defined
  curr_vector <- rep(1:roms_slices,runtime)
  runtime <- roms_slices*runtime                                ## counting full days
  
  ## loop over different time-slices
  for(irun in 1:runtime){                             
    
    ## assign current-speed/direction to the cells, this should be done differently
    romsparams$i_u <- all_i_u[,,,curr_vector[irun]]
    romsparams$i_v <- all_i_v[,,,curr_vector[irun]]
    romsparams$i_w <- all_i_w[,,,curr_vector[irun]]
    
#     i_u <<- all_i_u[,,,curr_vector[irun]]
#     i_v <<- all_i_v[,,,curr_vector[irun]]
#     i_w <<- all_i_w[,,,curr_vector[irun]]
    
    ## save an id for each particle to follow its path
    if(trajectories == TRUE) id_list[[irun]] <- id_vec
    
    ## run the particle-tracking for all floating particles
    if(domain == "3D"){
#       obj <- loopit_trackit_3D(pts = pts, romsobject = romsobject, 
#                                w_sink = speed, time = looping_time, parameters = params)
      obj <- trackit_3D(pts = pts, romsobject = romsobject, w_sink = speed, time = looping_time,
                        romsparams=romsparams)
      
    }else{
#       obj <- loopit_trackit_2D(pts = pts, romsobject = romsobject, w_sink = speed, time = looping_time)
      obj <- trackit_2D(pts = pts, romsobject = romsobject, w_sink = speed, time = looping_time,
                        romsparams=romsparams)
      
    }
      
    ## store the particles that stopped (settled)
    lon_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 1, obj$stopindex)])
    lat_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 2, obj$stopindex)])
    ## only for 3D
    if(domain == "3D")
      depth_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 3, obj$stopindex)])
     
    if(trajectories){
      ## store the cell-indices of each pts from each time-slice
      idx_list[[irun]] <- obj$indices
      idx_list_2D[[irun]] <- obj$indices_2D
      ## reduce the id_vec to new number of pts
      id_vec <- id_vec[obj$stopindex==0]
      ## create vector to check if list has some NULL in following if-statement
      NULL_test <- as.character(idx_list[[irun]])
      if(any(NULL_test == "NULL") == TRUE){
        fill_up_seq <- which(NULL_test == "NULL")
        for(ifill in fill_up_seq){
          idx_list[[irun]][[ifill]] <- matrix(NA, nrow = nrow(idx_list[[irun]][[1]]))
          idx_list_2D[[irun]][[ifill]] <- matrix(NA, nrow = nrow(idx_list_2D[[irun]][[1]]))
        }
      }
    } 
      
    ##re-assign coordinates of floating particles to re-run in "trackit"-function
    if (length(unlist(lon_list))!=nrow(pts_seeded)                                 ## check if all particles are settled
        & !is.null(nrow(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]]))    ## if there's only one particle left it bugs around...
        & length(lon_list)!=runtime){                                         ## if the run stops before all particles are settled, pts should not be overwritten
      if(domain == "3D"){
        ## Mike doesn't like this!! But it works... How to do it better?
        pts <- matrix(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]],ncol=3)
      } else pts <- matrix(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]],ncol=3)
    } else break
  }
  
  ## store stopping-locations of particles
  if(domain == "3D"){
    pend <- cbind(unlist(lon_list), unlist(lat_list), unlist(depth_list))
  }else pend <- cbind(unlist(lon_list), unlist(lat_list))
  
  if(nrow(pts)==1) pend <- pend[-nrow(pend),]
  
  message(paste0((dim(pts_seeded)[1]-dim(pend)[1])," particle(s) still floating"))
  
  ## store the output
  if(trajectories==TRUE){
    list(pts=pts, pend=pend, pnow=obj$pnow, stopindex=obj$stopindex, ptrack=obj$ptrack, lon_list=lon_list, idx_list=idx_list, idx_list_2D=idx_list_2D, id_list=id_list)
  }else {
    list(pts=pts, pend=pend, pnow=obj$pnow, stopindex=obj$stopindex, ptrack=obj$ptrack, lon_list=lon_list)
  }
}






#############################
## for the settlingmodel??

## this function allows to track a greater number of particles when an error in trackit occurs due to limited RAM
## this is done by looping a number of short trackings, each using the output of the previous short tracking as their input
## "time" needs to be defined as a fraction of the total "runtime" of the tracking

# loopit <- function(vertical_movement=0.001,time=86400,time_step=3600,in_days=F,runtime=4320000){
#   if(in_days==TRUE){
#     vertical_movement <- vertical_movement*60*60*24
#     time <- time*60*60*24
#     time_step <- time_step*60*60*24
#     runtime <- runtime*60*60*24
#   }
#   ## create lists to store all particles that settled at the end of each tracking-loop
#   long_list <- list()
#   lat_list <- list()
#   vertical_list <- list()
#   #if (missing(runtime)){
#   #      runtime <- ceiling(max(h)/vertical_movement)                                          ## no runtime defined
#   #} else runtime <- runtime                                                      ## runtime defined
#
#   for(irun in 1:(runtime/time)){
#     message(paste0(irun,".loop"))
#       ## run the particle-tracking for all floating particles
#       obj <- trackit(vertical_movement,time,time_step)
#       ## store the particles that stopped (settled)
#       long_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 1, obj$stopindex)])
#       lat_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 2, obj$stopindex)])
#       vertical_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 3, obj$stopindex)])
#       ##re-assign coordinates of floating particles to re-run in "trackit"-function
#       if (length(unlist(long_list))!=nrow(check)                                 ## check if all particles are settled
#           & !is.null(nrow(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]]))    ## if there's only one particle left it bugs around...
#           & length(long_list)!=runtime){                                         ## if the run stops before all particles are settled, pts should not be overwritten
#            ## Mike doesn't like this!! But it works... How to do it better?
#            pts <<- matrix(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]],ncol=3)
#       } else break
#   }
#   ## store stopping-locations of particles
#   pend <- cbind(unlist(long_list),unlist(lat_list),unlist(vertical_list))
#   ## store remaining unsettled particles
#   list(pts=pts,pend=pend,stopindex=obj$stopindex,ptrack=obj$ptrack,long_list=long_list)
# }
