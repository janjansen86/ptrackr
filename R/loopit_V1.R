
## this allows to track even a great number of particles all the way to the seafloor without producing errors due to limited RAM
## the loop runs in half-day intervals
## if no runtime is defined, the function will loop depending on the depth of the deepest cell and the sinking speed to allow each particle to possibly sink to the seafloor (2*max(h)/speed)
## if a runtime is defined, it is measured in days

#' Loopit
#'
#' Loop trackit_3D to follow particles through different consecutive ROMS-sclices. Looping can also increase performance when using very large number of particles by looping through shorter time steps.
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
#' pts <- create_points_pattern(surface_chl, multi=100)
#' h <- toyROMS$h
#' all_i_u <- toyROMS$i_u
#' all_i_v <- toyROMS$i_v
#' all_i_w <- toyROMS$i_w
#' 
#' track <- loopit(speed = 100, runtime = 50)
#' 
#' #plot(pts)
#' #points(track$pnow, col="red")
#' 
#' #library(rgl)
#' #plot3d(pts, zlim = c(-1500,1))
#' #plot3d(track$pnow, col="red", add=TRUE)

#' 

loopit <- function(speed, runtime, trajectories=FALSE){
  ## create lists to store all particles that settled at the end of each tracking-loop
  lon_list <- list()
  lat_list <- list()
  depth_list <- list()
  ## create lists to store the positions of each particle in each time-step
  if(missing(trajectories)){
    trajectories=FALSE
  } else if(trajectories==TRUE){
    idx_list <- list()
    idx_list_2D <- list()
    id_list <- list()
  }
  if (missing(runtime)){
    runtime <- ceiling(max(h)/speed)                  ## no runtime defined
  } else runtime <- runtime                           ## runtime defined
  curr_vector <- rep(1:4,runtime)
  runtime <- 4*runtime                                ## counting full days
  for(irun in 1:runtime){                             ## loop over different time-slices
    ## assign current-speed/direction to the cells, this should be done differently
    i_u <<- all_i_u[,,,curr_vector[irun]]
    i_v <<- all_i_v[,,,curr_vector[irun]]
    i_w <<- all_i_w[,,,curr_vector[irun]]
    
    ## save an id for each particle to follow its path
    if(trajectories==TRUE){
      id_list[[irun]] <- id_vec
    }
    
    ## run the particle-tracking for all floating particles
    obj <- trackit_3D(speed,days)
    ## store the particles that stopped (settled)
    lon_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 1, obj$stopindex)])
    lat_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 2, obj$stopindex)])
    depth_list[irun] <- list(obj$ptrack[cbind(seq(nrow(obj$ptrack)), 3, obj$stopindex)])
    
    if(trajectories==TRUE){
      ## store the cell-indices of each pts from each time-slice
      idx_list[[irun]] <- obj$indices
      idx_list_2D[[irun]] <- obj$indices_2D
      ## reduce the id_vec to new number of pts
      id_vec <- id_vec[obj$stopindex==0]
    }
    
    ##re-assign coordinates of floating particles to re-run in "trackit"-function
    if (length(unlist(lon_list))!=nrow(check)                                 ## check if all particles are settled
        & !is.null(nrow(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]]))    ## if there's only one particle left it bugs around...
        & length(lon_list)!=runtime){                                         ## if the run stops before all particles are settled, pts should not be overwritten
      ## Mike doesn't like this!! But it works... How to do it better?
      pts <<- matrix(obj$ptrack[obj$stopindex==0,,dim(obj$ptrack)[3]],ncol=3)
    } else break
  }
  ## store stopping-locations of particles
  pend <- cbind(unlist(lon_list), unlist(lat_list), unlist(depth_list))
  
  if(nrow(pts)==1){
    pend <- pend[-nrow(pend),]
  }
  
  ## store the output
  if(trajectories==TRUE){
    list(pts=pts, pend=pend, stopindex=obj$stopindex, ptrack=obj$ptrack, lon_list=lon_list, idx_list=idx_list, idx_list_2D=idx_list_2D, id_list=id_list)
  }else {
    list(pts=pts, pend=pend, stopindex=obj$stopindex, ptrack=obj$ptrack, lon_list=lon_list)
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
