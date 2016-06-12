#' Trackit 3D
#'
#' Function to track particles through a ROMS-field.
#'
#' the function needs an input for speed of the sinking particles (w_sink) and for time
#' due to the limitation of RAM available, time is restricted depending on the number of particles
#' (too long runs might give an error because the generated vector is too large)
#' I found half days work great
#'
#' @param pts input points
#' @param kdtree kd tree
#' @param w_sink sinking rate m/days
#' @param time total number of days to run the model
#' @param loop_trackit Default is FALSE, automatically turns TRUE when called from loopit_2D3D
#'
#' @return list(ptrack = ptrack, pnow = pnow, plast = plast, stopindex = stopindex, indices = indices, indices_2D = indices_2D)
#' @export
#' @examples 
#' data(surface_chl)
#' data(toyROMS)
#' pts <- create_points_pattern(surface_chl, multi=100)
#' track <- trackit_3D(pts = pts, romsobject = toyROMS)
#' 
#' ## checking the results
#' plot(pts)
#' points(track$pnow, col = "red")
#' 
#' # library(rgl)
#' # plot3d(pts, zlim = c(-1500,1))
#' # plot3d(track$pnow, col = "red", add = TRUE)
#' 
#' ## better:
#' library(rasterVis)
#' library(rgdal)
#' 
#' ra <- raster(nrow = 50, ncol = 50, ext = extent(surface_chl))
#' r_roms <- rasterize(x = cbind(as.vector(toyROMS$lon_u), as.vector(toyROMS$lat_u)), y = ra, field = as.vector(-toyROMS$h))
#' pr <- projectRaster(r_roms, crs = "+proj=laea +lon_0=137 +lat_0=-66")  #get the right projection (through the centre)
#' 
#' plot3D(pr, adjust = FALSE, zfac = 50)                    # plot bathymetry with 50x exaggerated depth
#' points <- matrix(NA, ncol=3, nrow=dim(track$ptrack)[1])  # get Tracking-points
#' for(i in seq_len(dim(track$ptrack)[1])){
#'   points[i,] <- track$ptrack[i,,track$stopindex[i]] 
#' }
#' pointsxy <- project(as.matrix(points[,1:2]), projection(pr))  #projection on Tracking-points
#' points3d(pointsxy[,1], pointsxy[,2], points[,3]*50)
#' 
#' ptsxy <- project(as.matrix(pts[,1:2]), projection(pr))  #projection on Tracking-points
#' points3d(ptsxy[,1], ptsxy[,2], pts[,3]*50, col = "red")

trackit_3D <- function(pts, romsobject, w_sink=100, time=50, romsparams, loop_trackit=FALSE){

  ## We need an id for each particle to follow individual tracks
  id_vec <- seq_len(nrow(pts))

  ## build a kdtree
  if(loop_trackit=TRUE){
    kdtree <- romsobject$kdtree
    kdxy <- romsobject$kdxy
  }else{
    sknn <- with(romsobject, setup_knn(lon_u, lat_u, hh))             # (lon_roms=lon_u, lat_roms=lat_u, depth_roms=hh)
    kdtree <- sknn$kdtree
    kdxy <- sknn$kdxy
  }

  ## assign current speeds and depth for each ROMS-cell (lat/lon position)
  if(missing(romsparams)){
    i_u <- romsobject$i_u
    i_v <- romsobject$i_v
    i_w <- romsobject$i_w
    h <- romsobject$h
  }else{
    i_u <- romsparams$i_u
    i_v <- romsparams$i_v
    i_w <- romsparams$i_w
    h <- romsparams$h
  }
#   i_u <- romsobject$i_u
#   i_v <- romsobject$i_v
#   i_w <- romsobject$i_w
#   h <- romsobject$h
  
  ## boundaries of the ROMS-area
  # roms_ext <- c(min(romsobject$lon_u), max(romsobject$lon_u), min(romsobject$lat_u), max(romsobject$lat_u))
  
  ## w_sink is m/days, time is days
  w_sink <- -w_sink/(60*60*24)                               ## sinking speed transformation into m/sec
  ntime <- time*24*2                                         ## days transformation into 0.5h-intervalls
  time_step <- 30*60                                         ## half hour time steps
  
  ## empty objects for the loop
  ptrack <- array(0, c(length(as.vector(pts))/3, 3, ntime))  ## create an empty array to store particle-tracks
  stopped <- rep(FALSE, length(as.vector(pts))/3)            ## create a stopping-vector
  stopindex <- rep(0, length(as.vector(pts))/3)              ## a vector to store indices of when particles stopped
  indices <- vector("list", ntime)                           ## a list of indices to store which 3D-cell a particle is in
  indices_2D <- vector("list", ntime)                        ## a list of indices to store which 2D-cell a particle is in                       
  
  pnow <- plast <- pts                ## copies of the starting points for updating in the loop
  
  for (itime in seq_len(ntime)) {
    ## index 1st nearest neighbour of trace points to grid points
    dmap <- kdtree$query(plast, k = 1, eps = 0)           ## one kdtree
    ## and to 2D space
    two_dim_pos <- kdxy$query(plast[,1:2], k = 1, eps = 0)

    ## store indices for tracing particle positions
    indices[[itime]] <- dmap$nn.idx
    indices_2D[[itime]] <- two_dim_pos$nn.idx

    ## extract component values from the vars
    thisu <- i_u[dmap$nn.idx]                             ## u-component of ROMS
    thisv <- i_v[dmap$nn.idx]                             ## v-component of ROMS
    thisw <- i_w[dmap$nn.idx]                             ## w-component of ROMS

    ## update this time step longitude, latitude, depth
    pnow[,1] <- plast[,1] + (thisu * time_step) / (1.852 * 60 * 1000 * cos(pnow[,2] * pi/180))
    pnow[,2] <- plast[,2] + (thisv * time_step) / (1.852 * 60 * 1000)
    pnow[,3] <- pmin(0, plast[,3])  + ((thisw + w_sink)* time_step )

#     stopped <- (pnow[,1] < roms_ext[1] | pnow[,1] > roms_ext[2] |
#                   pnow[,2] < roms_ext[3] | pnow[,2] > roms_ext[4]
#                 )
#     
    ##########---- only in trackit_3D:
    
    ## stopping conditions (hit the bottom)
    stopped <- pnow[,3] <= -h[two_dim_pos$nn.idx]
    stopindex[stopindex == 0 & stopped] <- itime
    
    ##########----
    
    ## assign stopping location of points to ptrack
    ptrack[,,itime] <- pnow
    plast <- pnow
    print(itime)
    if (all(stopped)) {
      message("exiting, all stopped")
      break;
    }
  }
  ptrack <- ptrack[,,seq(itime),drop=FALSE]
  list(ptrack = ptrack, pnow = pnow, plast = plast, stopindex = stopindex, indices = indices, indices_2D = indices_2D)
}



####################################
## for the settlingmodel?

## This function is called in the function loopit
## the function needs an input for speed of the sinking particles (vertical_movement), the time of the run and the resolution of the run (time_step)
## due to the limitation of RAM available, time is restricted depending on the number of particles
## (too long runs might give an error because the generated vector is too large)

## u, v and w components need to be defined (i_u, i_v, i_w)

## give some room to define stopping conditions

# trackit <- function(vertical_movement,time,time_step){ ## vertical_movement is m/s, time is s, steps is in s
#   ptrack <- array(0, c(length(as.vector(pts))/3, 3, time))  ## create an empty array to store particle-tracks
#   stopped <- rep(FALSE, length(as.vector(pts))/3)            ## create a stopping-vector
#   stopindex <- rep(0, length(as.vector(pts))/3)              ## a vector to store indices of when particles stopped
#   plast <- pts                                               ## copies of the starting points for updating in the loop
#   pnow <- plast
#   #plot(pts, pch = ".")                                      ## to plot particles after each timestep
#   for (itime in seq_len(time)) {
#     ## index 1st nearest neighbour of trace points to grid points
#     dmap <- kdtree$query(plast, k = 1, eps = 0)           ## one kdtree from script "..._readit.R"
#     ## extract component values from the vars
#     thisu <- i_u[dmap$nn.idx]                             ## u-component of ROMS
#     thisv <- i_v[dmap$nn.idx]                             ## v-component of ROMS
#     thisw <- i_w[dmap$nn.idx]                             ## w-component of ROMS
#
#     ## update this time step longitude, latitude, depth
#     pnow[,1] <- plast[,1] + (thisu * time_step) / (1.852 * 60 * 1000 * cos(pnow[,2] * pi/180))
#     pnow[,2] <- plast[,2] + (thisv * time_step) / (1.852 * 60 * 1000)
#     pnow[,3] <- pmin(0, plast[,3])  + ((thisw + vertical_movement)* time_step )
#
#     ## hit the bottom
#     stopped <- pnow[,3] <= -h[kdxy$query(pnow[,1:2], k = 1, eps = 0)$nn.idx]
#     stopindex[stopindex == 0 & stopped] <- itime
#     ptrack[,,itime] <- pnow
#     plast <- pnow
#     print(itime)
#     if (all(stopped)) {
#       message("exiting, all stopped")
#       break;
#     }
#   }
#   ptrack <- ptrack[,,seq(itime)]
#   list(ptrack = ptrack, pnow = pnow, plast = plast, stopindex = stopindex)
# }
#
