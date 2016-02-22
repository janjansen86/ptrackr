#' Loopit_Trackit 2D
#'
#' Function called by loopit to track particles through a 2D-ROMS-field.
#'
#' @param pts input points
#' @param kdtree kd tree
#' @param w_sink sinking rate m/days
#' @param time total number of days to run the model
#'
#' @return list(ptrack = ptrack, pnow = pnow, plast = plast, stopindex = stopindex, indices = indices, indices_2D = indices_2D)
#' @export

loopit_trackit_2D <- function(pts, romsobject, w_sink=100, time=50){
  
  roms_ext <- c(min(romsobject$lon_u), max(romsobject$lon_u), min(romsobject$lat_u), max(romsobject$lat_u))
  
  ## We need an id for each particle to follow individual tracks
  id_vec <- seq_len(nrow(pts))

  sknn <- with(romsobject, setup_knn(lon_u, lat_u, hh))             # (lon_roms=lon_u, lat_roms=lat_u, depth_roms=hh)
  kdtree <- sknn$kdtree
  kdxy <- sknn$kdxy
  
  ## its essentially the same as trackit_2D except these lines are commented here:
    #i_u <- romsobject$i_u
  #i_v <- romsobject$i_v
  #i_w <- romsobject$i_w
  #h <- romsobject$h
  
  ## w_sink is m/days, time is days
  w_sink <- -w_sink/(60*60*24)                               ## sinking speed transformation
  ntime <- time*24*2                                         ## days transformation
  time_step <- 30*60                                         ## half hour time steps
  ptrack <- array(0, c(length(as.vector(pts))/3, 3, ntime))  ## create an empty array to store particle-tracks
  stopped <- rep(FALSE, length(as.vector(pts))/3)            ## create a stopping-vector
  stopindex <- rep(0, length(as.vector(pts))/3)              ## a vector to store indices of when particles stopped
  ## copies of the starting points for updating in the loop
  pnow <- plast <- pts
  indices <- vector("list", ntime)                           ## a list of indices to store which 3D-cell a particle is in
  indices_2D <- vector("list", ntime)                        ## a list of indices to store which 2D-cell a particle is in

  for (itime in seq_len(ntime)) {
    ## index 1st nearest neighbour of trace points to grid points
    dmap <- kdtree$query(plast, k = 1, eps = 0)           ## one kdtree
    ## and to 2D space
    two_dim_pos <- kdxy$query(pnow[,1:2], k = 1, eps = 0)

    ## store indices for tracing particle positions
    indices[[itime]] <- dmap$nn.idx
    indices_2D[[itime]] <- two_dim_pos$nn.idx

    ## extract component values from the vars
    thisu <- i_u[dmap$nn.idx]                             ## u-component of ROMS
    thisv <- i_v[dmap$nn.idx]                             ## v-component of ROMS
    #thisw <- i_w[dmap$nn.idx]                             ## w-component of ROMS
    thish <- h[dmap$nn.idx]                               ## depth of ROMS-cell
    
    ## update this time step longitude, latitude, depth
    pnow[,1] <- plast[,1] + (thisu * time_step) / (1.852 * 60 * 1000 * cos(pnow[,2] * pi/180))
    pnow[,2] <- plast[,2] + (thisv * time_step) / (1.852 * 60 * 1000)
    #pnow[,3] <- pmin(0, plast[,3])  + ((thisw + w_sink)* time_step )
    
    ################################
    ## k returns the cell-index of each point   (find the nearest grid-point from lon_u/lat_u and return its index)
    tdp_idx <- two_dim_pos$nn.idx
    
    ## make sure particles are not travelling upwards
    #if (depth of pnow) < (depth of plast) then assign plast to pnow with no displacement
    uphill <- h[tdp_idx] > thish +30
    
    pnow[uphill==TRUE,] <- plast[uphill==TRUE,]
    
    ## how many particles are there for each cell-index    which(tabulate(k)!=0) selects the same cells as table(k), but keeps their index
    all_dens <- tabulate(tdp_idx)
    
    #score <- as.vector(all_dens)
    ##########################################
    ## speed and density for each cell:
    l <- unique(tdp_idx)  #indices of each "used" cell
    #speed in active cells, needs to be in squared for equation
    cell_chars <- data.frame(cbind(l,i_u[l] ^ 2 + i_v[l] ^ 2))
    #point-density in active cells
    cell_chars[,3] <- all_dens[l]
    #get u_div from observed and critical velocity 
    U_div <- 1-(cell_chars[,2] / Ucsq)
    ##no erosion:    
    U_div[U_div<0] <-0
    ## calculate number of points to settle for each cell (equation from McCave & Swift)
    cell_chars[,4] <- (testFunct(U_div, cell_chars[,3]))
    
    colnames(cell_chars) <- c("cell_index","velocity","n_pts_in_cell","n_pts_to_drop")
    
    ## stopping when outside the limits of the area 
    stopped <- (pnow[,1] >= roms_ext[1] | pnow[,1] <= roms_ext[2] | pnow[,2] >= roms_ext[3] | pnow[,2] <= roms_ext[4])
    
    ## forced settling out of the suspension:                      
    ## create an object to store points to be dropped
    #drop_pts <- rep(F,length(pnow)/2)                     
    
    ## qd: the quick and dirty solution
    point_chars <- data.frame(cbind(tdp_idx, cell_chars[match(tdp_idx, cell_chars[,1]),3:4]))
    point_chars[,4] <- point_chars[,3] / point_chars[,2]
    t_f <- runif(nrow(point_chars)) <= point_chars[,4]
    stopindex[(stopindex == 0 & stopped) | (stopindex == 0 & t_f)] <- itime
    ########################

    ## assign stopping location of points to ptrack 
    ptrack[,,itime] <- pnow
    plast <- pnow
    print(itime)
    if (all(stopindex!=0)) {
      message("exiting, all stopped")
      break 
    }
  }
  ptrack <- ptrack[,,seq(itime)]
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
