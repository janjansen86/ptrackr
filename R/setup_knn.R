## setup k-nearest-neighbour upfront for use in the particle tracking loops
## lookup to which ROMS-cell a particle belongs

#' Setup knn lookup to which ROMS-cell a particle belongs
#'
#' Setup k-nearest-neighbour upfront for use in the particle tracking loops.
#'
#' Details
#' @param lon_roms longitude values of ROMS grid
#' @param lat_roms latitude values of ROMS grid
#' @param depth_roms depth values of ROMS grid
#'
#' @return list, kdtree, kdxy
#' @export
setup_knn <- function(lon_roms=lon_u, lat_roms=lat_u, depth_roms=hh){
  library(nabor)
  ## another function: setup knn
  allxyz <- cbind(as.vector(lon_roms), as.vector(lat_roms), as.vector(depth_roms))
  ## needs to be a separate lookup for detecting when we hit the bottom
  xytop <- seq_len(prod(dim(lon_roms)))
  ## this is our 3D-search tree engine, built once upfront
  kdtree <- WKNND(allxyz)
  ## we store a second 2D-tree for the halting condition
  kdxy <-  WKNND(allxyz[xytop, 1:2 ]) ## knn(allxyz[xytop, ], xyzt[,1:2,itime - 1], k = 1)
  list(kdtree=kdtree, kdxy=kdxy)
}
