#' Generate points proportionally to surface chlorophyll-a values.
#'
#' Transform surface-chlorophyll values into a randomly distributed relative number of points on a plane depending on the value of each cell.
#' A multiplicator value is used to control the total number of points. A mixed layer depth is defined in the third column (particles are at the surface at zero, below surface at negative values)
#' @param raster input chlorophyll RasterLayer
#' @param multi  multiplicator to control number of points
#' @param mld mixed layer depth
#'
#' @return list with matrix of longitude, latitude, mixed layer depth
#' 
#' @examples  
#' data(surface_chl)
#' pts <- create_points_pattern(surface_chl, multi=1000)
#' plot(surface_chl)
#' points(pts,cex=0.1)
#' 
#' @export
#' 
create_points_pattern <- function(raster, multi = 1, mld = 0){
  ## requires:
  #require(geostatsp)
  #require(spatstat)
  chl_image <- as.im.RasterLayer(raster*multi)  #converts raster into an image with a multiplicator
  ptspattern <- rpoispp(chl_image)      #random poisson point process transform values into probability-density of points
  lon <- ptspattern$x                     #get locations of all points
  lat <- ptspattern$y
  ## store the initial point-locations and save them for comparison later
  cbind(lon,lat,as.vector(as.numeric(mld)))

}
