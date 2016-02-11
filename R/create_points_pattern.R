

#' Generate points proportionally to surface chlorophyll-a values.
#'
#' Transform surface-chlorophyll values into a randomly distributed relative number of points depending on the value of each cell
#' needs a chlorophyll raster and a multiplicator value to determine the number of points
#' mixed layer depth can be used to assign the third dimension (particles start from the surface at 0, from below surface at negative values)
#' @param a_chl_raster input chlorophyll RasterLayer
#' @param a_multiplicator  scale
#' @param a_mixed_layer_depth mixed layer depth
#'
#' @return list with matrix of longitude, latitude, mixed layer depth
#' @export
#'
create_points_pattern <- function(a_chl_raster, a_multiplicator = 1, a_mixed_layer_depth = 0){
  ## requires:
  library(geostatsp)
  library(spatstat)
  chl_image <- asImRaster(a_chl_raster*a_multiplicator)  #converts raster into an image with a multiplicator
  ptspattern <- rpoispp(chl_image)      #random poisson point process transform values into probability-density of points
  lon <- ptspattern$x                     #get locations of all points
  lat <- ptspattern$y
  ## store the initial point-locations and save them for comparison later
  cbind(lon,lat,as.vector(as.numeric(a_mixed_layer_depth)))

}
