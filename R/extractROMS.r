#' Extract single ROMS time step.
#'
#' See \link{toyROMS} for an example.
#' @param roms ROMS object
#' @param time index into 4th dim of i_u,i_v,i_w in ROMS object
#'
#' @return list single time slice ROMS
#' @export
#'
extractROMS <- function(roms, time = 1) {
  if (!length(time ) == 1) stop("time must be of length 1")
  if (time < 1 | time > dim(toyROMS$i_u)[4]) stop("time out of bounds")
  obj <- toyROMS  ## in-built data set
  obj$i_u <-obj$i_u[,,,time]
  obj$i_v <-obj$i_v[,,,time]
  obj$i_w <-obj$i_w[,,,time]

  obj
}