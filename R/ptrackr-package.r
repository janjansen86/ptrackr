#' ptrackr: a package to track individual particles through 2D- and 3D-space
#' 
#' The ptrackr package intends to provide the basic functions needed to store the positions and track the movement of individual particles through time and space. 
#' The package is written in the context of movement of particles in a velocity field taken from a Regional Oceanographic Modelling System (ROMS).
#' The main function to run particle-tracking is called "loopit_2D3D" (examples are provided below), with the other functions either being called from the main function or serving to setup the data for the modelling.
#' More details can be found in the individual desription files.
#'
#' @section overview:
#' The package is setup into 3 main functions. The functions that do the actual tracking of the particles are trackit_2D and trackit_3D for 2D and 3D respectively. 
#' These two functions work with user-specified time-steps (e.g. 30min) at which the current position of each particle and their movement during the next time-steps gets evaluated. 
#' The function loopit_2D3D loops one of the two trackit-functions over given time-intervalls (e.g. 6h, so in this example one loop of 6h with 30min time-steps would result in 12 rounds per loop). 
#' The advantage is that at the end of each loop those particles that have met a stopping criteria can be excluded from further loops, thereby greatly increasing the speed of the model-run.
#' 
#' @section basic 2D-example:
#' data(surface_chl); data(toyROMS)
#'  
#' pts_seeded <- create_points_pattern(surface_chl, multi=100)
#' 
#' run <- loopit_2D3D(pts_seeded = pts_seeded, roms_slices = 4, romsobject = toyROMS, speed = 100, runtime = 50, sedimentation = TRUE)
#' 
#' plot(pts_seeded); points(run$pend, col="red", cex=0.6); points(run$pts , col="blue", cex=0.6)
#' 
#' legend("bottomleft", pch=1, pt.cex=c(1,0.6,0.6), col=c("black","red","blue"), legend=c("seed-location","end-location","still moving"))
#' 
#' @section Authors:
#' Jan Jansen & Michael D. Sumner
#' 
#' Maintainer: Jan Jansen - Jan.Jansen(AT)utas.edu.au
#' 
#' DOI: 10.5281/zenodo.583336
#' 
#' URL: www.github.com/janjansen86/ptrackr
#'  
#' @docType package
#' @name ptrackr
NULL
