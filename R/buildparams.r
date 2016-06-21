#' buildparams
#'
#' Function to set up parameters for the sedimentation equation
#'
#' @param speed sinking speed
#' @param r particle radius
#'
#' @return list(p0 = p0, p1 = p1, cosO = cosO, g = g, K = K, E = E, r = r, Wd = Wd, Ucsq = Ucsq, testFunct = testFunct)
#' @export


buildparams <- function(speed, 
  ## from Jenkins & Bombosch (1995)
  p0 =1030,             #kg/m^3 seawater density
  p1 =1100,             #kg/m^3 Diatom density (so far a quick-look-up-average density from Ierland & Peperzak (1984))
  cosO =1,              #its 1 for 90degrees
  g =9.81,              #accelaration due to gravity
  K =0.0025,            #drag coefficient
  E =1,                #aspect ration of settling flocks (spherical = 1 ??)
  r =0.00016,           #particle-radius
  Ucsq =-(0.05*(p0-p1)*g*2*(1.5*E)^(1/3)*r)/(p0*K),
  Wd =speed,#/24/3600,
  testFunct =function(U_div,dens) 1800*-(p1*(dens)*Wd*cos(90)*(U_div)*(U_div))/p0){
  list(p0 = p0, p1 = p1, cosO = cosO, g = g, K = K, E = E, r = r, Wd = Wd, Ucsq = Ucsq, testFunct = testFunct)
}

