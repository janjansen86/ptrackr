#' buildparams
#'
#' Function to set up parameters for the sedimentation equation used in the trackit_2D function
#'
#' @param speed sinking speed
#' @param time_step_in_s the time-step in seconds
#' @param r particle radius
#'
#' @return list(p0 = p0, p1 = p1, cosO = cosO, g = g, K = K, E = E, r = r, Wd = Wd, Ucsq = Ucsq, SedFunct = SedFunct)
#' @export


buildparams <- function(speed,
                        time_step_in_s = time_step_in_s,
                        ## following the reasoning from Jenkins & Bombosch (1995) and McCave & Swift (1976)
                        p0 =1030,             #kg/m^3 seawater density
                        p1 =1100,             #kg/m^3 Diatom density (so far a quick-look-up-average density from Ierland & Peperzak (1984))
                        cosO =1,              #cos theta, which equals 1 for sedimentation on the flat seafloor
                        g =9.81,              #accelaration due to gravity
                        K =0.0025,            #drag coefficient
                        E =1,                 #aspect ration of settling flocks (spherical = 1)
                        r =0.00016,           #particle-radius
                        Ucsq =-(0.05*(p0-p1)*g*2*(1.5*E)^(1/3)*r)/(p0*K),
                        Wd =speed,#/24/3600,
                        SedFunct =function(U_div,dens) time_step_in_s*(p1*(dens)*Wd*cosO*(U_div)*(U_div))/p0){
  list(p0 = p0, p1 = p1, cosO = cosO, g = g, K = K, E = E, r = r, Wd = Wd, Ucsq = Ucsq, SedFunct = SedFunct, time_step_in_s = time_step_in_s)
}

