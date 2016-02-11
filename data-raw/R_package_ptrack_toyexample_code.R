###########################################################################################################################
## read in parameters from netcdf-files and grids
## adjust dimensions of variables
## setup kdtrees
###########################################################################################################################

## to create the example-data for the package (toyROMS)

# ## read in ROMS-file
# #e <- "misom016_grd.nc"
e <- "ToyROMS58.nc"
#e <- "ToyROMS59.nc"
#e <- "ToyROMS60.nc"
#e <- "ToyROMS61.nc"

nc_e <- nc_open(e)
## u,v,rho distinction is ignored below
lon_u <- ncvar_get(nc_e, "lon_u")
lat_u <- ncvar_get(nc_e, "lat_u")
angle <- ncvar_get(nc_e, "angle")
h <- ncvar_get(nc_e, "h")               # h is bathymetry, we need Cs function of w to determine every cell's depth
Cs_w <- ncvar_get(nc_e, "Cs_w")
i_u_raw <- ncvar_get(nc_e, "u", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
i_v_raw <- ncvar_get(nc_e, "v", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
i_w <- ncvar_get(nc_e, "w", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
nc_close(nc_e)

i_w <- i_w[,,-1]  # w is one vertical dimension more as its located on top and bottom of bins
Cs_w <- Cs_w[-32]       # same for Cs_w... this removes the uppermost layer
Cs_wr <- rev(Cs_w)      # top to bottom

## correction for angle in ROMS
i_u <- array(NA, dim=c(50,40,31))
i_v <- array(NA, dim=c(50,40,31))
for(idim in 1:31){
  i_u[,,idim] = i_u_raw[,,idim]*cos(angle)-i_v_raw[,,idim]*sin(angle)
  i_v[,,idim] = i_u_raw[,,idim]*sin(angle)+i_v_raw[,,idim]*cos(angle)
}

## populate the depths in one big array [lon lat w], remember h is one row extra in the u dim
hh <- array(rep(as.vector(h), length(Cs_wr)) * rep(Cs_wr, each = prod(dim(h))), c(dim(h), length(Cs_wr)))

# toyROMS_data <- list(lon_u,lat_u,h,hh,i_u,i_v,i_w)
# names(toyROMS_data) <- c("lon_u","lat_u","h","hh","i_u","i_v","i_w")
# save(toyROMS_data, file="toyROMS_data.Rdata")

# toyROMS58 <- list(lon_u,lat_u,h,hh,i_u,i_v,i_w)
# names(toyROMS58) <- c("lon_u","lat_u","h","hh","i_u1","i_v1","i_w1")
# save(toyROMS58, file="toyROMS58.Rdata")

# toyROMS59 <- list(lon_u,lat_u,h,hh,i_u,i_v,i_w)
# names(toyROMS59) <- c("lon_u","lat_u","h","hh","i_u2","i_v2","i_w2")
# save(toyROMS59, file="toyROMS59.Rdata")

# toyROMS60 <- list(lon_u,lat_u,h,hh,i_u,i_v,i_w)
# names(toyROMS60) <- c("lon_u","lat_u","h","hh","i_u3","i_v3","i_w3")
# save(toyROMS60, file="toyROMS60.Rdata")

# toyROMS61 <- list(lon_u,lat_u,h,hh,i_u,i_v,i_w)
# names(toyROMS61) <- c("lon_u","lat_u","h","hh","i_u4","i_v4","i_w4")
# save(toyROMS61, file="toyROMS61.Rdata")
