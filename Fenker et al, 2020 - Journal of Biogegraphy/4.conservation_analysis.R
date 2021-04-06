# This script is a modified version of the script available in Dan's Rosauer github,
# adapted for the paper Fenker et al, 2020 - published on Journal of Biogeography
# https://doi.org/10.5061/dryad.m0cfxpp05 

## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent

####     Jessica, addapted from         ####
####       Dan Rosauer 2016             ####
#### Australian National University     ####
#### dan.rosauer@anu.edu.au             ####

rm(list=ls())

library(raster)


################################################################################
#first, import the data

PE_all <- list.files(path = "/Users/jessica/Models/diversity_analysis/", pattern = "_logPE.asc")
e <- extent(-80, -33, -35, 8)
cerrado<-shapefile ("/Users/jessica/shapes/cerrado2004.shp")
us<- shapefile ("/Users/jessica/shapes/UC federais_BR/UCs_WWF_merge_PI_final.shp")
pi<- shapefile ("/Users/jessica/shapes/UC federais_BR/UCs_WWF_merge_US_final.shp")
indg<-shapefile ("/Users/jessica/shapes/terras indigenas_BR/ti_sirgas2000.shp")
  

################################################################################
#overlap the PE maps with the conservation units and indigenous lands
#if you want, repeat the same with the other diversity variables (e.g. PD)


#aling all results
resultsPE <- list()

for(i in 1:length(PE_all)) {
  r <-raster(PE_all[i]) # raster(files[i])
  #r <- sum(r, na.rm=T)
  re <- extend(r,e)
  rc <- crop(re, e)
  
  # commented for reproducible example      
  resultsPE[[i]] <- rc # rw <- writeRaster(rc, outfiles[i], overwrite=TRUE)
  # print(outfiles[i])
  
}

#sum all the results
PE_all <- stack(resultsPE)
PE_sum <- sum(PE_all[], na.rm=T)
PE_all_sum <- stackApply(PE_all, rep(1,nlayers(PE_all)), fun = sum, na.rm = T)

#plot area of interest (Cerrado)
quartz()
PE_cer.ras <- mask(PE_all_sum, cerrado)
plot(PE_cer.ras, xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
writeRaster(PE_cer.ras,filename="_PE Cerrado.asc")


## plot PE found only within PI reserves
##
quartz()
PE_pi.ras <- mask(PE_all, pi)
protectedPE <- sum(PE_pi.ras[], na.rm=T)
PE_pi.ras <- sum(PE_pi.ras, na.rm=T)
PE_pi.ras <- mask(PE_pi.ras, cerrado)
plot(PE_pi.ras, main="PE masked to PI reserves",xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
writeRaster(PE_pi.ras, filename = "_PE_uc_PI.asc")
cat("\nTotal PE", PE_sum, "\tProtected PE", protectedPE, "\tProportion protected", protectedPE / PE_sum)


## plot PE found only within US reserves
##
quartz()
PE_us.ras <- mask(PE_all, us)
protectedPE2 <- sum(PE_us.ras[], na.rm=T)
PE_us.ras <- sum(PE_us.ras, na.rm=T)
PE_us.ras <- mask(PE_us.ras, cerrado)
plot(PE_us.ras, main="PE masked to US reserves",xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
writeRaster(PE_us.ras, filename = "_PE_uc_US.asc")
cat("\nTotal PE", PE_sum, "\tProtected PE", protectedPE2, "\tProportion protected", protectedPE2 / PE_sum)


## plot PE found only in indigenous areas
##
quartz()
PE_indg.ras <- mask(PE_all, indg)
protected_indgPE <- sum(PE_indg.ras[], na.rm=T)
PE_indg.ras <- sum(PE_indg.ras, na.rm=T)
PE_indg.ras <- mask(PE_indg.ras, cerrado)
plot(PE_indg.ras, main="PE masked to indigenous areas",xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
writeRaster(PE_indg.ras, filename = "_PE_indg_nov.asc")
cat("\nTotal PE", PE_sum, "\tProtected indigenous PE", protected_indgPE, "\tProportion protected", protected_indgPE / PE_sum)

##total PE from conservation areas of US + PI
##
quartz()
PE_uc_total.ras <- sum(PE_pi.ras, PE_us.ras)
plot(PE_uc_total.ras, main="PE masked to PI and US areas",xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
protected_total <- sum(protectedPE,protectedPE2)
writeRaster(PE_uc_total.ras, filename = "_PE_uc_total_PI_US.asc")
cat("\nTotal PE", PE_sum, "\tProtected PE total", protected_total, "\tProportion protected", protected_total / PE_sum)


##total PE from uc total (PI and US) plus indigenous areas
##
quartz()
PE_all_total.ras <- sum(PE_uc_total.ras,PE_indg.ras)
plot(PE_all_total.ras, main="PE masked to PI and US and indigenous areas",xlim=c(-60,-40), ylim=c(-25,0))
plot(cerrado, add = T)
protected_all_total <- sum(protected_total,protected_indgPE)
writeRaster(PE_all_total.ras, filename = "_PE_uc_total_PI_US_indg.asc")
cat("\nTotal PE", PE_sum, "\tProtected PE total", protected_all_total, "\tProportion protected", protected_all_total / PE_sum)

#Done!
