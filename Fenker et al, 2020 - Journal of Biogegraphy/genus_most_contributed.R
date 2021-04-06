# This script is related to the publication
# Fenker et al, 2020 - Journal of Biogeography
# https://doi.org/10.5061/dryad.m0cfxpp05 

####      Jessica Fenker 2020        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####

## The aim of this script to figure out the genera that most contributed to the result
## and generate a nice plot

rm(list=ls())
library(raster)
library(rasterVis)

#define directories
working.dir   = '/Users/jessica/Models/'
template_ext ='/Users/jessica/Models/species_models/Ameiva/ameiva_median.asc'

setwd(working.dir)

#stack of the PE results for all genera
PE.stack <- stack("PE_ALLgenera_aligned_stack.grd")

#use stackApply to calculate the maximum value of each pixel
maxval  <- stackApply(PE.stack, rep(1,nlayers(PE.stack)), fun=max, overwrite=TRUE)
summary(maxval)

#now, use stackApply to calculate the sum
sum     <- stackApply(PE.stack, rep(1,nlayers(PE.stack)), fun=sum)

#divide the maximum value by the sum
maxprop <- maxval / sum
summary(maxprop)
writeRaster(maxprop,"highest_genus_PE_proportion.asc")

max_genus  <- which.max(PE.stack)
summary(max_genus)
summary(PE.stack)
table(max_genus@data@values)
writeRaster(max_genus,"max_genus_PE_proportion.asc")


#Some extra analysis
max_genus <- raster('/Users/jessica/Models/max_genus_PE_proportion.asc')
template.ras <- raster(template_ext)
max_genus  <- mask(max_genus, template.ras)
plot(max_genus)

writeRaster(max_genus,"max_genus_revised.asc", overwrite=TRUE)
stack_names <- data.frame(cbind(1:nlayers(PE.stack),names(PE.stack)))
names(stack_names) <- c("layer_num","layer_name")
write.csv(stack_names,"layer_names")

max_genus_above_med <- max_genus
max_genus_above_med[sum < 0.0000454] <- NA  # actually above log(10) which is well above median
writeRaster(max_genus_above_med,"max_genus_above_median.asc", overwrite=TRUE)

plot(max_genus_above_med)

# do something here to work out how high a location is across multiple genera
sum_minus_max <- sum - maxval
writeRaster(sum_minus_max,"PE_without_max_genus.asc", overwrite=TRUE)

log_sum_minus_max <- log(sum_minus_max)
writeRaster(log_sum_minus_max,"logPE_without_max_genus.asc", overwrite=TRUE)


median_PE     <- stackApply(PE.stack,rep(1,nlayers(PE.stack)),fun=median)
logmedian_PE  <- log(median_PE)
writeRaster(logmedian_PE,"log_median_PE.asc", overwrite=TRUE)





# If you want the values for a specific region


region_shapefile <- "/Users/jessica/shapes/cerrado2004.shp"
region_name      <- "Cerrado"

region.shp <- shapefile(region_shapefile)
PE.stack <- crop(PE.stack, region.shp )
PE.stack <- mask(PE.stack, region.shp)

maxval  <- stackApply(PE.stack, rep(1,nlayers(PE.stack)), fun=max)
sum     <- stackApply(PE.stack,rep(1,nlayers(PE.stack)), fun='sum')
plot(maxval)

median_PE <- median(sum[], na.rm=T)
max_genus  <- which.max(PE.stack)
max_genus_above_med <- max_genus
max_genus_above_med[sum < median_PE] <- NA
writeRaster(max_genus_above_med,"max_genus_above_median_Cerrado.asc", overwrite=TRUE)
writeRaster(max_genus,"max_genus_Cerrado.asc", overwrite=TRUE)


#Based on the table and values obtained, generate an improved image specifying the colors you want for each genera
quartz()
colors <-  colorRampPalette(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"))(255)
plot(max_genus, legend = FALSE, col = rev(colors))
legend("topleft", legend = c("Tropidurus","Phyllopezus", "Norops", "Enyalius", "Cercosaura", "Colobosaura"), fill = rev(colors))
plot(max_genus_above_med)

sum_minus_max <- sum - maxval
writeRaster(sum_minus_max,"PE_without_max_genus_Cerrado.asc", overwrite=TRUE)
plot(sum_minus_max)
log_sum_minus_max <- log(sum_minus_max)
writeRaster(log_sum_minus_max,"logPE_without_max_genus_Cerrado.asc", overwrite=TRUE)
plot(log_sum_minus_max)

without_max.stack <- PE.stack

for (i in 1:nrow(PE.stack)) {
  top_cells <- which(max_genus[]==i)
  without_max.stack[top_cells,i] <- 
}

head(max_genus[])
