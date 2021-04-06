# This script is related to the publication
# Fenker et al, 2020 - Journal of Biogeography
# https://doi.org/10.5061/dryad.m0cfxpp05 

####      Jessica Fenker 2020        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####


## This script will combine the PE results with conservation areas
## and indigenous lands in Cerrado 

rm(list=ls())

library(raster)

################################################################################

work.dir   = '/Users/jessica/Models/'
setwd(work.dir)

PE.stack <- stack('PE_ALLgenera_aligned_stack.grd')
region.shp <- shapefile("/Users/jessica/shapes/cerrado2004.shp")
class <- 'cerrado'
cerrado_classes <- data.frame(unique(region.shp$GovIndOth), stringsAsFactors = F)
names(cerrado_classes) <- 'class'
cerrado_sums <- data.frame(names(PE.stack), stringsAsFactors=F)
names(cerrado_sums) <- "genus_name"
cerrado_areas <-data.frame(cerrado_classes, stringsAsFactors=F)
names(cerrado_areas) <- "class"

for (class in cerrado_classes$class) {

  class_shape <- cerrado.shp[region.shp$GovIndOth==class,]
  class.stack <- crop(PE.stack, y = class_shape)
  class.stack <- mask(class.stack, mask = class_shape)

  class.array <- as.array(class.stack)
  colname   <- paste(class, 'sum', sep="_")
  cerrado_sums[, colname] <- apply(class.array, MARGIN = 3, FUN = sum, na.rm=T)

  colname   <- paste(class, 'mean', sep="_")
  cerrado_means[, colname] <- apply(class.array, MARGIN = 3, FUN = mean, na.rm=T)

  cerrado_areas$area[which(cerrado_classes == class)] <- length(which(!is.na(class.array[,,1])))

  cat(class, "\n")

}

# now get sums for the whole area

class.stack <- crop(PE.stack, y = region.shp)
class.stack <- mask(class.stack, mask = region.shp)

class.array <- as.array(class.stack)
colname   <- paste(class, 'sum', sep="_")
cerrado_sums[, colname] <- apply(class.array, MARGIN = 3, FUN = sum, na.rm=T)

colname   <- paste(class, 'mean', sep="_")
cerrado_means[, colname] <- apply(class.array, MARGIN = 3, FUN = mean, na.rm=T)

cerrado_areas[4,1] <- 'cerrado'
cerrado_areas$area[4] <- length(which(!is.na(class.array[,,1])))


cerrado_areas$PE_sum  <- apply(cerrado_sums[,2:5],  MARGIN=2, FUN=sum, na.rm=T)
cerrado_areas$PE_mean <- apply(cerrado_means[,2:5], MARGIN=2, FUN=sum, na.rm=T)

cerrado_areas_prop <- cerrado_areas
cerrado_areas_prop$area <- cerrado_areas_prop$area / cerrado_areas_prop$area[4]
cerrado_areas_prop$PE_sum <- cerrado_areas_prop$PE_sum / cerrado_areas_prop$PE_sum[4]

write.csv(cerrado_sums, "class_PE_sums.csv", row.names=F)
write.csv(cerrado_areas_prop, "class_PE_proportions.csv", row.names=F)

