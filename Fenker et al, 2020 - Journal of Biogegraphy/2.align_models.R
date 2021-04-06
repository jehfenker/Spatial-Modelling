# After ran the first script, I had some problems with models that were not aligned.
# This script solved the issue
# This script is related to the publication
# Fenker et al, 2020 - Journal of Biogeography
# https://doi.org/10.5061/dryad.m0cfxpp05 


library(raster)
library(stringr)

### Parameters ###
base.dir          <- "/Users/jessica/Models/"   # modify to the base directory for your lineage modelling

input.dir         <-  paste(base.dir, 'lineage_models/asc/', sep ='')          # location of existing lineage distribution models
output.dir        <-  paste(base.dir, 'lineage_models/asc_aligned/', sep='')   # now folder, where the aligned lineage distribution models will be saved
output_stack_name <-  "lineage_models"
new_extent        <-  extent(-80, -33, -35, 8)

new_only          <- FALSE  # if true, skip grids which are already in the output directory
file.pattern      <- '*.asc$'  #regex
### End Parameters ###

setwd(input.dir)

# check that the output directory exists and if not, create it
if (! dir.exists(output.dir)) {
  dir.create((output.dir))
  cat("\n\nCreated new directory for lineage models:", output.dir, "\n\n")
}

input_files   <-  list.files(path=input.dir, pattern=file.pattern, full.names=FALSE, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)
output_files  <-  list.files(path=output.dir, pattern=file.pattern, recursive=FALSE, ignore.case=TRUE, include.dirs=FALSE)

raster_names <- ""

for (tfile in input_files) {
  filepath = paste(input.dir,tfile,sep='')
  outname = paste(output.dir,tfile,sep="")

  if ((!tfile %in% output_files) | (! new_only)) {
    grid.ras = raster(tfile)
    grid_ext.ras = extend(grid.ras,new_extent) # extend to the union of current grid and new extent
    grid_ext.ras = crop(grid_ext.ras,new_extent) # crop back to new extent
    writeRaster(grid_ext.ras,outname,overwrite=TRUE, NAflag=-9999)
    cat("\nAligned asc raster written for",tfile)

    # make a vector of the layer names
    if (length(raster_names) == 1 & raster_names[1] == "") {
      raster_names <- outname
    } else {
      raster_names <- c(raster_names,outname)
    }

  } else {
    cat("\nSkipped",tfile)
  }
}

# now make a raster stack
raster_file <- list.files(path='/Users/jessica/Models/lineage_models/asc_aligned/', pattern =".asc", full.names=TRUE)

lin.stack <- stack(raster_file)
writeRaster(lin.stack, output_stack_name, overwrite=TRUE)

setwd(output.dir)

sum.ras     <- stackApply(lin.stack, rep(1,nlayers(lin.stack)), fun=sum)
writeRaster(sum.ras,"sum.asc", overwrite=TRUE)

stack_names <- data.frame(cbind(1:nlayers(lin.stack),names(lin.stack)))
names(stack_names) <- c("layer_num","layer_name")
stack_names$layer_name <- str_replace(stack_names$layer_name, pattern = "lin_model_", replacement = "")
write.csv(stack_names, "layer_names.csv", row.names = FALSE)

cat("\n\nFinished writing aligned model rasters to", output.dir, "\n\n")
