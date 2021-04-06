#### Copyright Jessica Fenker 2020   ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####


## Take MaxEnt models created in previous script and project to different time points.
## This script uses some elements from the previous script, called "SDMs_maxent-in-dismo_jess.R"
## You need to have the time periods that you are interested previously set up in your computer

#install necessary packages
install.packages("EnvStats")

#call the packages
library(raster)
library(dismo)
library(EnvStats)
library(dplyr)

#set file paths
output_dir <- "~/jessica/running_SDMs/SDM_results/" #change details to your computer folder
time_folders <- "~/jessica/running_SDMs/past_layers/" #say where your saved the past layers


variables <- c("bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")
time_periods <- c("001", "002", "003", "004", "005", "006", "007", "008", 
                  "009", "010", "011", "012", "014", "016", "018", "020", "022", 
                  "024", "026", "028", "030", "032", "034", "036", "038", "040", 
                  "042", "044", "048", "050", "052", "054", "056", "058", "060", 
                  "062", "064", "066", "068", "070", "072", "074", "076", "078", "080") # asnamed in your folder

#create function to project species
project_species <- function(species, time, output_dir, predictors){
  print(paste0("beginning ", species, " for ", time))
  species_folder <- paste0(output_dir, species)
  load(paste0(species_folder, "/", species,"_maxent.RData"))
  px <- dismo::predict(predictors, xm, ext=ext, progress='')
  writeRaster(px, filename=paste0(species_folder, "/", species, "_maxent_", time, ".asc"), format="ascii", overwrite=TRUE)
}

#go thru table and get only species where AUC more than 0.7
model_eval <- read.csv(paste0(output_dir, "maxent_model-evals_species.csv"))
head(model_eval)
species_list <- na.omit(model_eval)
species_list <- species_list[species_list$AUC > 0.7,] 
species_list <- species_list[,2]


#read in geological data
geo_data <- stack(bioclim_vars)
geo_data <- crop(geo_data, extent_AMT)
plot(geo_data)


#loop thru time periods
for (i in time_periods){
  
  time_dir <- paste0(time_folders, i)
  time_list <- list.files(time_dir, pattern=".asc", full.names=TRUE)
  clim_data <- stack(time_list)
  
  #crop stack to your area extension
  predictors <- crop(clim_data, extent_AMT)
  crs(predictors) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  #if you want, give names to your predictors
  pred_names <- c("bioclim_04","bioclim_10","bioclim_11","bioclim_12","bioclim_15","bioclim_16","bioclim_17")
  names(predictors) <- pred_names
  
  #loop thru species and read in model
  lapply(species_list, project_species, time=i, output_dir=output_dir, predictors=predictors)
  
}


##refugia
library(EnvStats)
species_folder <- "~/jessica/running_SDMs/SDM_results/"


for (i in species_unique){
  
  sp_dir <- paste0(species_folder, i)
  sp_list <- list.files(sp_dir, pattern=".asc", full.names=TRUE)
  sp_data <- stack(sp_list)
  
  crs(sp_data) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  #calculate mean, standard deviation, and geographic mean, saving the raster
  calc_mean <- calc(sp_data, mean)
  writeRaster(calc_mean, filename=paste0(species_folder, "/", i, "_mean_refugia", ".asc"), format="ascii", overwrite=TRUE)
  calc_sd <- calc(sp_data, sd)
  writeRaster(calc_sd, filename=paste0(species_folder, "/", i, "_sd_refugia", ".asc"), format="ascii", overwrite=TRUE)
  calc_geo_mean <- calc(sp_data, fun = geoMean)
  writeRaster(calc_geo_mean, filename=paste0(species_folder, "/", i, "_geo_mean_refugia", ".asc"), format="ascii", overwrite=TRUE)
  # plot and save a pdf 
  pdf(paste0(species_folder, "/", species, "_mean_refugia.pdf"), paper="a4r")
  plot(calc_mean, main='Mean values of Refugia')
  pdf(paste0(species_folder, "/", species, "_sd_refugia.pdf"), paper="a4r")
  plot(calc_sd, main='Standard Deviation values of Refugia')
  pdf(paste0(species_folder, "/", species, "_geo_mean_refugia.pdf"), paper="a4r")
  plot(calc_geo_mean, main='Geometric Mean values of Refugia')
  
  print(paste0(i, " finished"))
  
}

# Done!
