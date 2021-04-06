#######################################
####### SDM script - Jess Fenker ######
#######################################

#install necessary packages
install.packages("dismo")
install.packages("rJava")
install.packages("svMisc")

#call the packages
library(dismo)
library(sp)
library(rgeos)
library(raster)
library(rJava)
library(svMisc)
library(EnvStats)
library(dplyr)
library(ncdf4)
library(lattice)
library(ggplot2)

#if you didnt create a project,set your working directory
#setwd("~/Box/running_SDMs")


#To use maxent you must first download the program from 
# https://biodiversityinformatics.amnh.org/open_source/maxent/
#Put the file 'maxent.jar' in the 'java' folder of the 'dismo' package. 
#That is the folder returned by system.file("java", package="dismo").
#Then, go to Finder and in the upper menu select Go/Go to Folder... and copy the file path

#check if maxent is in the folder
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='') 
if (file.exists(jar)) {
  cat("can continue, maxent is available")} else {
  cat('cannot run this because maxent is not available')}


##############################
## GET ENVIRONMENTAL LAYERS ##
#############################
file <- "LateQuaternary_Environment.nc"

var_names <- c('npp', 'BIO1', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13','BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19')  
###
for (var_name in var_names) {
  
  # Create raster stack
  x1 <- brick(file, varname = var_name)
  ext <- c(112, 155, -45, -8) #extension of Australia
  x1 <- crop(x1, ext)
  
  #create netcdf file
  writeRaster(x = x1, 
              filename = paste0(var_name, '_out.nc'),
              overwrite = TRUE, 
              format = 'CDF')
  }

# this command will save a netcdf file of each of your variable in your directory folder


##############################
###  ADDING SPECIES DATA ####
#############################


#read species data and crop to area 
species_data_all <-read.csv("dipes_coord.csv")
head(species_data_all)
species_data_unique <-read.csv("dipes_species.csv")
species_unique <- unique(species_data_unique$species)
#for fast visualization of your data
plot(species_data_all$x, species_data_all$y)

##############################
### RUNNING PRESENT SDM  ####
#############################


time.dir <- "~/Box/running_SDMs/"
var_names <- c('npp', 'BIO1', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10', 'BIO11', 'BIO12', 'BIO13','BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19')  


x.times.list <- list()
for (j in 1:72){
  x.list <- list()
  for (i in 1:length(var_names)){
    x.tmp <- brick(paste0(var_names[i], "_out.nc"))
    x.list[[i]] <- x.tmp[[j]]
  }
  names(x.list) <- var_names
  x.times.list[[j]] <- stack(x.list)
  names(x.times.list)[j] <- names(x.tmp)[j]
}

#put all the present variables (time= xX0) in the stack
predictors<-stack(x.times.list$X0)
extent_AMT <- c(122, 142, -22, -7) #extent of your working area
predictors <- crop(predictors,extent_AMT)
plot(predictors)

#set up a folder where your results will go
output_file <- "~/Box/SDM/running_SDMs/SDM_results/"

#set up dataframe for model evaluation data
model_eval <- data.frame(species=character(), presences=numeric(), absences=numeric(), AUC=numeric(), cor=numeric())

#loop through species and fit models
for (species in species_unique[1:length(species_unique)]){
  
  print(paste0("beginning ", species, ", background "))
  
  #get species data
  species_data <- subset(species_data_all, full_name %in% species) 
  species_data <- species_data[,1:2]
  coordinates(species_data)<-~x+y #transform dataframe to sp object
  crs(species_data) <- crs(predictors) #projection
  
  #crop background for each species
  original <- predictors
  buffer <- gBuffer(species_data, width = 3) #set up your buffer size; it depends on your species (width=degrees)
  background <- mask(original,buffer) # cut the layers to buffer sizes
  mask <- background[[2]]
  ext <- extent(background)
  
  print(paste0("beginning ", species, ", species number ", match(species, species_unique), " of ", length(species_unique)))
  
  if (length(species_data$x) > 5)
  {
    #create training and testing set
    group <- kfold(species_data, 5)
    pres_train <- species_data[group != 1, ]
    pres_test <- species_data[group == 1, ] #20% of data for testing
    
    #background data for training and testing set
    backg <- randomPoints(mask, n=10) ### set according to your dataset
    colnames(backg) = c('lon', 'lat')
    group <- kfold(backg, 2)
    backg_train <- backg[group != 1, ]
    backg_test <- backg[group == 1, ]
    
    #check it doesn't have a lot of NAs
    x <- extract(predictors, pres_train)
    x <- na.omit(x)
    na_count <- nrow(x)/length(pres_train$x) #need to be 0.5 or more
    
    #maxent
    xm <- maxent(background, pres_train) #factors=categorical variables. If has 'undefined columns' error make sure that the names of the factors are correct
    e <- evaluate(pres_test, backg_test, xm, background)
    px <- predict(background, xm, ext=ext, progress='')
    plot(px)
    
    #create folder and output results
    output_dir <- paste0(output_file, species)
    dir.create(output_dir)
    save(xm, file=paste0(output_dir, "/", species, "_maxent.RData"))
    writeRaster(px, filename=paste0(output_dir, "/", species, "_maxent_000",".asc"), format="ascii")
    
    #write plots to pdf
    pdf(paste0(output_dir, "/", species, "_maxent-summary.pdf"), paper="a4r")
    
    plot(px, main='Maxent, raw values')
    tr <- threshold(e, 'spec_sens')
    plot(px > tr, main='presence/absence')
    points(pres_train, pch='+')
    
    plot(xm)
    response(xm)
    
    dev.off()
    
    
    #add model data to df
    model_data <- list(species, e@np, e@na, e@auc, e@cor)
    model_eval <- rbind(model_eval, model_data, stringsAsFactors=FALSE)
    
  }
  else{
    model_data <- list(species, length(species_data$x), 0, NA, NA)
    model_eval <- rbind(model_eval, model_data, stringsAsFactors=FALSE)
    
  }
  
  
}



#check and output model evaluation table
colnames(model_eval) <- c("species", "presence", "absence", "AUC", "correlation")
head(model_eval)
write.csv(model_eval, file=paste0(output_file, "maxent_model-evals_species.csv"))




############
# paleo
#set file paths etc

output_dir <- "~/Box/SDM/running_SDMs/SDM_results/" #change details to your computer folder

time_periods <- c("X.80000","X.78000","X.76000","X.74000","X.72000","X.70000" ,"X.68000","X.66000","X.64000","X.62000","X.60000","X.58000","X.56000","X.54000",
                  "X.52000","X.50000","X.48000","X.46000","X.44000","X.42000","X.40000","X.38000","X.36000","X.34000","X.32000","X.30000","X.28000","X.26000",
                  "X.24000","X.22000","X.21000","X.20000","X.19000","X.18000","X.17000","X.16000","X.15000","X.14000","X.13000","X.12000","X.11000","X.10000",
                  "X.9000","X.8000","X.7000", "X.6000","X.5000","X.4000","X.3000","X.2000","X.1000") 

#create function to project species
project_species <- function(species, time, output_dir, predictors){
  print(paste0("beginning ", species, " for ", time))
  species_folder <- paste0(output_dir, species)
  load(paste0(species_folder, "/", species,"_maxent.RData"))
  px <- predict(predictors, xm, ext=ext, progress='')
  writeRaster(px, filename=paste0(species_folder, "/", species, "_maxent_", time, ".asc"), format="ascii", overwrite=TRUE)
}

#go thru table and get only species where AUC more than 0.7
model_eval <- read.csv(paste0(output_dir, "maxent_model-evals_species.csv"))
head(model_eval)
species_list <- na.omit(model_eval)
species_list <- species_list[species_list$AUC > 0.7,] 
species_list <- species_list[,2]



for (i in time_periods){

  time_dir <- x.times.list[[i]]
  clim_data <- stack(time_dir)
  #values(clim_data[[9]])[values(clim_data[[9]]) > 0] = 1
  
  #loop thru species
  #read in model
  lapply(species_list, project_species, time=i, output_dir=output_dir, predictors=clim_data)
  
}


##refugia
species_unique <- species_unique[c(3,6,7,8)]
species_folder <- "~/Box/SDM/running_SDMs/SDM_results/"

for (i in species_unique){
  
  print(paste0("beginning ", i))
  
  sp_dir <- paste0(species_folder, i)
  sp_list <- list.files(sp_dir, pattern=".asc", full.names=TRUE)
  sp_data <- stack(sp_list)
  
  #crop stack to extent_kimb
  crs(sp_data) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  #loop thru species
  #read in model
  calc_mean <- calc(sp_data, mean)
  writeRaster(calc_mean, filename=paste0(species_folder, "/", i, "_mean_refugia", ".asc"), format="ascii", overwrite=TRUE)
  calc_sd <- calc(sp_data, sd)
  writeRaster(calc_sd, filename=paste0(species_folder, "/", i, "_sd_refugia", ".asc"), format="ascii", overwrite=TRUE)
  calc_geo_mean <- calc(sp_data, fun = geoMean)
  writeRaster(calc_geo_mean, filename=paste0(species_folder, "/", i, "_geo_mean_refugia", ".asc"), format="ascii", overwrite=TRUE)
  pdf(paste0(species_folder, "/", species, "_mean_refugia.pdf"), paper="a4r")
  plot(calc_mean, main='Mean values of Refugia')
  pdf(paste0(species_folder, "/", species, "_sd_refugia.pdf"), paper="a4r")
  plot(calc_sd, main='Standard Deviation values of Refugia')
  pdf(paste0(species_folder, "/", species, "_geo_mean_refugia.pdf"), paper="a4r")
  plot(calc_geo_mean, main='Geometric Mean values of Refugia')
  
  print(paste0(i, " finished"))
  
}
