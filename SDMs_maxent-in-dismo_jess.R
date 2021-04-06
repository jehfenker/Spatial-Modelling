
####      Jessica Fenker 2020        ####
#### Australian National University  ####
####    jehfenker@gmail.com          ####

## This script can be use to run MaxEnt in R, simultaniusly 
## for more than one species 

#install necessary packages
packages_to_install <- c("dismo","sp","rgeos","raster","rJava", "svMisc")
install.packages(packages_to_install)

#call the packages
library(dismo)
library(sp)
library(rgeos)
library(raster)
library(rJava)
library(svMisc)

#if you didnt create a project,set your working directory
#setwd("~/jessica/running_SDMs")


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

#read species data, including distribution coordinates, and crop to area 
species_data_all <-read.csv("dipes_coord.csv")
head(species_data_all)

species_unique <- unique(species_data_all$species)
#for fast visualization of your data
plot(species_data_all$x, species_data_all$y)


#read in rasters with current climatic variables
bioclim_vars <- list.files("~/jessica/bioclim/000/asc_only", pattern=".asc", full.names=TRUE)
bioclim_vars

predictors <- stack(bioclim_vars) # put all the variables together; can take some time
names(predictors)
#set the projection of your layers
crs(predictors) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
#extension of your working area
extent_AMT <- c(122, 142, -22, -7)
predictors <- crop(predictors,extent_AMT)
plot(predictors)


#set up the output file
output_file <- "~/jessica/running_SDMs/SDM_results/"

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
    backg <- randomPoints(mask, n=1000)
    colnames(backg) = c('lon', 'lat')
    group <- kfold(backg, 5)
    backg_train <- backg[group != 1, ]
    backg_test <- backg[group == 1, ]
    
    #check it doesn't have a lot of NAs
    x <- extract(predictors, pres_train)
    x <- na.omit(x)
    na_count <- nrow(x)/length(pres_train$x) #need to be 0.5 or more
    
    #maxent
    xm <- maxent(background, pres_train, factors='calcrete') #factors=categorical variables. If has 'undefined columns' error make sure that the names of the factors are correct
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



#check and create an output model evaluation table
colnames(model_eval) <- c("species", "presence", "absence", "AUC", "correlation")
head(model_eval)
write.csv(model_eval, file=paste0(output_file, "maxent_model-evals_species.csv"))

#Done!
