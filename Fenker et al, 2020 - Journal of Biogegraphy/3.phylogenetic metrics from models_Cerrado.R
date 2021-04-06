# This script is a modified version of the script available in Dan's Rosauer github,
# adapted for the paper Fenker et al, 2020 - published on Journal of Biogeography
# https://doi.org/10.5061/dryad.m0cfxpp05 

## this script calculates richness and endemism from modelled suitability surfaces
## it requires all the model grids to have the same extent

####       Dan Rosauer 2016             ####
#### Australian National University     ####
#### dan.rosauer@anu.edu.au             ####

rm(list=ls())

library(SDMTools)
library(raster)
library(ape)
library(phylobase)
library(foreach)
library(parallel)
library(ggplot2)

rm(list=ls())
source("/Users/jessica/Code_Dan/phylogenetic endemism.r")

################################################################################
#first define some functions

map_raster = function(raster, output_file, title) {

  p         <- rasterToPoints(raster)
  p         <- data.frame(p)
  names(p) <- c("x", "y", "Model")
  colour_gradient <- scale_fill_gradientn(colours = rainbow(15), values=p$model)
  #colour_gradient <- scale_fill_gradient2(low="white", mid="yellow", high="red",
  #                                        limits=c(min(p$Model),max(p$Model)), midpoint=quantile(p$Model, 0.75), space='Lab')
  m <- ggplot(data=p) + geom_tile(aes(x, y, fill=Model)) + coord_equal() + labs(x=NULL, y=NULL) + colour_gradient

  # delete a previous file if needed
  if (file.exists(output_file)) {
    file.remove(output_file)
    cat("Previous", output_file, "removed...")
  }

  m <- m + ggtitle(title)
  m <- m + theme(axis.title=element_text(face="bold", size="18"))
  m <- m + theme(axis.text=element_text(face="bold", size="14"))
  m <- m + theme(plot.title=element_text(face="bold", size="24"))
  m <- m + xlab("longitude") + ylab("latitude")

  png(output_file, width=image.width, height=image.height)
  print(m)
  dev.off()
  cat("\tand new version saved.\n")
  m <- NULL
}

################################################################################

max.rows        <- 10000000
core_count      <- 2 # number of cores to use for parallel steps
write_matrices  <- FALSE

# size in pixels for maps
image.width=1400
image.height=1400

#define directories
base.dir        <- '/Users/jessica/Models/PEanalysis/'
input.dir       <- '~/Models/lineage_models/asc_aligned/'
output.dir      <- 'Diversity'
file_pattern    <- 'lin_model_Ameiva'

template_grid   <- 'lin_model_Ameiva_ameiva_clade1.asc'
group_lin_file  <- 'Diversity/group_lineage_list.csv'

#tree details  - this works for one genus at a time
tree.file       <- 'Models/trees/Ameiva_tree.tre'
outgroup        <- ''
preface         <- 'lin_model_Ameiva_'  # this is text to be removed from the model names before matching to the tree tips. May vary depending on how the models were named

genus           <- 'Ameiva'
output_prefix   <- "Ameiva_"
threshold       <- 0.01  # this is not a species level threshold, but one used for each lineage model

####  end of parameters  ####

setwd(base.dir)
files <- list.files(path = input.dir, pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

#setwd(input.dir)
#template.asc = read.asc.gz(files[1])
template_grid <- paste(base.dir, input.dir, template_grid, sep="")

template.asc = read.asc(template_grid)
model_rows=nrow(template.asc)
model_cols=ncol(template.asc)

# the original version, excluding NA cells
pos <- as.data.frame(which(is.finite(template.asc),arr.ind=TRUE)) # get all points that have data

cat("\nLoading model rasters in parallel\n")

# cl <- makeCluster(core_count)
# registerDoParallel(cl)

setwd(input.dir)

pos_par <- foreach (j=1:length(files), .combine=cbind, .packages='SDMTools') %do% {
#for (j in 1:length(files)) {
 tfile <- files[j]
  pos_temp <- pos
  checkname = unlist(strsplit(tfile,".",fixed=T))
  if (checkname[length(checkname)]=="asc") {   # only accept filenames ending in .asc

    cat(j)

    tasc = read.asc(tfile)                                #read in the data
    dataname=gsub(".asc",'',tfile)
    newname <- tolower(gsub(preface, "", dataname))

    pos_temp[newname] <- tasc[cbind(pos_temp$row, pos_temp$col)]
    pos_temp[(which(pos_temp[newname]< threshold)), newname] <- 0    # set values below the threshold to 0
    pos_temp[(which(is.na(pos_temp[newname]))), newname]     <- 0    # set the nulls to 0
    cat("\n", j, newname, "loaded")
    newcol <- data.frame(pos_temp[, newname])
    newcol <- round(newcol,3)
    names(newcol) <- newname
    newcol
  }
}

pos <- cbind(pos,pos_par)
rm(pos_par, pos_temp)

setwd(base.dir)
group_lin_list <-read.csv(group_lin_file, stringsAsFactors=F)
group_lin_list$lineage <- tolower(group_lin_list$lineage)
group_lin_list <- group_lin_list[group_lin_list$genus == genus,]

# read in the tree from one of two formats
tree <- read.nexus(tree.file)
#tree <- read.tree(tree.file)

if (outgroup != '') {
  tree <- drop.tip(tree,tip=which(tree$tip.label==outgroup))
}
tree <- phylo4(tree)
labels(tree) <- tolower(labels(tree))

# ensure that the tree tips match the model names

model.names <- names(pos)
model.names <- model.names[-(1:2)] # names of all columns except the 1st two which are row, col
preface <- tolower(preface)
model.names <- tolower(gsub(preface, "", model.names))
model.groups <- data.frame(model.groups=vector("character",nTips(tree)),stringsAsFactors=F)

for (i in 1:nTips(tree)) {
  cat(labels(tree)[i],"\n")
  row <- group_lin_list[group_lin_list$lineage==labels(tree)[i],]
  if (nrow(row) > 0) {
    new_tip_name <- tolower(paste(row$genus, row$model_group, row$lineage, sep="_"))
    labels(tree)[i] <- new_tip_name
    model.groups[i,1] <- as.character(row$model_group)
  }
}

tree <- phylo4d(tree,tip.data=model.groups)

tree.names  <- as.character(labels(tree)[1:nTips(tree)])
matched.names <- intersect(model.names,tree.names)
matched.tips  <- which(labels(tree,"tip") %in% matched.names)
cat("\nNot in tree names:",setdiff(model.names,tree.names),"\n")
cat("\nNot in model names:",setdiff(tree.names,model.names),"\n")

# a subtree containing only the tips for which there is a corresponding model
subtree <- subset(tree,tips.include=matched.tips)

if (write_matrices) {
  cat("\nWriting the site x lineage matrix to file\n")
  write.csv(pos,paste(output.dir,"sites_x_lineage.csv",sep=''))
}

# limit the occurrence table to lineages which match the tree
names(pos) <- tolower(gsub(preface, "", names(pos)))
matching_pos_columns <- which(names(pos) %in% matched.names)
matching_pos_columns <- unique(c(1, 2, matching_pos_columns))  # ensure that row and column are included
pos <- pos[, matching_pos_columns]

cat("\n\nRemoving unoccupied cells\n")
cat("Before:",nrow(pos),"\n")
rowsums <- apply(pos[,3:ncol(pos)],1,sum,na.rm=T)
pos <- pos[which(rowsums>0),]
rm(rowsums)
cat("After:",nrow(pos),"\n")
max.rows <- min(max.rows,nrow(pos))

gc()
# send the tree and distribution data to the function which calculates PE and other metrics
result <- calc_PE_from_models(subtree, pos[1:max.rows, which(names(pos) %in% matched.names)], core_count = core_count)
gc()

cat("\nDiversity calculations completed. Now writing outputs to file.")

pos_output <- cbind(pos[1:max.rows,1:2],result)
pos_output <- pos_output[,-3] # omit the site column

# add lat and long columns
cellsize <- attr(template.asc,"cellsize")
ymin <- attr(template.asc,"yll")
xmin <- attr(template.asc,"xll")

x <- ((pos_output$row - 1) * cellsize) + xmin
y <- ((pos_output$col - 1) * cellsize) + ymin
pos_output <- cbind(pos_output[,1:2],x,y,pos_output[,-(1:2)])


# optional - add residual columns
PD_SR_mod <- lm(pos_output$PD~pos_output$SR)
pos_output$PD_SR_resid <- PD_SR_mod$residuals
PD_SR_loglog_mod <- lm(log(pos_output$PD)~log(pos_output$SR),subset=which(!is.infinite(log(pos_output$SR))))
pos_output_log1 <- cbind(pos_output[which(!is.infinite(log(pos_output$SR))),],PD_SR_loglog_mod$residuals)
output_filename1 <- paste(output_prefix,"PD_SR_loglog_resid.asc",sep="")
dataframe2asc(pos_output_log1[,c(4,3,10)],paste(output_prefix,"PD_SR_loglog_resid.asc",sep=""))

####
PE_WE_mod <- lm(pos_output$PE~pos_output$WE)
pos_output$PE_WE_resid <- PE_WE_mod$residuals
PE_WE_loglog_mod <- lm(log(pos_output$PE)~log(pos_output$WE),subset=which(!is.infinite(log(pos_output$WE))))
pos_output_log <- cbind(pos_output[which(!is.infinite(log(pos_output$WE))),],PE_WE_loglog_mod$residuals)
output_filename <- paste(output_prefix,"PE_WE_loglog_resid.asc",sep="")
dataframe2asc(pos_output_log[,c(4,3,10)],paste(output_prefix,"PE_WE_loglog_resid.asc",sep=""))

#logPE
pos_output$logPE <- log(pos_output$PE)

#creat csv
filenames <- c(paste(output_prefix,"PE.asc",sep=""), paste(output_prefix, "PD.asc", sep=""), paste(output_prefix, "WE.asc", sep=""), paste(output_prefix, "SR.asc", sep=""), paste(output_prefix, "logPE.asc", sep="") , paste(output_prefix, "PE_WE_resid.asc", sep=""),paste(output_prefix, "PD_SR_resid.asc", sep=""))

columns_to_use <- c(5:11)
dir <- paste(base.dir, output.dir, sep="")
dataframe2asc(pos_output[, c(4, 3, columns_to_use)], filenames, dir)

write.csv(pos_output, paste(dir, "/", output_prefix, "scores.csv", sep=""), row.names=FALSE)

# stopCluster

setwd(dir)

cerrado <- shapefile("/Users/jessica/shapes/cerrado2004.shp")

# make some output images
PE.ras <- raster(filenames[1])
PE.ras <- mask(PE.ras, cerrado)
map_filename <- paste(output_prefix, "PE.png", sep="")
map_raster(PE.ras, map_filename, paste(output_prefix, "PE"))

# log PE
PElog.ras <- raster(filenames[5])
PElog.ras <- mask(PElog.ras, cerrado)
map_filename <- paste(output_prefix, "log_PE.png", sep="")
map_raster(PElog.ras, map_filename, paste(output_prefix, "log PE"))

PD.ras <- raster(filenames[2])
PD.ras <- mask(PD.ras, cerrado)
map_filename <- paste(output_prefix, "PD.png", sep="")
map_raster(PD.ras, map_filename, paste(output_prefix, "PD"))

WE.ras <- raster(filenames[3])
WE.ras <- mask(WE.ras, cerrado)
map_filename <- paste(output_prefix, "WE.png", sep="")
map_raster(WE.ras, map_filename, paste(output_prefix, "WE"))

SR.ras <- raster(filenames[4])
SR.ras <- mask(SR.ras, cerrado)
map_filename <- paste(output_prefix, "SR.png", sep="")
map_raster(SR.ras, map_filename, paste(output_prefix, "SR"))

PD_SR_resid.ras <- raster(filenames[7])
PD_SR_resid.ras <- mask(PD_SR_resid.ras, cerrado)
map_filename <- paste(output_prefix, "PD_SR_resid.png", sep="")
map_raster(PD_SR_resid.ras, map_filename, paste(output_prefix, "PD_SR_resid"))

PE_WE_resid.ras <- raster(filenames[6])
PE_WE_resid.ras <- mask(PE_WE_resid.ras, cerrado)
map_filename <- paste(output_prefix, "PE_WE_resid.png", sep="")
map_raster(PE_WE_resid.ras, map_filename, paste(output_prefix, "PE_WE_resid"))

map_filename <- paste(output_prefix, "_4maps.png", sep="")
png(map_filename, 1600, 1200)
par(mfrow=c(2,2),mar=c(3,4,3,2))
plot(PE.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="PE", col=rainbow(25,start=0.1,end=1))
plot(PD.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="PD", col=rainbow(25,start=0.1,end=1))
plot(WE.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="WE", col=rainbow(25,start=0.1,end=1))
plot(SR.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="SR", col=rainbow(25,start=0.1,end=1))
plot(PD_SR_resid.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="PD_SR_resid", col=rainbow(25,start=0.1,end=1))
plot(PE_WE_resid.ras, xlim = c(-60,-42), ylim = c(-25,-2), main="PE_WE_resid", col=rainbow(25,start=0.1,end=1))

dev.off()

