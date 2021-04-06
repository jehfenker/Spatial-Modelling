#### Copyright Dan Rosauer 2016         ####
#### Australian National University     ####
#### September 2012 - November 2016     ####
#### dan.rosauer@anu.edu.au             ####

library(ggplot2)

# load sizes of google maps.  This could be hard coded to remove need for file
gmap_dimensions <- read.csv("/Users/moritzlab/Dropbox/JéLeo/PE/analysis/Code originals Dan/gmap_zoom_dimension.csv")

strip_illegal_characters <- function(instring) {
  illegals <- c("[\"?\\/:*?<>|]")
  outstring <- str_replace_all(instring, illegals, "")
  return(outstring)
}

choose_gmap_zoom <- function(width, height) {
  wide_enough <- which(gmap_dimensions$xwidth > width)
  high_enough <- which(gmap_dimensions$yheight > height)
  big_enough  <- intersect(wide_enough, high_enough)
  zoom        <- max(big_enough)
  return(zoom)
}

ggplot_bounds <- function(g_plot) {
  # this function takes a ggplot object and returns its x and y limits
  # as a vector of xmin, xmax, ymin, ymax

  #check for correct argument type
  if (("ggplot" %in% class(g_plot)) == FALSE) {
    warning("Function ggplot_bounds requires a ggplot object.\n")
    return(NA)
  }

  m_build <- ggplot_build(g_plot)
  ranges <- m_build[[2]]$ranges
  ranges <- ranges[[1]]
  x_range <- ranges$x.range
  y_range <- ranges$y.range
  result <- c(x_range, y_range)
  names(result) <- c("xmin", "xmax", "ymin", "ymax")
  return(result)
}


map_raster = function(raster, output_file, title, xlimits=NA, ylimits=NA, poly=NA, filetype='png', image.width=2100, image.height=1800) {

  p         <- rasterToPoints(raster)
  p         <- data.frame(p)
  names(p) <- c("x", "y", "Model")
  colour_gradient <- scale_fill_gradientn(colours = rainbow(15), values=p$model)
  #colour_gradient <- scale_fill_gradient2(low="white", mid="yellow", high="red",
  #                                        limits=c(min(p$Model),max(p$Model)), midpoint=quantile(p$Model, 0.75), space='Lab')
  my.aes <- aes(x, y, fill=Model)

  m <- ggplot(data=p) + geom_tile(my.aes) + coord_equal() + labs(x=NULL, y=NULL) + colour_gradient

  #m <- m + ggplot(data=shape2, aes(long, lat, group=group)) + geom_polygon()

  if (! is.na(xlimits[1])) {
    m <- m + xlim(xlimits)
  }

  if (! is.na(ylimits[1])) {
    m <- m + ylim(ylimits)
  }

  if (class(poly)=="SpatialPolygonsDataFrame") {  # this is not working correctly!
    poly <- fortify(poly)
    m <- m + aes(data=poly) + geom_polygon() + geom_path(color="black") + xlim(xlimits)
  }

  # delete a previous file if needed
  if (file.exists(output_file)) {
    file.remove(output_file)
    cat("Previous", output_file, "removed\n")
  }

  m <- m + ggtitle(title)
  m <- m + theme(axis.title=element_text(face="bold", size="18"))
  m <- m + theme(axis.text=element_text(face="bold", size="14"))
  m <- m + theme(plot.title=element_text(face="bold", size="24"))
  m <- m + xlab("longitude") + ylab("latitude")

  if (filetype=='png') {
    png(output_file, width=image.width, height=image.height)
  } else if (filetype=='pdf') {
    pdf(output_file, width=image.width, height=image.height)
  }
  print(m)
  dev.off()
  m <- NULL
}

