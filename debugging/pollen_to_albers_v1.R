pollen_to_albers <- function(pollen_ts) {
  library(sp)
  
  # Remove rows with missing coordinates
  pollen_ts <- pollen_ts[!is.na(pollen_ts$lat) & !is.na(pollen_ts$long), ]
  
  if(nrow(pollen_ts) == 0) stop("No valid coordinates to transform.")
  
  centers_pol <- data.frame(x = pollen_ts$long, y = pollen_ts$lat)
  coordinates(centers_pol) <- ~ x + y
  proj4string(centers_pol) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  centers_polA <- spTransform(centers_pol, CRS(SRS_string = "EPSG:3175"))
  
  centers_polA_df <- as.data.frame(centers_polA)
  colnames(centers_polA_df) <- c("x", "y")
  
  pollen_ts$x <- centers_polA_df$x
  pollen_ts$y <- centers_polA_df$y
  
  return(pollen_ts)
}
