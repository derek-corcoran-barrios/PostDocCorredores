#Mean temp
library(raster)
library(rasterVis)

e <- readRDS("e.rds")
now <- getData('worldclim', var='bio', res=10)
now <- crop(now, e)
ac50 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=50)
ac50 <- crop(ac50, e)
ac70 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70)
ac70 <- crop(ac70, e)

b <- setValues(ac50[[1]], NA)
bio.stack <- stack( mget( rep( "b" , 8 ) ) )
bio.stack[[1]] <- now[[1]]
bio.stack[[6]] <- ac50[[1]]
bio.stack[[8]] <- ac70[[1]]
bio.stack <- approxNA(bio.stack)
bio.stack <- (bio.stack/10)

brks <- round(seq(floor(cellStats(bio.stack[[1]], stat = "min", na.rm = TRUE)), ceiling(cellStats(bio.stack[[8]], stat = "max", na.rm = TRUE)), length.out = 10), 0)
nb <- length(brks)-1 
colors <- rev(heat.colors(nb))
years <- as.character(seq(2000, 2070, by = 10))


library(animation)
saveGIF(for(i in 1:8){plot(bio.stack[[i]], col = colors, breaks = brks, main = paste("Mean temp", years[i]))}, movie.name = "Mean_temp.gif", img.name = "Rplot", convert = "convert", clean = TRUE)

#######