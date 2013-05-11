library(maptools)
library(CircStats)
library(RColorBrewer)


diffAscii <- readAsciiGrid("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/differenceMaps/diff_hiie2.asc")
diff <- diffAscii@data[[1]]
diffCode <- rep(NA, length(diff))

diffCode[which(diff < -1.5)] <- 1
diffCode[which(diff > -1.499999 & diff < -1)] <- 2
diffCode[which(diff > -0.999999 & diff < -0.5)] <- 3
diffCode[which(diff > -0.499999 & diff < 0.5)] <- 4
diffCode[which(diff > 0.500001 & diff < 1)] <- 5
diffCode[which(diff > 1.00001 & diff < 1.5)] <- 6
diffCode[which(diff > 1.500001)] <- 7

diffCodeAscii <- diffAscii
diffCodeAscii@data[[1]] <- diffCode

col <- c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#67A9CF", "#2166AC")


codeRaster <- raster(diffCodeAscii)
country.map <- readShapeSpatial("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/country_line.shp")
png("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/differenceMaps/redToBlue.png", width=5000, height=5000)
plot(codeRaster, col=col, axes=FALSE, legend=FALSE)
plot(country.map, col="black", lwd=6, add=TRUE)
dev.off()
