library(maptools)
library(CircStats)
library(RColorBrewer)

percNaturalAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/percConvVsPercMoving/percent_natural.asc"
hiiAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/resistanceData/hiie2.asc"
countryMapFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/country_line.shp"
percNaturalAscii <- readAsciiGrid(percNaturalAsciiFile)
hiiAscii <- readAsciiGrid(hiiAsciiFile)
country.map <- readShapeSpatial(countryMapFile)
percNatural <- percNaturalAscii@data[[1]]
hii <- hiiAscii@data[[1]]
percConverted <- 1 - percNatural

load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hiie2.RData")
percMoving <- magnitude / richness
percMoving[which(is.nan(percMoving))] <- NA
percMoving[which(percMoving == Inf)] <- NA


####################################
####
####    Quantiles

percMovingBreaks <- as.numeric(quantile(percMoving, probs=seq(0,1,1/3), na.rm=TRUE)[2:3])
percConvertedBreaks <- as.numeric(quantile(percConverted, probs=seq(0,1,1/3), na.rm=TRUE)[2:3])

percMovingTerc <- rep(NA, length(percMoving))
percConvertedTerc <- rep(NA, length(percConverted))

percMovingTerc[which(percMoving <= percMovingBreaks[1])] <- 1  # 1st percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[1] & percMoving <= percMovingBreaks[2])] <- 2 # 2nd percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[2])] <- 3  # 3rd percent moving tercile

percConvertedTerc[which(percConverted <= percConvertedBreaks[1])] <- 1  # 1st percent moving tercile
percConvertedTerc[which(percConverted > percConvertedBreaks[1] & percConverted <= percConvertedBreaks[2])] <- 2 # 2nd percent moving tercile
percConvertedTerc[which(percConverted > percConvertedBreaks[2])] <- 3  # 3rd percent moving tercile

doubleTercile <- rep(NA, length(percConverted))  # Coded as below

# ------------
#  4 | 3 | 2 |   ^
# ------------   |
#  5 | 0 | 1 |   |  % Moving
# ------------   |
#  6 | 7 | 8 |   |
# ------------
#   ------>
#  % Converted

#  color codes

doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 1)] <- 6
doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 2)] <- 7
doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 3)] <- 8
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 1)] <- 5
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 2)] <- 0
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 3)] <- 1
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 1)] <- 4
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 2)] <- 3
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 3)] <- 2

llCol <- brewer.pal(3, "PuOr")[2] # almost white
percConvertedCol <- brewer.pal(7, "PuOr")[5:7] # light orange to dark orange
percMovingCol <- brewer.pal(7, "PuOr")[3:1] # light purple to dark purple
cenCol <- "#C59982"
urCol <- "#4D0000"

col <- c(cenCol, percConvertedCol[3], urCol, percMovingCol[3:1], llCol, percConvertedCol[1:2])

codeAscii <- percNaturalAscii
codeAscii@data[[1]] <- doubleTercile

codeRaster <- raster(codeAscii)
png("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/percConvVsPercMoving/quantiles.png", width=5000, height=5000)
# image(codeAscii, col=col)
# plot(codeAscii, col=col)
plot(codeRaster, col=col, axes=FALSE, legend=FALSE)
plot(country.map, col="black", lwd=6, add=TRUE)
dev.off()

####################################
####
####    Percent converted on 20th and 40th percentiles

percMovingBreaks <- as.numeric(quantile(percMoving, probs=seq(0,1,1/3), na.rm=TRUE)[2:3])
percConvertedBreaks <- as.numeric(quantile(percConverted, probs=c(0,0.2,0.4), na.rm=TRUE)[2:3])

percMovingTerc <- rep(NA, length(percMoving))
percConvertedTerc <- rep(NA, length(percConverted))

percMovingTerc[which(percMoving <= percMovingBreaks[1])] <- 1  # 1st percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[1] & percMoving <= percMovingBreaks[2])] <- 2 # 2nd percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[2])] <- 3  # 3rd percent moving tercile

percConvertedTerc[which(percConverted <= percConvertedBreaks[1])] <- 1  # 1st percent moving tercile
percConvertedTerc[which(percConverted > percConvertedBreaks[1] & percConverted <= percConvertedBreaks[2])] <- 2 # 2nd percent moving tercile
percConvertedTerc[which(percConverted > percConvertedBreaks[2])] <- 3  # 3rd percent moving tercile

doubleTercile <- rep(NA, length(percConverted))  # Coded as below


doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 1)] <- 6
doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 2)] <- 7
doubleTercile[which(percMovingTerc == 1 & percConvertedTerc == 3)] <- 8
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 1)] <- 5
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 2)] <- 0
doubleTercile[which(percMovingTerc == 2 & percConvertedTerc == 3)] <- 1
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 1)] <- 4
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 2)] <- 3
doubleTercile[which(percMovingTerc == 3 & percConvertedTerc == 3)] <- 2

llCol <- brewer.pal(3, "PuOr")[2] # almost white
percConvertedCol <- brewer.pal(7, "PuOr")[5:7] # light orange to dark orange
percMovingCol <- brewer.pal(7, "PuOr")[3:1] # light purple to dark purple
cenCol <- "#C59982"
urCol <- "#4D0000"

col <- c(cenCol, percConvertedCol[3], urCol, percMovingCol[3:1], llCol, percConvertedCol[1:2])

codeAscii <- percNaturalAscii
codeAscii@data[[1]] <- doubleTercile

codeRaster <- raster(codeAscii)
png("D:/Data/ruesch/sppMovement/percConvVsPercMoving/20_40_converted.png", width=5000, height=5000, restoreConsole=TRUE)
# image(codeAscii, col=col)
# plot(codeAscii, col=col)
plot(codeRaster, col=col, axes=FALSE, legend=FALSE)
plot(country.map, col="black", lwd=6, add=TRUE)
dev.off()


####################################
####
####    Quantiles based on HII

percMovingBreaks <- round(as.numeric(quantile(percMoving, probs=seq(0,1,1/3), na.rm=TRUE)[2:3]))
hiiBreaks <- as.numeric(quantile(hii, probs=seq(0,1,1/3), na.rm=TRUE)[2:3])

percMovingTerc <- rep(NA, length(percMoving))
hiiTerc <- rep(NA, length(hii))

percMovingTerc[which(percMoving <= percMovingBreaks[1])] <- 1  # 1st percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[1] & percMoving <= percMovingBreaks[2])] <- 2 # 2nd percent moving tercile
percMovingTerc[which(percMoving > percMovingBreaks[2])] <- 3  # 3rd percent moving tercile

hiiTerc[which(hii <= hiiBreaks[1])] <- 1  # 1st percent moving tercile
hiiTerc[which(hii > hiiBreaks[1] & hii <= hiiBreaks[2])] <- 2 # 2nd percent moving tercile
hiiTerc[which(hii > hiiBreaks[2])] <- 3  # 3rd percent moving tercile

doubleTercile <- rep(NA, length(hii))  # Coded as below


doubleTercile[which(percMovingTerc == 1 & hiiTerc == 1)] <- 6
doubleTercile[which(percMovingTerc == 1 & hiiTerc == 2)] <- 7
doubleTercile[which(percMovingTerc == 1 & hiiTerc == 3)] <- 8
doubleTercile[which(percMovingTerc == 2 & hiiTerc == 1)] <- 5
doubleTercile[which(percMovingTerc == 2 & hiiTerc == 2)] <- 0
doubleTercile[which(percMovingTerc == 2 & hiiTerc == 3)] <- 1
doubleTercile[which(percMovingTerc == 3 & hiiTerc == 1)] <- 4
doubleTercile[which(percMovingTerc == 3 & hiiTerc == 2)] <- 3
doubleTercile[which(percMovingTerc == 3 & hiiTerc == 3)] <- 2

llCol <- brewer.pal(3, "PuOr")[2] # almost white
hiiCol <- brewer.pal(7, "PuOr")[5:7] # light orange to dark orange
percMovingCol <- brewer.pal(7, "PuOr")[3:1] # light purple to dark purple
cenCol <- "#C59982"
urCol <- "#4D0000"

col <- c(cenCol, hiiCol[3], urCol, percMovingCol[3:1], llCol, hiiCol[1:2])

codeAscii <- percNaturalAscii
codeAscii@data[[1]] <- doubleTercile

codeRaster <- raster(codeAscii)
png("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/percConvVsPercMoving/hiie2_quantiles.png", width=5000, height=5000)
# image(codeAscii, col=col)
# plot(codeAscii, col=col)
plot(codeRaster, col=col, axes=FALSE, legend=FALSE)
plot(country.map, col="black", lwd=6, add=TRUE)
dev.off()

