source("H:/docs/UW_Data/direction_spp_movement/circuitscape/direction_mapping_functions032211_circuitscape.r")
library(maptools)
library(raster)

load("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hiie2.RData")
cellSize <<- 50000

scaleFactor = 0.75
arrowWidthScalar = 3
arrowHeadSize = 0.17
lenBins=c(0,5,11,18,34,59,100)

imageDir <- "H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/fancy/insets"
imageFile <- paste(imageDir, "/hiie2_fancyInset_south.png", sep="")
png(file=imageFile, bg="transparent", width=4000, height=4000)
xlimsouth <- c(-6800000,-5500000); ylimsouth <- c(-4400000,-3100000)
breaks = mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=scaleFactor
	, arrowWidthScalar=arrowWidthScalar
	, arrowHeadSize=arrowHeadSize
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
	, xlim=c(xlimsouth[1]-50000*10, xlimsouth[2]+50000*10)
	, ylim=c(ylimsouth[1]-50000*10, ylimsouth[2]+50000*10)
	
)
dev.off()

a = raster(nrows=4000, ncols=4000, xmn=xlimsouth[1]-50000*10, xmx=xlimsouth[2]+50000*10, ymn=ylimsouth[1]-50000*10, ymx=ylimsouth[2]+50000*10)
a = raster(imageFile)
extent(a) = extent(xmn=xlimsouth[1]-50000*10, xmx=xlimsouth[2]+50000*10, ymn=ylimsouth[1]-50000*10, ymx=ylimsouth[2]+50000*10)



imageDir <- "H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/fancy/insets"
imageFile <- paste(imageDir, "/hiie2_fancyInset_north.png", sep="")
png(file=imageFile, bg="transparent", width=4000, height=4000)
xlimnorth <- c(-8900000, -7600000);  ylimnorth <- c(3700000, 5000000)
breaks = mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=scaleFactor
	, arrowWidthScalar=arrowWidthScalar
	, arrowHeadSize=arrowHeadSize
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
	, xlim=c(xlimnorth[1]-50000*10, xlimnorth[2]+50000*10)
	, ylim=c(ylimnorth[1]-50000*10, ylimnorth[2]+50000*10)
	
)
dev.off()
