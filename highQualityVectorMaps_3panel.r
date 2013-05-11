source("H:/docs/UW_Data/direction_spp_movement/circuitscape/direction_mapping_functions032211_circuitscape.r")

load("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hiie2.RData")
cellSize <<- 50000

scaleFactor = 1.25
arrowWidthScalar = 4
arrowHeadSize = 0.23
lenBins=c(0,5,11,18,34,59,100)

imageDir <- "H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps"
imageFile <- paste(imageDir, "/3panel.tif", sep="")
tiff(file=imageFile, bg="white", width=4000, height=4000, units="px")
par(mfrow=c(2,2))
xlim1 <- c(-9350000,-8050000); ylim1 <- c(1400000,2700000)
xlim2 <- c(-11500000,-10200000); ylim2 <- c(4600000,5900000)
xlim3 <- c(-5100000,-3800000); ylim3 <- c(-3100000,-1800000)
mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=0.5
	, arrowWidthScalar=0.65
	, arrowHeadSize=0.04
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
)
polygon(c(xlim1[1], xlim1[2], xlim1[2], xlim1[1])
	, c(ylim1[2], ylim1[2], ylim1[1], ylim1[1])
	, border="gray40", lwd=8)
polygon(c(xlim2[1], xlim2[2], xlim2[2], xlim2[1])
	, c(ylim2[2], ylim2[2], ylim2[1], ylim2[1])
	, border="gray40", lwd=8)
polygon(c(xlim3[1], xlim3[2], xlim3[2], xlim3[1])
	, c(ylim3[2], ylim3[2], ylim3[1], ylim3[1])
	, border="gray40", lwd=8)
mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=scaleFactor
	, arrowWidthScalar=arrowWidthScalar
	, arrowHeadSize=arrowHeadSize
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
	, xlim=c(xlim1[1], xlim1[2])
	, ylim=c(ylim1[1], ylim1[2])
)
mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=scaleFactor
	, arrowWidthScalar=arrowWidthScalar
	, arrowHeadSize=arrowHeadSize
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
	, xlim=c(xlim2[1], xlim2[2])
	, ylim=c(ylim2[1], ylim2[2])
)
mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=scaleFactor
	, arrowWidthScalar=arrowWidthScalar
	, arrowHeadSize=arrowHeadSize
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
	, xlim=c(xlim3[1], xlim3[2])
	, ylim=c(ylim3[1], ylim3[2])
)
dev.off()








