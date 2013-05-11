
#mapVectors is a function that draws vectors according to an angle, a variable defining the length of the vector, and a variable defining the color of the variable.
mapVectorsScaleArrowWidth_Classes <- function
	(angle=direction
	, lenVar=magnitude
	, scaleFactor=0.4
	, arrowHeadSize=0.05
	, arrowWidthScalar=1.5
	, color="mediumpurple4"
	, lenBins=seq(0,20,by=10/6)
	) {
	print("begin mapping")
	cellSize=50000
	xy = coordinates(angle)
	x = xy[,1]
	y = xy[,2]
	angle = getValues(angle)
	lenVar = getValues(lenVar)
	lenVector <- rep(NA, times=length(angle))
	lwdVector <- rep(NA, times=length(angle))
	lenClass <- classIntervals(lenVar, n=6, style="quantile", dataPrecision=0)
	lenVector <- findInterval(lenVar, lenClass$brks)
	lenVector[lenVar == 0] = 0
	breaks = sort(unique(lenClass$brks))
	cellSizeScaleFactor <- sqrt((cellSize/2)^2 + (cellSize/2)^2)
	x0 <- as.double(x); y0 <- as.double(y)
	x1 <- (cos(angle) * lenVector * cellSizeScaleFactor*scaleFactor) + x0
	y1 <- (sin(angle) * lenVector * cellSizeScaleFactor*scaleFactor) + y0
	for (i in 1:length(angle)) {
		if (!any(is.na(x0[i]), is.na(y0[i]), is.na(x1[i]), is.na(y1[i]))) {
			if (lenVector[i] == 0) {
				points(x0[i], y0[i], col=color, pch=16)
			} else {
				arrows(x0[i], y0[i], x1[i], y1[i], col=color, lwd=lenVector[i]*arrowWidthScalar
					, length=arrowHeadSize*lenVector[i])
			}
		}
	}
	return(breaks)
}
install.packages(c("raster", "RColorBrewer", "classInt"))
library(raster)
library(RColorBrewer)
library(classInt)
pal <- c("#BEBEBE", brewer.pal(7, "Blues")[7:3])

rangeShift = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332.asc")
hiie2 = raster("H:/docs/UW_Data/direction_spp_movement/resistanceData/hiie2.asc")
hiie2Vals = getValues(hiie2)
colClass = classIntervals(hiie2Vals, n=5, style="fisher", dataPrecision=2)
colClass = findInterval(hiie2Vals, colClass$brks)
hiie2 = setValues(hiie2, colClass)
hiie2[hiie2 == 6] = 5
hiie2[rangeShift == 1] = 0

nocostAngleContract = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_contract_nocost_angle.asc")
nocostMagnitudeContract = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_contract_nocost_magnitude.asc")
nocostAngleExpand = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_expand_nocost_angle.asc")
nocostMagnitudeExpand = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_expand_nocost_magnitude.asc")

hiie2AngleContract = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_contract_hiie2_angle.asc")
hiie2MagnitudeContract = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_contract_hiie2_magnitude.asc")
hiie2AngleExpand = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_expand_hiie2_angle.asc")
hiie2MagnitudeExpand = raster("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/data/A1332_expand_hiie2_magnitude.asc")

xlim1 = c(-7.9e6,-3.5e6); ylim1 = c(-6.1e6, -1.7e6)
xlim2 = c(-6.9e6,-4.5e6); ylim2 = c(-4.5e6,-2.1e6)

tiff(file= "H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/schematicFigure/schematic.tif"
	, width=4000
	, height=4000
	, units="px"
)
par(mfrow=c(2,2), mar=c(0,0,0,0))
plot(rangeShift, col=c("grey", "mediumpurple4", "green3", "orange2"), xlim=xlim1, ylim=ylim1, axes=FALSE, legend=FALSE, box=FALSE)
plot(rangeShift, col=c("grey", "mediumpurple4", "green3", "orange2"), xlim=xlim2, ylim=ylim2, axes=FALSE, legend=FALSE, box=FALSE)
nocostAngleContract[rangeShift != 2] = NA
nocostMagnitudeContract[rangeShift != 2] = NA
mapVectorsScaleArrowWidth_Classes(angle=nocostAngleContract
	, lenVar=nocostMagnitudeContract
	, color="white"
)
nocostAngleExpand[rangeShift != 3] = NA
nocostMagnitudeExpand[rangeShift != 3] = NA
mapVectorsScaleArrowWidth_Classes(angle=nocostAngleExpand
	, lenVar=nocostMagnitudeExpand
	, color="white"
)
plot(hiie2, col=pal, xlim=xlim2, ylim=ylim2, axes=FALSE, legend=FALSE, box=FALSE)
plot(hiie2, col=pal, xlim=xlim2, ylim=ylim2, axes=FALSE, legend=FALSE, box=FALSE)
hiie2AngleContract[rangeShift != 2] = NA
hiie2MagnitudeContract[rangeShift != 2] = NA
mapVectorsScaleArrowWidth_Classes(angle=hiie2AngleContract
	, lenVar=hiie2MagnitudeContract
	, color="white"
)
hiie2AngleExpand[rangeShift != 3] = NA
hiie2MagnitudeExpand[rangeShift != 3] = NA
mapVectorsScaleArrowWidth_Classes(angle=hiie2AngleExpand
	, lenVar=hiie2MagnitudeExpand
	, color="white"
)
dev.off()



