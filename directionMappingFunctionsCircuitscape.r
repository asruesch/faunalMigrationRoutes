#Post-processing script for calculating circular statistics for all species range shifts
#Inputs: One ESRI ASCII data for each species denoting the direction of expected movement
#Written by:  Aaron Ruesch, asruesch@gmail.com
#Date: 7/6/2010

#Description: Populates set of directions for each species movement within each cell. Calculates circular statistics on set of directions including mean, variance
#initializeVariables <- function(animal="all") {

library(CircStats)
library(maptools)
library(RColorBrewer)
library(TeachingDemos)
library(classInt)

loadData <- function(
	gcms = c("cccma_cgcm3.1_t47_run1", "cnrm-cm3_run1", "gfdl-cm2.0_run1"
	, "gfdl-cm2.1_run1", "giss-er_run1", "inm-cm3.0_run1", "miroc3.2_medres_run1"
	, "mri_cgcm2.3.2a_run1", "ncar_ccsm3.0_run1", "ukmo-hadcm3_run1")
	, animal = c("birds", "amphibs", "mammals")
	, shiftType = c("contract", "expand")
	, resistance = "hii"
	, root = "Z:/circuitscape/out"
  	, rangeRoot = "Z:/asciigrids"
  	, latlonDir = "D:/Data/ruesch/sppMovement"
  	, vectorObjectFile = "Z:/circuitscape/vectorObjects/mammals.RData") 
	
	{
	scenario <- "a2"

	ncols <- 283
	nrows <- 269
	npixels <- ncols * nrows
	cellSize <- 50000
	
	magnitudeByGcm <- NULL
	directionByGcm <- NULL
	varianceByGcm <- NULL
	richnessByGcm <- NULL
	currentRichnessRatioByGcm <- NULL
	movementRichnessRatioByGcm <- NULL
	
	#need to iterate this next block for all gcms, then resummarize across gcms at end
	for (gcm in gcms) {
		angleMatrix <- NULL
    		currentMatrix <- NULL
    		presenceMatrix <- NULL
    		for (animalClass in animal) {
			print(paste(gcm, "  :::   ", animalClass))
			decade <- list.files(paste(root, "/", resistance, "/", animalClass, "/", gcm, "/", scenario, "/", sep=""))[1]
			directionDirectory <- paste(root, "/", resistance, "/", animalClass, "/", gcm, "/", scenario, "/", decade, "/", sep="")
			print(directionDirectory)
      			decade <- list.files(paste(rangeRoot, "/", animalClass, "/", gcm, "/", scenario, "/", sep=""))[3]
      			rangeDirectory <- paste(rangeRoot, "/", animalClass, "/", gcm, "/", scenario, "/", decade, "/", sep="")
			setwd(directionDirectory)
			if (length(shiftType) == 2) {
        			anglePattern = "(contract|expand)*angle\\.asc$"
        			currentPattern = "(contract|expand)*magnitude\\.asc$"
      			} else {
        			anglePattern = paste(shiftType, "_angle\\.asc$", sep="")
			        currentPattern = paste(shiftType, "_magnitude\\.asc$", sep="")
      			}
      			rangePattern = "[ABM][1-9]*\\.asc$"
			angleFileList = list.files(pattern = anglePattern)
			currentFileList = list.files(pattern = currentPattern)
			rangeFileList = list.files(rangeDirectory, pattern=rangePattern, full.names=TRUE)
			###### TESTING ###########
			#	angleFileList = angleFileList[1:10]
			#	currentFileList = currentFileList[1:10]
			#	rangeFileList = rangeFileList[1:10]
			##########################
			#Start Loop -- This loop populates a list where each pixel has a vector of directions associated with it
			for (angleFile in angleFileList) {
				print(paste(animalClass, "   :::   ", gcm, "   :::   ", scenario, "   :::   ", decade, "   :::   ", angleFile, sep=""))
				#read in ascii file, integerize, format NoData, and convert to degrees, 0 thru 360, starting at East
				asciiData <- scan(file=angleFile, what="numeric", skip = 6, quiet=TRUE, na.strings = "-9999.000000")
				angleRadians <- as.numeric(asciiData)
        			angleMatrix <- cbind(angleMatrix, angleRadians)
			}
			for (currentFile in currentFileList) {
				print(paste(animalClass, "   :::   ", gcm, "   :::   ", scenario, "   :::   ", decade, "   :::   ", currentFile, sep=""))
				#read in ascii file, integerize, format NoData, and convert to degrees, 0 thru 360, starting at East
				asciiData <- scan(file=currentFile, what="numeric", skip = 6, quiet=TRUE, na.strings = "-9999.000000")
				current <- as.numeric(asciiData)
				currentMatrix <- cbind(currentMatrix, current)
			}
			for (rangeFile in rangeFileList) {
				print(paste(animalClass, "   :::   ", gcm, "   :::   ", scenario, "   :::   ", decade, "   :::   ", rangeFile, sep=""))
				asciiData <- scan(file=rangeFile, what="integer", skip = 6, quiet=TRUE, na.strings = "-9999")
				presence = as.numeric(asciiData)
				presence[!(presence == 2 | presence == 4)] = 0
				presence[!(presence == 0)] = 1
				presenceMatrix = cbind(presenceMatrix, presence)
			}
		}
		allNA = apply(angleMatrix, 1, function(x) {all(is.na(x))})
		o = rowSums(sin(angleMatrix) * currentMatrix, na.rm=TRUE)
		a = rowSums(cos(angleMatrix) * currentMatrix, na.rm=TRUE)
		o[allNA] = NA
		a[allNA] = NA
		direction <- atan2(o,a)
		magnitude <- sqrt(o^2 + a^2)
		variance <- apply(angleMatrix, 1, function(x) { circ.disp(x[!is.na(x)])$var })
		currentSpeciesRichness <- rowSums(presenceMatrix, na.rm=TRUE)
		movementRichnessRatio <- magnitude / currentSpeciesRichness
		magnitudeByGcm <- cbind(magnitudeByGcm, magnitude)
		directionByGcm <- cbind(directionByGcm, direction)
		varianceByGcm <- cbind(varianceByGcm, variance)
		richnessByGcm <- cbind(richnessByGcm, currentSpeciesRichness)
		movementRichnessRatioByGcm <- cbind(movementRichnessRatioByGcm, movementRichnessRatio)
	}
	# END GCM LOOP
	#circular statistics across GCMs
	allNA - apply(directionByGcm, 1, function(x) {all(is.na(x))})
	o = rowSums(sin(directionByGcm) * magnitudeByGcm, na.rm=TRUE)
	a = rowSums(cos(directionByGcm) * magnitudeByGcm, na.rm=TRUE)
	o[allNA] = NA
	a[allNA] = NA
	directionMA = atan2(o,a)
	magnitudeMA = sqrt(o^2 + a^2)
	varianceMA = rowMeans(varianceByGcm)
	richnessMA = rowMeans(richnessByGcm)
	richnessMA[richnessMA == 0] = NA
	movementRichnessRatioMA = magnitudeMA / richnessMA

	magnitude <- magnitudeMA
	direction <- directionMA
	variance <- varianceMA
	richness <- richnessMA
	movementRichnessRatio <- movementRichnessRatioMA

	y <- scan(file=paste(latlonDir, "lat.asc", sep="/"), what="integer", skip = 6, quiet=TRUE, na.strings = "-9999")
	y <- as.integer(y) #y[movementIndeces])
	x <- scan(file=paste(latlonDir, "lon.asc", sep="/"), what="integer", skip = 6, quiet=TRUE, na.strings = "-9999")
	x <- as.integer(x) #x[movementIndeces])
	save(x,y,magnitude,direction,variance,richness,movementRichnessRatio,file=vectorObjectFile)
}


##  MAPPING FUNCTIONS

#mapVectors is a function that draws vectors according to an angle
#a variable defining the length of the vector, and a variable defining
#the color of the variable.

#mapVectors is a function that draws vectors according to an angle, a variable defining the length of the vector, and a variable defining the color of the variable.
mapVectorsScaleArrowWidth_Classes <- function
	(angle=direction
	, lenVar=magnitude
	, colorVar=1-variance
	, scaleFactor=3
	, arrowHeadSize=0.25
	, arrowWidthScalar=5
	, colorBrewerScheme="BuPu"
	, xlim = c(-16025344, -3375344)
	, ylim = c(-5922616, 6927384)
	, autoBins=T
	, colorBins=seq(0,1,by=1/6)
	, lenBins=seq(0,20,by=10/6)
	, boundaries = "H:/docs/UW_Data/direction_spp_movement/country_line.shp"
	) {
	print("begin mapping")

	colorVector <- rep(NA, times=length(angle))
	lenVector <- rep(NA, times=length(angle))
	lwdVector <- rep(NA, times=length(angle))
	myPallette <- brewer.pal(6, colorBrewerScheme)
	#bin color variable into quantiles  
	if (autoBins) {
		colClass <- classIntervals(colorVar, n=6, style="equal", dataPrecision=2)
		colorVector <- findColours(colClass, myPallette)
		lenClass <- classIntervals(lenVar, n=6, style="fisher", dataPrecision=0)
		lenVector <- findInterval(lenVar, lenClass$brks)
		breaks = list(sort(unique(colClass$brks)), sort(unique(lenClass$brks)))
	} else {
		colorBins[7] = max(colorVar, na.rm=T)
		colClass = classIntervals(colorVar, n=6, style="fixed", fixedBreaks=colorBins)
		colorVector = findColours(colClass, myPallette)
		lenBins[7] = max(lenVar, na.rm=T)
		lenClass = classIntervals(lenVar, n=6, style="fixed", fixedBreaks=lenBins)
		lenVector = findInterval(lenVar, lenClass$brks)
		breaks = list(sort(unique(colClass$brks)), sort(unique(lenClass$brks)))
	}
	cellSizeScaleFactor <- sqrt((cellSize/2)^2 + (cellSize/2)^2)
	x0 <- as.double(x); y0 <- as.double(y)
	x1 <- (cos(direction) * lenVector * cellSizeScaleFactor*scaleFactor) + x0
	y1 <- (sin(direction) * lenVector * cellSizeScaleFactor*scaleFactor) + y0
	#sort colors so darkest colors draw last
	arrowTable <- data.frame(colorVar=colorVar, colorVector=colorVector
		, lenVector=lenVector, x0=x0, y0=y0, x1=x1, y1=y1)
	arrowTableSort <- arrowTable[order(colorVar),]
	colorVector <- as.character(arrowTableSort$colorVector)
	lenVector <- arrowTableSort$lenVector
	x0 <- arrowTableSort$x0
	x1 <- arrowTableSort$x1
	y0 <- arrowTableSort$y0
	y1 <- arrowTableSort$y1 
	#draw basemap
	print("reading feature class")
	country.map <- readShapeSpatial(boundaries)
	print("feature class loaded")
	# png(file=fileName, bg="transparent", width=4000, height=4000, restoreConsole=TRUE)
	plot(country.map, col='black', xlim=xlim, ylim=ylim)
	#draw arrows
	for (i in 1:length(angle)) {
		if (!any(is.na(x0[i]), is.na(y0[i]), is.na(x1[i]), is.na(y1[i]))) {
			arrows(x0[i], y0[i], x1[i], y1[i], col=colorVector[i], lwd=lenVector[i]*arrowWidthScalar
				, length=arrowHeadSize*lenVector[i])
		}
	}
	return(breaks)
}

