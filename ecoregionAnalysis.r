library(maptools)
library(CircStats)

# define names of covariate files
biomeAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/ecoregionAnalysis/NSbiome.asc"
percNaturalAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/percConvVsPercMoving/percent_natural.asc"
elevRangeAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/ecoregionAnalysis/elevation_range.asc"
hiiAsciiFile <- "/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/percConvVsPercMoving/hii.asc"
# read covariate ascii files
biome <- readAsciiGrid(biomeAsciiFile)@data[[1]]
percNatural <- readAsciiGrid(percNaturalAsciiFile)@data[[1]]
elevRange <- readAsciiGrid(elevRangeAsciiFile)@data[[1]]
hii <- readAsciiGrid(hiiAsciiFile)@data[[1]]
# summarize by biome
percNatSumm <- aggregate(percNatural, list(biome=biome), mean, na.rm=TRUE)
elevRangeSumm <- aggregate(elevRange, list(biome=biome), mean, na.rm=TRUE)
hiiSumm <- aggregate(hii, list(biome=biome), mean, na.rm=TRUE)
n <- NULL
for (i in sort(unique(biome))) { 
	n <- c(n, length(which(biome == i)))
}
summTable <- cbind(sort(unique(biome)), n, percNatSumm[,2], elevRangeSumm[,2], hiiSumm[,2])

# richness <- readAsciiGrid("D:/Data/ruesch/sppMovement/richnessMaps/total.asc")@data[[1]]

load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/nocost.RData")
nocostMag <- magnitude/richness
nocostMag[which(is.nan(nocostMag))] <- NA
nocostMag[which(nocostMag == Inf)] <- NA
load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hii.RData")
hiiMag <- magnitude/richness
hiiMag[which(is.nan(hiiMag))] <- NA
hiiMag[which(hiiMag == Inf)] <- NA
diff <- hiiMag - nocostMag
diffSumm <- aggregate(diff, list(biome=biome), mean, na.rm=TRUE)
summTable <- cbind(summTable, diffSumm[,2])

load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/nocost.RData")
percMoving <- magnitude / richness
percMoving[which(is.nan(percMoving))] <- NA
percMoving[which(percMoving == Inf)] <- NA
percMovingSumm <- aggregate(percMoving, list(biome=biome), mean, na.rm=TRUE)
totMovingSumm <- aggregate(magnitude, list(biome=biome), mean, na.rm=TRUE)
direction[is.na(percMoving)] <- NA
dirSumm <- NULL
varSumm <- NULL
for (i in sort(unique(biome))) {
	pInd <- which(biome == i)
	dirBiome <- direction[pInd]
	dirBiome <- dirBiome[!is.na(dirBiome)]
	dirMean <- circ.mean(dirBiome)
	dirSumm <- c(dirSumm, dirMean)
	dirVar <- circ.disp(dirBiome)$var
	varSumm <- c(varSumm, dirVar)
}
summTable <- cbind(summTable, percMovingSumm[,2], totMovingSumm[,2], varSumm, dirSumm)


load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hii.RData")
percMoving <- magnitude / richness
percMoving[which(is.nan(percMoving))] <- NA
percMoving[which(percMoving == Inf)] <- NA
percMovingSumm <- aggregate(percMoving, list(biome=biome), mean, na.rm=TRUE)
totMovingSumm <- aggregate(magnitude, list(biome=biome), mean, na.rm=TRUE)
direction[is.na(percMoving)] <- NA
dirSumm <- NULL
varSumm <- NULL
for (i in sort(unique(biome))) {
	pInd <- which(biome == i)
	dirBiome <- direction[pInd]
	dirBiome <- dirBiome[!is.na(dirBiome)]
	dirMean <- circ.mean(dirBiome)
	dirSumm <- c(dirSumm, dirMean)
	dirVar <- circ.disp(dirBiome)$var
	varSumm <- c(varSumm, dirVar)
}
summTable <- cbind(summTable, percMovingSumm[,2], totMovingSumm[,2], varSumm, dirSumm)

load("/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/hiie2.RData")
percMoving <- magnitude / richness
percMoving[which(is.nan(percMoving))] <- NA
percMoving[which(percMoving == Inf)] <- NA
percMovingSumm <- aggregate(percMoving, list(biome=biome), mean, na.rm=TRUE)
totMovingSumm <- aggregate(magnitude, list(biome=biome), mean, na.rm=TRUE)
direction[is.na(percMoving)] <- NA
dirSumm <- NULL
varSumm <- NULL
for (i in sort(unique(biome))) {
	pInd <- which(biome == i)
	dirBiome <- direction[pInd]
	dirBiome <- dirBiome[!is.na(dirBiome)]
	dirMean <- circ.mean(dirBiome)
	dirSumm <- c(dirSumm, dirMean)
	dirVar <- circ.disp(dirBiome)$var
	varSumm <- c(varSumm, dirVar)
}
summTable <- cbind(summTable, percMovingSumm[,2], totMovingSumm[,2], varSumm, dirSumm)
summTable <- as.data.frame(summTable)

names(summTable) <- c("biome", "nPixels", "percNatural", "elevRange", "hii", "diffPercMoving", "nocostPercMoving", "nocostTotMoving", "nocostVar", "nocostDir", "hiiPercMoving", "hiitotMoving", "hiiVar", "hiiDir", "hiie2PercMoving", "hiie2TotMoving", "hiie2Var", "hiie2Dir")

write.csv(summTable, file="/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/ecoregionAnalysis/ecoregionAnalysis_NS.csv")
