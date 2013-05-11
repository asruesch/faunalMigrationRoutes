library(raster)
library(rgdal)
crs = CRS("+proj=eqc +lat_ts=30 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs")
# Data objects
objectsDir ="/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects"
outDir="/media/aaronRueschBackup/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps/rasters"
rasterOptions(overwrite=T)

makeTifs <- function(prefix, outDir, crs, direction, magnitude, richness, variance) {
  ext = extent(min(x) - 25000, max(x) + 25000, min(y) - 25000, max(y) + 25000)
  template = raster(ext=ext, nrows=269, ncols=283, crs=crs)
  variance[is.nan(variance)] = NA
  magnitude[is.na(variance)] = NA
  direction[is.na(magnitude)] = NA
  richness[is.na(magnitude)] = NA
  direction = matrix(direction, nrow=269, ncol=283, byrow=T)
  magnitude = matrix(magnitude, nrow=269, ncol=283, byrow=T)
  richness = matrix(richness, nrow=269, ncol=283, byrow=T)
  variance = matrix(variance, nrow=269, ncol=283, byrow=T)
  direction = raster(direction, template=template)
  magnitude = raster(magnitude, template=template)
  richness = raster(richness, template=template)
  variance = raster(variance, template=template)
  outDirectionFile = paste(outDir, "/", prefix, "_direction.tif", sep="")
  outMagnitudeFile = paste(outDir, "/", prefix, "_magnitude.tif", sep="")
  outRichnessFile = paste(outDir, "/", prefix, "_richness.tif", sep="")
  outVarianceFile = paste(outDir, "/", prefix, "_variance.tif", sep="")
  writeRaster(direction,outDirectionFile,format="ascii",dataType="FLT4S")
  writeRaster(magnitude,outMagnitudeFile,format="ascii",dataType="FLT4S")
  writeRaster(richness,outRichnessFile,format="ascii",dataType="FLT4S")
  writeRaster(variance,outVarianceFile,format="ascii",dataType="FLT4S")
}

amphibsDta = paste(objectsDir, "/amphibs_hiie2.RData", sep="/")
load(amphibsDta)
makeTifs(prefix="amphibs", outDir, crs, direction, magnitude, richness, variance)

mammalsDta = paste(objectsDir, "/mammals_hiie2.RData", sep="/")
load(mammalsDta)
makeTifs(prefix="mammals", outDir, crs, direction, magnitude, richness, variance)
  
birdsDta = paste(objectsDir, "/birds_hiie2.RData", sep="/")
load(birdsDta)
makeTifs(prefix="birds", outDir, crs, direction, magnitude, richness, variance)

conOnlyDta = paste(objectsDir, "/contract_hiie2.RData", sep="/")
load(conOnlyDta)
makeTifs(prefix="contract", outDir, crs, direction, magnitude, richness, variance)

expOnlyDta = paste(objectsDir, "/expand_hiie2.RData", sep="/")
load(expOnlyDta)
makeTifs(prefix="expand", outDir, crs, direction, magnitude, richness, variance)
  
hiiDta = paste(objectsDir, "/hii.RData", sep="/")
load(hiiDta)
makeTifs(prefix="hii", outDir, crs, direction, magnitude, richness, variance)

hiie2Dta = paste(objectsDir, "/hiie2.RData", sep="/")
load(hiie2Dta)
makeTifs(prefix="hiie2", outDir, crs, direction, magnitude, richness, variance)
  
nocostDta = paste(objectsDir, "/nocost.RData", sep="/")
load(nocostDta)
makeTifs(prefix="nocost", outDir, crs, direction, magnitude, richness, variance)
  
  
  
  