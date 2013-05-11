source("H:/docs/UW_Data/direction_spp_movement/circuitscape/direction_mapping_functions032211_circuitscape.r")
library(maptools)

load("H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorObjects/nocost.RData")
cellSize <<- 50000

allClasses = c(0,5,11,18,34,59,100)
taxaClasses = c(0,5,11,18,34,59,100)
shiftTypeClasses = c(0,2,4,8,16,32,100)

lenBins = allClasses

imageDir <- "H:/docs/UW_Data/direction_spp_movement/circuitscape/vectorMaps"
imageFile <- paste(imageDir, "/nocost.tif", sep="")
tiff(file=imageFile, bg="white", width=4000, height=4000, units="px")
xlimsouth <- c(-6800000,-5500000); ylimsouth <- c(-4400000,-3100000)
xInset <- -14000000; ySInset <- -3200000
xlimnorth <- c(-8900000, -7600000);  ylimnorth <- c(3700000, 5000000)
yNInset <- 1600000
breaks = mapVectorsScaleArrowWidth_Classes(angle=direction
	, lenVar=magnitude / richness
	, colorVar=1-variance
	, scaleFactor=0.6
	, arrowWidthScalar=1
	, arrowHeadSize=0.04
	, colorBrewerScheme="YlOrRd"
	, autoBins=F
	, colorBins=seq(0,1,by=1/6)
	, lenBins=lenBins
)
subplot(
  mapVectorsScaleArrowWidth_Classes(angle=direction
    , lenVar=magnitude/richness
    , colorVar=1-variance
    , scaleFactor=0.75
    , arrowHeadSize=0.17
    , arrowWidthScalar=3
    , colorBrewerScheme="YlOrRd"
    , autoBins=FALSE
    , colorBins=seq(0,1,by=1/6)
    , lenBins=lenBins
    , xlim=c(xlimsouth[1]+50000*3, xlimsouth[2]-50000*3)
    , ylim=c(ylimsouth[1]+50000*3, ylimsouth[2]-50000*3)
  ), x=xInset, y=ySInset,size=c(18,18)
)
subplot(
  mapVectorsScaleArrowWidth_Classes(angle=direction
    , lenVar=magnitude/richness
    , colorVar=1-variance
    , scaleFactor=0.75
    , arrowHeadSize=0.17
    , arrowWidthScalar=3
    , colorBrewerScheme="YlOrRd"
    , autoBins=FALSE
    , colorBins=seq(0,1,by=1/6)
    , lenBins=lenBins
    , xlim=c(xlimnorth[1]+50000*3, xlimnorth[2]-50000*3)
    , ylim=c(ylimnorth[1]+50000*3, ylimnorth[2]-50000*3)
  ), x=xInset, y=yNInset,size=c(18,18)
)
polygon(c(xInset-2330000, xInset+2330000, xInset+2330000, xInset-2330000)
        , c(ySInset+2330000, ySInset+2330000, ySInset-2330000, ySInset-2330000)
        , border="gray40", lwd=8)
polygon(c(xInset-2330000, xInset+2330000, xInset+2330000, xInset-2330000)
        , c(yNInset+2330000, yNInset+2330000, yNInset-2330000, yNInset-2330000)
        , border="gray40", lwd=8)
polygon(c(xlimsouth[1], xlimsouth[2], xlimsouth[2], xlimsouth[1])
        , c(ylimsouth[2], ylimsouth[2], ylimsouth[1], ylimsouth[1])
        , border="gray40", lwd=8)
polygon(c(xlimnorth[1], xlimnorth[2], xlimnorth[2], xlimnorth[1])
        , c(ylimnorth[2], ylimnorth[2], ylimnorth[1], ylimnorth[1])
        , border="gray40", lwd=8)

xTail <- rep(-9.725^7, times=6)
yTail <- seq(from=-11000, to=-3*10^6, length.out=6)
xHead <- xTail + seq((1/6),1,by=1/6) * 50000 * 5
for (i in 1:6) {
  arrows(xTail[i], yTail[i], xHead[i], yTail[i], col="black"
         , length=(seq((1/6),1,by=1/6)*0.75)[i], lwd=(seq((1/6),1,by=1/6)*20)[i]
  )
}
legendCol <- brewer.pal(6, "YlOrRd")
legendLabs = cbind(round(breaks[[1]][1:5], 2), round(breaks[[1]][2:6], 2))
legendLabs = apply(legendLabs, 1, function(x) {paste(x[1], x[2], sep="-")})
legendColVal = c(legendLabs, paste(round(breaks[[1]][6], 2), "-", round(breaks[[1]][7], 2), sep=""))
legendLabs = cbind(round(breaks[[2]][1:5], 1), round(breaks[[2]][2:6], 1))
legendLabs = apply(legendLabs, 1, function(x) {paste(x[1], x[2], sep="-")})
legendLenVal = c(legendLabs, paste(">", round(breaks[[2]][6], 0), sep=""))

for (i in 1:6) {
  polygon(c(-1.1*10^7-150000, -1.1*10^7+150000, -1.1*10^7+150000, -1.1*10^7-150000)
          , c(yTail[i]+150000, yTail[i]+150000, yTail[i]-150000, yTail[i]-150000)
          , border="black", lwd=8, col=legendCol[7-i])
  text(-1.1*10^7+300000, yTail[i], legendColVal[i], cex=6, pos=4)
  text(-1.1*10^7+1700000, yTail[i], legendLenVal[i], cex=6, pos=4)
}
text(-1.1*10^7-150000, yTail[1]+600000, "Variance in\ndirection", cex=6, pos=4, font=2)
text(-1.1*10^7+1700000, yTail[1]+600000, "Movements\n/ richness", cex=6, pos=4, font=2)

dev.off()
