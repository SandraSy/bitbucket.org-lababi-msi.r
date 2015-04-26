# This script generates mass images for a selected m/z range and m/z tolerance

# Clear Workspace
rm(list=ls())

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")

# Library for 2D maps
#install.packages("latticeExtra")
library(latticeExtra)

# Load the imzML data
filename <- file.path("..", "data", "ltpmsi-chilli.imzML")

imagespectra <- importImzMl(filename, centroided=TRUE)

centermz<-40 #in m/z, start of scanning
scanmaxmz<-540 #in m/z, end of scanning
scansteps<-0.2 #in m/z, should be smaller than scantolerance
scantolerance<-0.4 #in m/z

## trim spectra to m/z range of interest, to reduce data size
imagespectra <- trim(imagespectra, range=c(centermz, scanmaxmz))

## generate intensity matrix
intmatrix <- intensityMatrix(imagespectra)
mzs <- as.double(colnames(intmatrix))

#determine size and position of image
imagepos <- sapply(imagespectra, function(x)metaData(x)$imaging$pos)
ranges <- apply(imagepos, MARGIN=1, FUN=range)
## adjust imagepos and rotate them to use them as index in the loop
imagepos <- t(imagepos - ranges[1, ] + 1)

n <- apply(ranges, MARGIN=2, FUN=diff)+1
curintmat <- matrix(0, nrow=n[1], ncol=n[2])

## precalculate scan ranges
centermz <- seq(from=centermz, to=scanmaxmz, by=scansteps)
lowermass <- centermz-scantolerance
highermass <- centermz+scantolerance

## precalculate indices (columns)
leftIdx <- findInterval(lowermass, mzs, all.inside=TRUE)
rightIdx <- findInterval(highermass, mzs, all.inside=TRUE)

col <- rainbow(100, start=1/5)

## loop through intervals
for (i in seq(along=centermz)) {
  # This will be name and title of the .png file
  # To sort the files according to the m/z mass, a prefix is used
  # The actual pngfiletitle is printed on the console for monitoring the progress of the script
  pngfilename <- sprintf("%08.1f_%.1f.png", centermz[i], scantolerance)
  pngfiletitle <- sprintf("mz %.1f+/-%.1f", centermz[i], scantolerance)

  curintmat[imagepos] <- rowSums(intmatrix[,leftIdx[i]:rightIdx[i]], na.rm=TRUE)

  message(pngfilename)
  png(pngfilename)
  print(levelplot(curintmat,
                  main=pngfiletitle, xlab="x/ pixel", ylab="y/ pixel",
                  scales=list(draw=FALSE), contour=TRUE, pretty=TRUE,
                  col.regions=col))
  dev.off()
}
