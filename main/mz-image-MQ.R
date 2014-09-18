# This script generates an image of a selected m/z range
# Please enter the desired mass range in  the line mzrangeindex <- (specmatrix[,1] > 820 & specmatrix[,1] < 840)
# Further important parameters are:
# * Location of the imzMLConverter (.jaddClassPath)
# * Location of imzML file

# Clear Workspace
rm(list=ls())
#set working directory to source file
setwd("~/SRC/MSI.R/main")

# Load Java
#install.packages("rJava")
#library(rJava)
#.jinit()

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")

# Library for 2D maps
#install.packages("latticeExtra")
library(latticeExtra)

# Load the imzML data
#filename <- "/home/rob/SRC/MSI.R/SampleData/test.imzML"
filename <- "/home/rob/MS-DATA/LTP-MSI-chilli/ltpmsi-chilli.imzML"

imagespectra <- importImzMl(filename, centroided=TRUE)
#plotImsSlice(imagespectra)
elementsimzML<-length(imagespectra)

#determine size and position of image
imagepos <- sapply(imagespectra, function(x)metaData(x)$imaging$pos)
minx<-min(imagepos[1,])
maxx<-max(imagepos[1,])
miny<-min(imagepos[2,])
maxy<-max(imagepos[2,])

xrange<-maxx-minx+1
yrange<-maxy-miny+1

intmatrix <- matrix(0,yrange,xrange)
colnames(intmatrix) <- seq(minx,maxx)
rownames(intmatrix) <- seq(miny,maxy)

specno <-1

while (specno<=elementsimzML)
{

mzs <- mass(imagespectra[[specno]])
counts <- intensity(imagespectra[[specno]])
specmatrix <- cbind(mzs, counts)

mzrangeindex <- (specmatrix[,1] > 305 & specmatrix[,1] < 310)
decisionvalue<- sum(mzrangeindex)
if (decisionvalue=="0"){mzpixintensity<-0}
if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-mzrangeints['counts']}
if (decisionvalue>"1"){
mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-sum(mzrangeints[,2])}

specposition<-metaData(imagespectra[[specno]])$imaging$pos
specpositionx<-paste(specposition[1])
specpositiony<-paste(specposition[2])

intmatrix[specpositiony,specpositionx]<-mzpixintensity
print(specno)
specno=specno+1
}

png("mzimage.png")
print(levelplot(intmatrix,main="Intensity m/z",xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
dev.off()

print("mz image saved to file.")
