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
library(rJava)
.jinit()

# Library for 2D maps
#install.packages("latticeExtra")
library(latticeExtra)

# You need the imzMLConverter: http://www.cs.bham.ac.uk/~ibs/imzMLConverter/
.jaddClassPath(path="/home/rob/SRC/MSI.R/imzMLConverter/imzMLConverter.jar")

# Load the imzML data
filename <- "/home/rob/SRC/MSI.R/SampleData/test.imzML"
#filename <- "/home/rob/MS-DATA/LTP-MSI-chilli/ltpms-chilli.imzML"
imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(filename)
imagewidth <- J(imzML, 'getWidth')
imageheight <- J(imzML, 'getHeight')

intmatrix <- matrix(0,imagewidth,imageheight)
colnames(intmatrix) <- seq(1:imageheight)

x <- 1
y <- 1
specno <-1

while (y<=imageheight)
{  
x <- 1
while (x<=imagewidth)
{
spectrum <- J(imzML, 'getSpectrum', as.integer(x), as.integer(y))
mzs <- J(spectrum, 'getmzArray')
counts <- J(spectrum, 'getIntensityArray')
specmatrix <- cbind(mzs, counts)

mzrangeindex <- (specmatrix[,1] > 820 & specmatrix[,1] < 840)
decisionvalue<- sum(mzrangeindex)
if (decisionvalue=="0"){mzpixintensity<-0}
if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-mzrangeints['counts']}
if (decisionvalue>"1"){
mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-sum(mzrangeints[,2])}
intmatrix[x,y]<-mzpixintensity
print(specno)
x=x+1
specno=specno+1
}
y=y+1
}

png("mzimage.png")
print(levelplot(intmatrix,main="Intensity m/z",xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
dev.off()

print("mz image saved to file.")
