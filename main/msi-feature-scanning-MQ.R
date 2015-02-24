# This script generates mass images for a selected m/z range and m/z tolerance

# Clear Workspace
rm(list=ls())
#set working directory to source file
setwd("~/SRC/MSI.R/main")

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


centermz<-0.1 #in m/z, start of scanning
scanmaxmz<-1000 #in m/z, end of scanning
scansteps<-0.2 #in m/z, should be smaller than scantolerance
scantolerance<-0.4 #in m/z

while (centermz<=scanmaxmz)
{
  # This will be name and title of the .png file
  # To sort the files according to the m/z mass, a prefix is used
  # The actual pngfiletitle is printed on the console for monitoring the progress of the script
  prefix<-""
  if (centermz<1000000){prefix<-"0"}
  if (centermz<100000){prefix<-"00"}
  if (centermz<10000){prefix<-"000"}
  if (centermz<1000){prefix<-"0000"}
  if (centermz<100){prefix<-"00000"}
  if (centermz<10){prefix<-"000000"}
  pngfilename<-paste(prefix,toString(centermz),"_",toString(scantolerance),".png",sep="")
  pngfiletitle<-paste("m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
  print(pngfiletitle)

specno <-1

while (specno<=elementsimzML)
{
  
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  specmatrix <- cbind(mzs, counts)
  
  lowermass<-centermz-scantolerance
  highermass<-centermz+scantolerance
  
  mzrangeindex <- (specmatrix[,1] > lowermass & specmatrix[,1] < highermass)
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

png(pngfilename)
#if you don't want contour lines, change to contour=FALSE
#rainbow colors
print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = rainbow(100,start=1/5)))
#terrain (map-like) colors
#print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
#greyscale
#print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = grey.colors(100)))
dev.off()

print("mz image saved to file.")
centermz=centermz+scansteps
}
