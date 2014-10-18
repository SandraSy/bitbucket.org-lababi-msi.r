# This script generates an .tiff image of a selected m/z range
# as well as an image and a .mzML file for the spectrum with the 
# highest intensity for this mz.
# Please enter the desired mz in centermz and the tolerance

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

specno <-1

centermz<-137.1325 #in m/z, mass trace of interest
scantolerance<-0.4 #in m/z

lowermass<-centermz-scantolerance
highermass<-centermz+scantolerance

while (specno<=elementsimzML)
{

mzs <- mass(imagespectra[[specno]])
counts <- intensity(imagespectra[[specno]])
specmatrix <- cbind(mzs, counts)

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


# printing spectrum. format options could be e.g. png, eps, tiff, pdf. You can control the resolution with res
prefix<-""
if (centermz<1000000){prefix<-"0"}
if (centermz<100000){prefix<-"00"}
if (centermz<10000){prefix<-"000"}
if (centermz<1000){prefix<-"0000"}
if (centermz<100){prefix<-"00000"}
if (centermz<10){prefix<-"000000"}
tifffilename<-paste(prefix,toString(centermz),"_",toString(scantolerance),".tiff",sep="")
tifffiletitle<-paste("m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
tiff(filename=tifffilename,res=1200,compression="lzw",height=200,width=200,units="mm")
#if you don't want contour lines, change to contour=FALSE
#rainbow colors
print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = rainbow(100,start=1/5)))
#terrain (map-like) colors
#print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
#greyscale
#print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = grey.colors(100)))
dev.off()

# search for spectrum with highest intensity for this m/z
maxintind = which(intmatrix == max(intmatrix), arr.ind = TRUE)
#max(intmatrix)
ymaxint=rownames(intmatrix)[maxintind[,1]]
xmaxint=colnames(intmatrix)[maxintind[,2]]
specposition['y']<-as.numeric(ymaxint)
specposition['x']<-as.numeric(xmaxint)
xs<-imagepos[1,]==specposition['x']
ys<-imagepos[2,]==specposition['y']
specno<-which(xs&ys)

#extract spectrum data
mzs <- mass(imagespectra[[specno]])
counts <- intensity(imagespectra[[specno]])
#plot to tiff
# OUTPUT
tiffspecname<-paste(toString(specposition['x']),"x_",toString(specposition['y']),"y_maxInt_",prefix,toString(centermz),"_",toString(scantolerance),".tiff",sep="")
tiffspectitle<-paste("x=",toString(specposition['x']),", y=",toString(specposition['y'])," m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
mzMLspecname<-paste(toString(specposition['x']),"x_",toString(specposition['y']),"y_maxInt_",prefix,toString(centermz),"_",toString(scantolerance),".mzML",sep="")

tiff(filename=tiffspecname,res=1200,compression="lzw",height=200,width=200,units="mm")
plot(mzs,counts,xlab="m/z",ylab="intensity",main=tiffspectitle,"h")
dev.off()
#export .mzML data
s <- list(createMassSpectrum(mass=mzs, intensity=counts))
exportMzMl(s[[1]], force=TRUE, file=mzMLspecname)

print("mz image and spectra saved to files.")
