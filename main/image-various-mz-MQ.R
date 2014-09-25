#this R script displays various m/z traces with colours

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

# Image R

specno <-1

centermz<-306.4 #in m/z, mass trace of interest
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

R<-intmatrix

col.R <- colorRampPalette(c('white','red'))(100) 
tiff(filename="R.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(R,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.R))
dev.off()

# Image G

specno <-1

centermz<-63.1 #in m/z, mass trace of interest
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

G<-intmatrix

col.G <- colorRampPalette(c('white','green'))(100) 
tiff(filename="G.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(G,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.G))
dev.off()

# Image B

specno <-1

centermz<-62.1 #in m/z, mass trace of interest
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

B<-intmatrix

col.B <- colorRampPalette(c('white','blue'))(100) 
tiff(filename="B.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(B,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.B))
dev.off()

#col.D <- colorRampPalette(c('white','blue'))(100) 
#D<-print(levelplot(intmatrix,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.D))

R<-t(R)
G<-t(G)
B<-t(B)

#define cut-offs

minR<-0.1*max(R)
minG<-0.1*max(G)
minB<-0.1*max(B)

R[ R < minR ] <- NA
G[ G < minG ] <- NA
B[ B < minB ] <- NA

png(filename="RGB-plot.png")
image(G, col=colorRampPalette(c('white','darkgreen'))(15)) 
image(R, add=TRUE, col=colorRampPalette(c('white','darkred'))(15)) 
image(B, add=TRUE, col=colorRampPalette(c('white','darkblue'))(15))
dev.off()
