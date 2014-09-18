# This script generates mass images for a selected m/z range and m/z tolerance
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

centermz<-0 #in m/z, start of scanning
scanmaxmz<-1000 #in m/z, end of scanning
scansteps<-5 #in m/z, should be smaller than scantolerance
scantolerance<-10 #in m/z

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
intmatrix[x,y]<-mzpixintensity
print(specno)
x=x+1
specno=specno+1
}
y=y+1
}

png(pngfilename)
print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
dev.off()

print("mz image saved to file.")
centermz=centermz+scansteps
}
