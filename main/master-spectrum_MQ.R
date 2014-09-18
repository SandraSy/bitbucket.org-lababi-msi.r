# This script summarizes the x (e.g. 100) most intense peaks/ point in one single spectrum
# The increasing number indicates the number of processed spectra
# The final spectrum is printed into a publication quality .tiff file
# Important parameters are:
# * Location of the imzMLConverter (.jaddClassPath)
# * Location of imzML file
# You can change the number of selected peaks with (this would be the setting for 100 signals):
#if (length(newspec)>"200"){        //2x the number of peaks 
#  newspec <- newspec[1:100,]}      //1:x strongest signals


# Clear Workspace
rm(list=ls())
# set working directory to source file
setwd("~/SRC/MSI.R/main")

# Load Java
#install.packages("rJava")
library(rJava)
.jinit()

# You need the imzMLConverter: http://www.cs.bham.ac.uk/~ibs/imzMLConverter/
.jaddClassPath(path="/home/rob/SRC/MSI.R/imzMLConverter/imzMLConverter.jar")

# Load the imzML data
filename <- "/home/rob/SRC/MSI.R/SampleData/test.imzML"
#filename <- "/home/rob/MS-DATA/LTP-MSI-chilli/ltpms-chilli.imzML"
imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(filename)
imagewidth <- J(imzML, 'getWidth')
imageheight <- J(imzML, 'getHeight')

x <- 1
y <- 1
specno <-1

while (y<=imageheight)
#while (y<=5)
    
{  
  x <- 1
while (x<=imagewidth)
#  while (x<=5)
  {
    spectrum <- J(imzML, 'getSpectrum', as.integer(x), as.integer(y))
    mzs <- J(spectrum, 'getmzArray')
    counts <- J(spectrum, 'getIntensityArray')
    newspec  <- cbind(mzs, counts)
    newspec <- newspec[order(-counts),]
    if (length(newspec)>"20"){
    newspec <- newspec[1:10,]}
    
    mzs<-newspec[,1]
    counts<-newspec[,2]
    
    if(exists("mzsum")=="TRUE"){mzsum<-c(mzsum,mzs)}
    if(exists("countsum")=="TRUE"){countsum<-c(countsum,counts)}
    if(exists("mzsum")=="FALSE"){mzsum<-mzs}
    if(exists("countsum")=="FALSE"){countsum<-counts}
    specmatrix <- cbind(mzsum, countsum)  
    print(specno)
  x=x+1
  specno=specno+1
  }
  y=y+1
}

# printing spectrum. format options could be e.g. png, eps, tiff, pdf. You can control the resolution with res
tiff(filename="masterspectrum.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
plot(specmatrix[,1],specmatrix[,2],xlab="m/z",ylab="intensity",main="Pseudo image spectrum","h")
dev.off()

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")
s <- list(createMassSpectrum(mass=specmatrix[,1], intensity=specmatrix[,2]))
## export a single spectrum
exportMzMl(s[[1]], file="master-spectrum.mzML")

print("image spectrum saved to graphic and .mzML.")
