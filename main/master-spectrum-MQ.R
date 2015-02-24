# This script summarizes the x (e.g. 100) most intense peaks/ point of each spectrum in one single spectrum
# The increasing number indicates the number of processed spectra
# The final spectrum is exported into a .mzML file
# A density distribution map is printed to a publication quality .tiff

# You can change the number of selected peaks with (this would be the setting for 100 signals):
#if (length(newspec)>"200"){        //2x the number of peaks 
#  newspec <- newspec[1:100,]}      //1:x strongest signals

# Clear Workspace
rm(list=ls())
# set working directory to source file
setwd("~/SRC/MSI.R/main")

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")

# Load the imzML data
#filename <- "/home/rob/SRC/MSI.R/SampleData/test.imzML"
filename <- "/home/rob/MS-DATA/LTP-MSI-chilli/ltpmsi-chilli.imzML"

imagespectra <- importImzMl(filename, centroided=TRUE)
#plotImsSlice(imagespectra)
elementsimzML<-length(imagespectra)

specno <-1

while (specno<=elementsimzML)
{
    mzs <- mass(imagespectra[[specno]])
    counts <- intensity(imagespectra[[specno]])
    newspec  <- cbind(mzs, counts)
    newspec <- newspec[order(-counts),]
    if (length(newspec)>"200"){
    newspec <- newspec[1:100,]}
    
    mzs<-newspec[,1]
    counts<-newspec[,2]
    
    if(exists("mzsum")=="TRUE"){mzsum<-c(mzsum,mzs)}
    if(exists("countsum")=="TRUE"){countsum<-c(countsum,counts)}
    if(exists("mzsum")=="FALSE"){mzsum<-mzs}
    if(exists("countsum")=="FALSE"){countsum<-counts}
    specmatrix <- cbind(mzsum, countsum)  
    print(specno)
  specno=specno+1
}

mzsumdensity<-density(mzsum)

# printing spectrum. format options could be e.g. png, eps, tiff, pdf. You can control the resolution with res
tiff(filename="spectrumdensity.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
plot(mzsumdensity,,xlab="m/z",ylab="density",main="m/z density distribution")
polygon(mzsumdensity, col="grey") 
dev.off()

## export a single spectrum
s <- list(createMassSpectrum(mass=specmatrix[,1], intensity=specmatrix[,2]))
exportMzMl(s[[1]], force=TRUE, file="masterspectrum.mzML")

print("image spectrum saved to graphic and .mzML.")
