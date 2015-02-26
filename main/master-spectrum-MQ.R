# This script summarizes the x (e.g. 100) most intense peaks/ point of each spectrum in one single spectrum
# The increasing number indicates the number of processed spectra
# The final spectrum is exported into a .mzML file and printed to a .tiff
# A density distribution map is printed to a publication quality .tiff
# Further, peaks are picked and exported to a .csv text file and added to a zoom spectrum
# Please adjust the paramenters for peak picking (e.g. signal/noise) and spectrum zoom below

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

mzsumdensity<-density(mzsum, bw=0.01)

# printing spectrum. format options could be e.g. png, eps, tiff, pdf. You can control the resolution with res
tiff(filename="spectrumdensity.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
plot(mzsumdensity,,xlab="m/z",ylab="density",main="m/z density distribution")
polygon(mzsumdensity, col="grey") 
dev.off()

# aggregate spectrum data
specmatrix<-aggregate(countsum ~ mzsum, data=specmatrix, FUN=sum)
specmatrix<-specmatrix[ order(specmatrix[,1]), ]
s <- list(createMassSpectrum(mass=specmatrix[,1], intensity=specmatrix[,2]))

# print master spectrum
tiff(filename="masterspectrum.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
plot(s[[1]], xlab="m/z")
dev.off()

# export a single mzML spectrum
exportMzMl(s[[1]], force=TRUE, file="masterspectrum.mzML")

# peak picking and export to a peaklist.csv file
peaks <-  detectPeaks(s, SNR = 3)

tiff(filename="spectrum-with-peaks.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
plot(s[[1]],xlim=c(40, 550), xlab="m/z")
points(peaks[[1]],pch=18)
dev.off()

# export of peaklist
peaklist<-intensityMatrix(peaks)
tpeaklist<-t(peaklist)
colnames(tpeaklist)<-c("intensity")
write.csv(tpeaklist, file="peaklist.csv")

print("Files generated: Pseudo master spectrum (tiff and mzML), spectrum density (tiff), spectrum-zoom with peaks (tiff, peaklist.csv)")
