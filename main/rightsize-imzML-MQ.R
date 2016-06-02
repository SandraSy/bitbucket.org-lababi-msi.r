# Clear Workspace
rm(list=ls())

# Load rmsi library
library("rmsi")

# Load the imzML data
filename <- "/home/rob/MS-DATA/MSI-urinary-bladder/HR2MSI mouse urinary bladder S096.imzML"

imagespectra <- importImzMl(filename, centroided=TRUE)

elementsimzML<-length(imagespectra)

imagexy <- sapply(imagespectra, function(x)metaData(x)$imaging$pos)
imagexy <- t(imagexy)
colnames(imagexy) <- c("x","y")

specno <-1

imagemz <- list(1:elementsimzML)

while (specno<=elementsimzML)
{
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  singlespec <- createMassSpectrum(mass=mzs, intensity=counts)
  imagemz[[specno]] <- singlespec
  singlespec <- detectPeaks(singlespec, SNR=10)
  print(specno)
  specno=specno+1
}

imagemz <- topN(imagemz, n=30)

exportImzMl(imagemz, path="bladder.imzML", force=TRUE, coordinates=imagexy)

