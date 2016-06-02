# This script reads an .imzML file and prints the number of spectra

# Clear Workspace
rm(list=ls())

# Load rmsi library
library("rmsi")

# Load the imzML data
filename <- "mq.imzML"

imagespectra <- importImzMl(filename, centroided=TRUE)

elementsimzML<-length(imagespectra)

print("Number of spectra in tested .imzML:")
print(elementsimzML)