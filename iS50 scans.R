#iS50 scans


library(spectrolab)
library(readr)
library(dplyr)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50/ES/iS50_ESA")

################################################################################
#read spectra
################################################################################




spec1 = read.csv("ESA20-1.csv", header = F)
spec2 = t(spec1)
colnames(spec2) <- spec2[1,]
spec3 = spec2[-1,]

colnames(spec1) <- c("wavelength", "reflectance")

read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}

spec4 = read.tcsv("ESA20-1.csv")
