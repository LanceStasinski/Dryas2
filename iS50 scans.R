#iS50 scans


library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50/ES/iS50_ESA")

################################################################################
#format csv function
################################################################################
#format csvs output by iS50 into those that can be read as spectra
format.spec.csv = function(x){
  spec1 = t(read.csv(x, head = F))
  colnames(spec1) <- spec1[1,]
  spec2 = t(as.data.frame(spec1[-1,]))
  rownames(spec2) <- x
  return(spec2)
}

#convert csvs to spectra then smooth spectra
csv.to.spec = function(files){
  csv.stack = do.call(rbind, lapply(files, format.spec.csv))
  spectra = smooth(as.spectra(csv.stack))
}

################################################################################
#create spectra from csvs
################################################################################
#ES
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50/ES/iS50_ESA")
es = list.files()
es_iS50 = csv.to.spec(es)

#TM
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50/TM/iS50_TM")
tm = list.files()
tm_iS50 = csv.to.spec(tm)

#WDB
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50/WDB/iS50_WDB")
wdb = list.files()
wdb_iS50 = csv.to.spec(wdb)

################################################################################
#Add Meta Data
################################################################################
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Scans_raw/iS50")

es_meta = read.csv(file = "ES/es_iS50.csv", header = TRUE, stringsAsFactors = FALSE)
meta(es_iS50) <- es_meta

tm_meta = read.csv("TM/tm_iS50.csv", header = TRUE, stringsAsFactors = FALSE)
meta(tm_iS50) <- tm_meta

wdb_meta = read.csv("WDB/wdb_iS50.csv", header = TRUE, stringsAsFactors = FALSE)
meta(wdb_iS50) <- wdb_meta

################################################################################
#combine
################################################################################

spec_iS50 = Reduce(spectrolab::combine, list(es_iS50, tm_iS50, wdb_iS50))

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Clean-up/Clean_spectra")
saveRDS(spec_iS50, "spec_iS50.rds")


ispec = readRDS('spec_iS50.rds')
data = meta(ispec)
ispec.m = as.matrix(ispec)
colnames(ispec.m) <- 10000000/as.integer(colnames(ispec.m))
ispec.m2 = ispec.m[, order(ncol(ispec.m):1)]

iS50_spec = as.spectra(ispec.m2)
meta(iS50_spec) <- data
plot(iS50_spec)
saveRDS(iS50_spec, "spec_iS50.rds")
