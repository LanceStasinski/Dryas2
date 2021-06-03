################################################################################
#Clean up scans
#Data can be obtained from shared Google Drive folder:
#FM_Polar_Studies > Dryas_Spectral_Analyses > Scans_raw

################################################################################
#Load Packages
################################################################################

library("spectrolab")

################################################################################
#Prerequisite functions
################################################################################
#Read in spectra and metadata. Remove scans names that do not match row names in 
#the metadata. Add metadata to the spectra
add_meta <- function(spectra_path, metadata_path){
  spectra_raw = read_spectra(path = spectra_path, format = "sed")
  metadata = read.csv(file = metadata_path, header = TRUE, 
                      stringsAsFactors = FALSE)
  spec_m = as.matrix(spectra_raw)
  rownames(spec_m) = sub(".sed", "", rownames(spec_m))
  remove = setdiff(rownames(spec_m), metadata$scan.name)
  spec_m = spec_m[!rownames(spec_m) %in% remove, ]
  spectra_raw = as_spectra(spec_m)
  meta(spectra_raw) = metadata
  return(spectra_raw)
}

#take the mean reflectance for each individual plant
trim.spectra = function(spectra){
  avg_spec = aggregate(spectra, by = meta(spectra)$Name, mean, try_keep_txt(mean))
  return(avg_spec)
}

################################################################################
#Primary functions
################################################################################
#This function adds the metadata, cuts spectra to wavelength 400:2400, removes
#reflectance values greater than 1, reduces the data to the mean for each 
#individual and smooths the spectra.
thebigclean <- function(spectra_path, metadata_path){
  meta.spectra = add_meta(spectra_path, metadata_path)
  spectra_cut = meta.spectra[, 400:2400]
  spec1 = spectra_cut[!rowSums(spectra_cut > 1),]
  spec2 = trim.spectra(spec1)
  clean_spectra = smooth(spec2)
  return(clean_spectra)
}


################################################################################
#Set working directory to folder containing downloaded spectral data
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
#Twelve Mile
################################################################################

tm_path = "Scans_raw/Dry Scans 2/Twelve_Mile/12mile-dryas-08082019"
tm_meta = "Scans_raw/Dry Scans 2/Twelve_Mile/tm_pops.csv"
tm_clean1 = thebigclean(tm_path, tm_meta)
tm_clean = tm_clean1[!meta(tm_clean1)$Clone == "Yes",]


################################################################################
#Bison Gulch
################################################################################

bg_path = "Scans_raw/Dry Scans 2/Bison_Gulch/bisonGulch_dry"
bg_meta = "Scans_raw/Dry Scans 2/Bison_Gulch/bg_pops.csv"
bg_clean = thebigclean(bg_path, bg_meta)

################################################################################
#Eagle Summit
################################################################################    

es_path = "Scans_raw/Dry Scans 2/Eagle_Summit/es-dry"
es_meta = "Scans_raw/Dry Scans 2/Eagle_Summit/es_pops.csv"
es_clean1 = thebigclean(es_path, es_meta)

es_clean2 = es_clean1[!meta(es_clean1)$Clone == "Yes",]
es_clean = es_clean2[!meta(es_clean2)$Notes == "not sequenced",]


################################################################################
#Murphy Dome B
################################################################################

mdb_path = "Scans_raw/Dry Scans 2/Murphy_Dome_B/murphyB-dry"
mdb_meta = "Scans_raw/Dry Scans 2/Murphy_Dome_B/mdb_pops.csv"
mdb_clean = thebigclean(mdb_path, mdb_meta)


################################################################################
#Wickersham Dome A
################################################################################    

wda_path = "Scans_raw/Dry Scans 2/Wickersham_A/wickershamdome-dry"
wda_meta = "Scans_raw/Dry Scans 2/Wickersham_A/wda_pops.csv"
wda_clean = thebigclean(wda_path, wda_meta)


################################################################################
#Wickersham Dome B
################################################################################    

wdb_path = "Scans_raw/Dry Scans 2/Wickersham_B/wdb-dry"
wdb_meta = "Scans_raw/Dry Scans 2/Wickersham_B/wdb_pops.csv"
wdb_clean = thebigclean(wdb_path, wdb_meta)



################################################################################
#combine
################################################################################

#clean_all is the object used for most analyses
clean_all = Reduce(combine, list(tm_clean, es_clean, wdb_clean, mdb_clean, 
                                 wda_clean, bg_clean))
saveRDS(clean_all, "Clean-up/Clean_spectra/clean_all.rds")


################################################################################
################################################################################
#Wet Scans
################################################################################
#Eagle Summit
es_w_path = "Scans_raw/Wet Scans/es_wet"
es_w_meta = "Scans_raw/Wet Scans/es_wet/es_pops_wet.csv"
es_w_clean1 = thebigclean(es_w_path, es_w_meta)
es_w_clean2 = es_w_clean1[!meta(es_w_clean1)$Clone == "Yes",]
es_w_clean = es_w_clean2[!meta(es_w_clean2)$Notes == "not sequenced",]
vn_es_w = normalize(es_w_clean)

#Wickersham Dome A
wda_w_path = "Scans_raw/Wet Scans/wda_wet"
wda_w_meta = "Scans_raw/Wet Scans/wda_wet/wda_pops_wet.csv"
wda_w_clean = thebigclean(wda_w_path, wda_w_meta)
vn_wda_w = normalize(wda_w_clean)

#Wickerhamd Dome B
wdb_w_path = "Scans_raw/Wet Scans/wdb_wet"
wdb_w_meta = "Scans_raw/Wet Scans/wdb_wet/wdb_pops_wet.csv"
wdb_w_clean = thebigclean(wdb_w_path, wdb_w_meta)
vn_wdb_w = normalize(wdb_w_clean)

#all wet
Clean_all_w = Reduce(combine, list(wdb_w_clean, wda_w_clean, es_w_clean))
vn_all_w = normalize(Clean_all_w)
saveRDS(Clean_all_w, "Clean-up/Clean_spectra/clean_all_w.rds")
saveRDS(vn_all_w, "Clean-up/Vector_normalized/vn_all_w.rds")

