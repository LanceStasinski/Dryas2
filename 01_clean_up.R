################################################################################
#Clean up scans
#Data can be obtained from shared Google Drive folder:
#FM_Polar_Studies > Dryas_Spectral_Analyses > Scans_raw

################################################################################
#Load Packages
################################################################################

library("spectrolab")

################################################################################
#Load Functions
################################################################################

add_meta <- function(spectra_path, metadata_path){
  spectra_raw = read_spectra(path = spectra_path, format = "sed")
  metadata = read.csv(file = metadata_path, header = TRUE, stringsAsFactors = FALSE)
  meta(spectra_raw) <- metadata
  return(spectra_raw)
}


#This function adds the metadata, cuts spectra to wavelength 400:2400, removes
#reflectance values greater than 1, reduces the data to the mean reflectance
# of each individual plant, and smooths the spectra

thebigclean <- function(spectra_path, metadata_path){
  meta.spectra = add_meta(spectra_path, metadata_path)
  spectra_cut = meta.spectra[, 400:2400]
  spec1 <- spectra_cut[!rowSums(spectra_cut > 1),]
  metadata = meta(spec1)
  ag.spec = spectrolab:::aggregate.spectra(spec1,
                                                 metadata$Name, mean,
                                                 try_keep_txt(mean))
  clean_spectra = smooth(ag.spec)
  return(clean_spectra)
}


################################################################################
#Set working directory to folder containing downloaded spectral data
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#Twelve Mile
################################################################################

tm_path = "Scans_raw/Twelve_Mile/12mile-dryas-08082019"
tm_meta = "Scans_raw/Twelve_Mile/twelvemile.csv"
tm_clean = thebigclean(tm_path, tm_meta)
tm_vn = normalize(tm_clean)

################################################################################
#Bison Gulch
################################################################################

bg_path = "Scans_raw/Bison_Gulch/bisonGulch_dry"
bg_meta = "Scans_raw/Bison_Gulch/bisongulch2.csv"
bg_clean = thebigclean(bg_path, bg_meta)
bg_vn = normalize(bg_clean)

################################################################################
#Eagle Summit
################################################################################    

es_path = "Scans_raw/Eagle_Summit/es-dry"
es_meta = "Scans_raw/Eagle_Summit/eaglesummit1.csv"
es_clean = thebigclean(es_path, es_meta)
es_vn = normalize(es_clean)

################################################################################
#Murphy Dome B
################################################################################

mdb_path = "Scans_raw/Murphy_Dome_B/murphyB-dry"
mdb_meta = "Scans_raw/Murphy_Dome_B/murphydomeb1.csv"
mdb_clean = thebigclean(mdb_path, mdb_meta)
mdb_vn = normalize(mdb_clean)

################################################################################
#Wickersham Dome A
################################################################################    

wda_path = "Scans_raw/Wickersham_A/wickershamdome-dry"
wda_meta = "Scans_raw//Wickersham_A/wickershama1.csv"
wda_clean = thebigclean(wda_path, wda_meta)
wda_vn = normalize(wda_clean)

################################################################################
#Wickersham Dome B
################################################################################    

wdb_path = "Scans_raw/Wickersham_B/wdb-dry"
wdb_meta = "Scans_raw/Wickersham_B/wdb_revised.csv"
wdb_clean = thebigclean(wdb_path, wdb_meta)
wdb_vn = normalize(wdb_clean)

################################################################################
#combine
################################################################################

clean1 = combine(tm_clean, es_clean)
clean_big3 = combine(clean1, wdb_clean)
clean2 = combine(clean_big3, mdb_clean)
clean3 = combine(clean2, wda_clean)
clean_all = combine(clean3, bg_clean)
all_vn = normalize(clean_all)

vn_big3 = normalize(clean_big3)

################################################################################
#separate by three primary sites
################################################################################

meta3 = meta(clean_big3)
ala3 = clean_big3[meta3$Species == "alaskensis"]
oct3 = clean_big3[meta3$Species == "octopetala"]

big3.no_hybrids = combine(ala3, oct3)

vn_big3.no_hybrids = normalize(big3.no_hybrids)

################################################################################
#Save Data
################################################################################

#Save Cleaned Spectra
saveRDS(tm_clean, "Clean-up/Clean_spectra/tm_clean.rds")
saveRDS(bg_clean, "Clean-up/Clean_spectra/bg_clean.rds")
saveRDS(es_clean, "Clean-up/Clean_spectra/es_clean.rds")
saveRDS(mdb_clean, "Clean-up/Clean_spectra/mdb_clean.rds")
saveRDS(wda_clean, "Clean-up/Clean_spectra/wda_clean.rds")
saveRDS(wdb_clean, "Clean-up/Clean_spectra/wdb_clean.rds")
saveRDS(clean_all, "Clean-up/Clean_spectra/clean_all.rds")
saveRDS(clean_big3, "Clean-up/Clean_spectra/clean_big3.rds")
saveRDS(big3.no_hybrids, "Clean-up/Clean_spectra/big3.no_hybrids.rds")

#Save Vector Normalized Spectra
saveRDS(tm_vn, "Clean-up/Vector_normalized/tm_vn.rds")
saveRDS(bg_vn, "Clean-up/Vector_normalized/bg_vn.rds")
saveRDS(es_vn, "Clean-up/Vector_normalized/es_vn.rds")
saveRDS(mdb_vn, "Clean-up/Vector_normalized/mdb_vn.rds")
saveRDS(wda_vn, "Clean-up/Vector_normalized/wda_vn.rds")
saveRDS(wdb_vn, "Clean-up/Vector_normalized/wdb_vn.rds")
saveRDS(all_vn, "Clean-up/Vector_normalized/all_vn.rds")
saveRDS(vn_big3, "Clean-up/Vector_normalized/vn_big3.rds")
saveRDS(vn_big3.no_hybrids, "Clean-up/Vector_normalized/vn_big3.no_hybrids.rds")


