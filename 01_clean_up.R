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
#add metadata to raw spectra
add_meta <- function(spectra_path, metadata_path){
  spectra_raw = read_spectra(path = spectra_path, format = "sed")
  metadata = read.csv(file = metadata_path, header = TRUE, 
                      stringsAsFactors = FALSE)
  meta(spectra_raw) <- metadata
  return(spectra_raw)
}



################################################################################
#Primary functions
################################################################################
#This function adds the metadata, cuts spectra to wavelength 400:2400, removes
#reflectance values greater than 1, reduces the data to the 4 measurements
#that are closest to the mean for each individual and smooths the spectra.
thebigclean <- function(spectra_path, metadata_path){
  meta.spectra = add_meta(spectra_path, metadata_path)
  spectra_cut = meta.spectra[, 400:2400]
  spec1 = spectra_cut[!rowSums(spectra_cut > 1),]
  clean_spectra = smooth(spec1)
  return(clean_spectra)
}


################################################################################
#Set working directory to folder containing downloaded spectral data
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
#Twelve Mile
################################################################################

tm_path = "Scans_raw/Dry Scans/Twelve_Mile/12mile-dryas-08082019"
tm_meta = "Scans_raw/Dry Scans/Twelve_Mile/tm_pops.csv"
tm_clean1 = thebigclean(tm_path, tm_meta)
tm_clean = tm_clean1[!meta(tm_clean1)$Clone == "Yes",]


################################################################################
#Bison Gulch
################################################################################

bg_path = "Scans_raw/Dry Scans/Bison_Gulch/bisonGulch_dry"
bg_meta = "Scans_raw/Dry Scans/Bison_Gulch/bg_pops.csv"
bg_clean = thebigclean(bg_path, bg_meta)


################################################################################
#Eagle Summit
################################################################################    

es_path = "Scans_raw/Dry Scans/Eagle_Summit/es-dry"
es_meta = "Scans_raw/Dry Scans/Eagle_Summit/es_pops.csv"
es_clean1 = thebigclean(es_path, es_meta)

es_clean2 = es_clean1[!meta(es_clean1)$Clone == "Yes",]
es_clean = es_clean2[!meta(es_clean2)$Notes == "not sequenced",]


################################################################################
#Murphy Dome B
################################################################################

mdb_path = "Scans_raw/Dry Scans/Murphy_Dome_B/murphyB-dry"
mdb_meta = "Scans_raw/Dry Scans/Murphy_Dome_B/mdb_pops.csv"
mdb_clean = thebigclean(mdb_path, mdb_meta)


################################################################################
#Wickersham Dome A
################################################################################    

wda_path = "Scans_raw/Dry Scans/Wickersham_A/wickershamdome-dry"
wda_meta = "Scans_raw/Dry Scans/Wickersham_A/wda_pops.csv"
wda_clean = thebigclean(wda_path, wda_meta)


################################################################################
#Wickersham Dome B
################################################################################    

wdb_path = "Scans_raw/Dry Scans/Wickersham_B/wdb-dry"
wdb_meta = "Scans_raw/Dry Scans/Wickersham_B/wdb_pops.csv"
wdb_clean = thebigclean(wdb_path, wdb_meta)



################################################################################
#combine
################################################################################
#clean_all is the object used for most analyses
clean_all = Reduce(combine, list(tm_clean, es_clean, wdb_clean, mdb_clean, 
                                 wda_clean, bg_clean))
#write.csv(meta(clean_all), 'meta_6.csv')
#correct populations in excel
meta = read.csv('meta_6.csv', stringsAsFactors = F)
meta(clean_all) = meta
saveRDS(clean_all, "Clean-up/Clean_spectra/clean_all_6scans.rds")

#code for averaging spectra by layers
l1 = clean_all[meta(clean_all)$layer == 1,]
l1.avg = aggregate(l1, by = meta(l1)$Name, mean, try_keep_txt(mean))

l2 = clean_all[meta(clean_all)$layer == 2,]
l2.avg = aggregate(l2, by = meta(l2)$Name, mean, try_keep_txt(mean))

l3 = clean_all[meta(clean_all)$layer == 3,]
l3.avg = aggregate(l3, by = meta(l3)$Name, mean, try_keep_txt(mean))

clean_all2 = Reduce(combine, list(l1.avg, l2.avg, l3.avg))
write.csv(meta(clean_all2), "meta_3lyr.csv")
#fix populations in excel
meta_new = read.csv("meta_3lyr.csv", stringsAsFactors = F)
meta(clean_all2) = meta_new

saveRDS(clean_all2, "Clean-up/Clean_spectra/clean_all_3lyr.rds")

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

