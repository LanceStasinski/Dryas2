#Random
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(dplyr)
library(randomForest)
library(cluster)

setwd("~/GitHub/Dryas2")
################################################################################
#Data setup 
################################################################################

#data
spec_all = readRDS("Data/clean_all.rds")

#add new population delineations
s.m = as_spectra(as.matrix(spec_all))
meta(s.m) = read.csv('Data/metadata_2.csv', stringsAsFactors = F)
spec_all = s.m

#remove any NaN values
##Note: this removes hybrids
spec_all = spec_all[!meta(spec_all)$GenePop_ID == "NaN",]

spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm to reduce redundant data
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
spec_df = as.data.frame(spec_mat)

#add important classes
spec_df$Species_ID = spec_all.df$Species_ID
spec_df$Location = spec_all.df$Location
spec_df$GenePop_ID  = spec_all.df$GenePop_ID
spec_df$sp_loc = spec_all.df$sp_loc

################################################################################
#Unsupervised Random Forest 
################################################################################

s.kappa <- c()
p.kappa <- c()
l.kappa <- c()
sl.kappa <- c()


for(i in 1:100){
#model
rf = randomForest(x = spec_df[,-c(202,203,204,205)], ntree= 2000, proximity = T)


#predict classes from proximity matrix extracted from model
prox = rf$proximity

#species
pam2 = pam(prox, 2)
pred.sp = cbind(pam2$clustering, spec_df$Species_ID)
sp.t = table(pred.sp[,2],pred.sp[,1])
colnames(sp.t) = rownames(sp.t)
s.cm = confusionMatrix(sp.t)
s.kap = assign(paste0("s.kap",i), s.cm$overall[2])
s.kappa <- append(s.kappa, get('s.kap'))

#population
pam6 = pam(prox, 6)
pred.pop = cbind(pam6$clustering, spec_df$GenePop_ID)
pop.t = table(pred.pop[,2],pred.pop[,1])
colnames(pop.t) = rownames(pop.t)
p.cm = confusionMatrix(pop.t)
p.kap = assign(paste0("p.kap",i), p.cm$overall[2])
p.kappa <- append(p.kappa, get('p.kap'))

#location
pam6 = pam(prox, 6)
pred.loc = cbind(pam6$clustering, spec_df$Location)
loc.t = table(pred.loc[,2],pred.loc[,1])
colnames(loc.t) = rownames(loc.t)
l.cm = confusionMatrix(loc.t)
l.kap = assign(paste0("l.kap",i), l.cm$overall[2])
l.kappa <- append(l.kappa, get('l.kap'))

#species+location
pam9 = pam(prox, 9)
pred.sploc = cbind(pam9$clustering, spec_df$sp_loc)
sploc.t = table(pred.sploc[,2],pred.sploc[,1])
colnames(sploc.t) = rownames(sploc.t)
sl.cm = confusionMatrix(sploc.t)
sl.kap = assign(paste0("sl.kap",i), sl.cm$overall[2])
sl.kappa <- append(sl.kappa, get('sl.kap'))

}

#Results
df = data.frame(row.names = c("Species", "Population", "Location", "Species+Loc"))

s.mean = mean(s.kappa)
p.mean = mean(p.kappa)
l.mean = mean(l.kappa)
sl.mean = mean(sl.kappa)
k.mean = c(s.mean, p.mean, l.mean, sl.mean)
df = cbind(df, k.mean)

s.sd = sd(s.kappa)
p.sd = sd(p.kappa)
l.sd = sd(l.kappa)
sl.sd = sd(sl.kappa)
k.sd = c(s.sd, p.sd, l.sd, sl.sd)
df = cbind(df, k.sd)

colnames(df) = c("Kappa Mean", "Kappa SD")
write.csv(df, file="rf_confusion_matrix_output/kappa_stats.csv")











