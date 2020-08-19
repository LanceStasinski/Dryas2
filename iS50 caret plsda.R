#is50 PLSDA
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
ispec = readRDS('Clean-up/Clean_spectra/spec_iS50.rds')
spec_all = ispec
spec_all = spec_all[!meta(spec_all)$sp_loc == "NaN"]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX"]

spec_all.m = as.matrix(spec_all)
spec_mat = spec_all.m
spec_all.df = as.data.frame(spec_all)


#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$sp_loc)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$sp_loc"] <- "sp_loc"


for(i in 1:10){
  
  set.seed(i)
  
  inTrain <- caret::createDataPartition(
    y = spec_df$sp_loc,
    p = .8,
    list = FALSE
  )
  
  training <- spec_df[inTrain,]
  testing <- spec_df[-inTrain,]
  
  #tune model
  ctrl <- trainControl(
    method = "repeatedcv", 
    number = 10,
    repeats = 3)
  
  
  plsFit <- train(
    sp_loc ~ .,
    data = training,
    method = "pls",
    preProc = c("center", "scale"),
    trControl = ctrl,
    tuneLength = 100)
  
  assign(paste0('plsFit', i), plsFit)
  
  loadings = plsFit$finalModel$loadings
  loadings.m = as.matrix(loadings)
  class(loadings.m) <- 'matrix'
  assign(paste0('lm',i), loadings.m)
  
  comp1 = loadings.m[,1]
  assign(paste0('comp1_',i), comp1)
  comp2 = loadings.m[,2]
  assign(paste0('comp2_',i), comp2)
  comp3 = loadings.m[,3]
  assign(paste0('comp3_',i), comp3)
  
  
  
  #test model
  plsClasses <- predict(plsFit, newdata = testing)
  
  
  #Confusion matrices
  cm = confusionMatrix(data = plsClasses, testing$sp_loc)
  
  
  
  cm.m = as.matrix(cm)
  
  assign(paste0("cm", i), cm.m)
}

cm.total = (cm1 + cm2 + cm3 + cm4 + cm5 + cm6 + cm7 + cm8+ cm9 + cm10)/10
cm.total = as.matrix(cm)
cm.total = t(cm.total)
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('DA_es', 'DA_tm', 'DA_wdb', 'DO_es', 'DO_tm', 'DO_wdb', 
                        'DX_wdb' )

par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
corrplot(cm.total, is.corr = F, method = 'color', addCoef.col = 'darkorange2',
         tl.srt = 90, tl.offset = 1.5, number.digits = 2, tl.cex = .75,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Reference", side = 2, line = -5, cex = 1.5)
mtext("Prediction", side = 3, cex = 1.5, at = 2, line = 2)
