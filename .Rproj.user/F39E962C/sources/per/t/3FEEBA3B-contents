#Caret PLSDA
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(mlbench)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
# Fit PLS_DA model all dry
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Species_ID)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Species_ID"] <- "Species_ID"


#Partition Data

set.seed(1)

inTrain <- caret::createDataPartition(
  y = spec_df$Species_ID,
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
  Species_ID ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  trControl = ctrl,
  tuneLength = 200)

plsFit

#test model
plsClasses <- predict(plsFit, newdata = testing)

confusionMatrix(data = plsClasses, testing$Species_ID)







cm.m = as.matrix(cm)





assign(paste0("cm", i), cm.m)





cm.total = (cm1 + cm2 + cm3 + cm4 + cm5 + cm6 + cm7 + cm8+ cm9 + cm10)/10
