################################################################################
#Set up
################################################################################

library(tidyverse)
library(caret)
library(pls)
library(spectrolab)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

spectra = readRDS("Clean-up/Vector_normalized/all_vn.rds")
es = spectra[meta(spectra)$Location == "Eagle Summit",]
tm = spectra[meta(spectra)$Location == "Twelve Mile",]
wdb = spectra[meta(spectra)$Location == "Wickersham Dome B",]
spectra = Reduce(spectrolab::combine, list(es, tm, wdb))
spectra = spectra[!meta(spectra)$DA == "NaN",]

spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)
spectra.m = cbind(DA = spectra.df$DA, spectra.m)
spectra.m = cbind(Location = as.factor(spectra.df$Location), spectra.m)

################################################################################
#Training and testing sets - all
################################################################################

set.seed(12)
training.samples <- 
  createDataPartition(spectra.m[,"DA"], p = 0.6, list = FALSE)

train.data  <- spectra.m[training.samples, ]
test.data <- spectra.m[-training.samples, ]

################################################################################
#PLS R - all
################################################################################

set.seed(12)
model <- train(
  DA~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 30),
  tuneLength = 50
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune

summary(model$finalModel)

# Make predictions
predictions <- model %>% predict(test.data)
# Model performance metrics
results = data.frame(
  RMSE = caret::RMSE(predictions, test.data[,"DA"]),
  Rsquare = caret::R2(predictions, test.data[,"DA"])
)

write.csv(results, "Figures/DA_plsr.csv")
