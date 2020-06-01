################################################################################
#Set up
################################################################################

library(tidyverse)
library(caret)
library(pls)
library(spectrolab)


setwd("~/Dryas Research/Data")

spectra = readRDS("sm_all.rds")
all.df = as.data.frame(spectra)
all.m = as.matrix(spectra)
all.m1 = cbind( LMA = all.df$LMA, all.m)
all.m2 = cbind(Location = as.factor(all.df$Location), all.m1)

keep2 = !is.na(all.m2[,"LMA"])
all.m4 = all.m2[keep2,]
all.m5 = all.m4["Location" = c(2, 5, 7),]

#Remove NA values
keep = !is.na(all.m1[,"LMA"])
all.m3 = all.m1[keep,]

################################################################################
#Training and testing sets - all
################################################################################

set.seed(12)
training.samples <- 
  createDataPartition(all.m3[,"LMA"], p = 0.8, list = FALSE)

train.data  <- all.m3[training.samples, ]
test.data <- all.m3[-training.samples, ]

################################################################################
#PLS R - all
################################################################################

set.seed(12)
model <- train(
  LMA~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 20),
  tuneLength = 20
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
data.frame(
  RMSE = caret::RMSE(predictions, test.data[,"LMA"]),
  Rsquare = caret::R2(predictions, test.data[,"LMA"])
)


################################################################################
#Training and testing sets - 3 sites
################################################################################

set.seed(78)
training.samples <- 
  createDataPartition(all.m5[,"LMA"], p = 0.8, list = FALSE)

train.data  <- all.m5[training.samples, ]
test.data <- all.m5[-training.samples, ]

################################################################################
#PLS R - all
################################################################################

set.seed(78)
model <- train(
  LMA~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 40),
  tuneLength = 40
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
data.frame(
  RMSE = caret::RMSE(predictions, test.data[,"LMA"]),
  Rsquare = caret::R2(predictions, test.data[,"LMA"])
)


################################################################################
#PCR
################################################################################

set.seed(78)
model <- train(
  LMA~., data = train.data, method = "pcr",
  scale = TRUE,
  trControl = trainControl("cv", number = 50),
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
data.frame(
  RMSE = caret::RMSE(predictions, test.data[,"LMA"]),
  Rsquare = caret::R2(predictions, test.data[,"LMA"])
)






