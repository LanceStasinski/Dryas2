################################################################################
#Set up
################################################################################

library(tidyverse)
library(caret)
library(pls)
library(spectrolab)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spectra = readRDS("Clean-up/Clean_spectra/clean_all.rds")
spectra = spectra[!meta(spectra)$DA == "NaN",]

es = spectra[meta(spectra)$Location == "Eagle Summit",]
tm = spectra[meta(spectra)$Location == "Twelve Mile",]
wdb = spectra[meta(spectra)$Location == "Wickersham Dome B",]
spectra = Reduce(spectrolab::combine, list(es, tm, wdb))


spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)
spec_df = as.data.frame(spectra.m)
spec_df = cbind(spec_df, spectra.df$DA)
spec_df = cbind(spec_df, spectra.df$Location)
colnames(spec_df)[colnames(spec_df) == "spectra.df$DA"] <- "DA"
colnames(spec_df)[colnames(spec_df) == "spectra.df$Location"] <- "Location"
spec_df$Location <- as.factor(spec_df$Location)

spec_df2 = cbind(spec_df, spectra.df$Species_ID)
colnames(spec_df2)[colnames(spec_df2) == "spectra.df$Species_ID"] <- "Species_ID"

################################################################################
#Training and testing sets - all
################################################################################

set.seed(12)
training.samples <- 
  createDataPartition(spec_df$DA, p = 0.8, list = FALSE)

train.data  <- spec_df[training.samples, ]
test.data <- spec_df[-training.samples, ]

train.plot <-spec_df2[training.samples, ]
test.plot <- spec_df2[-training.samples, ]
test.plot$color = "black"
test.plot$color[test.plot$Species_ID == "DO"]="red"
test.plot$color[test.plot$Species_ID == "DA"]="blue"

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

write.csv(results, "Figures/PLSR/DA_plsr20_3sites_location.csv")

################################################################################
#Plot
################################################################################

par(mfrow = c(1,1))
plot(test.plot$DA, predictions, xlab = "Actual", ylab = "Predicted", xlim = c(-.3,1.3), ylim = c(-.3, 1.3),
     main = expression("Predicting Proportion of "*italic("Dryas alaskensis")*" ancestry"))
lines(x = c(-2, 2), y = c(-2, 2))
points(test.plot$DA, predictions, col = test.plot$color, pch = 16)
legend("bottomright", inset = 0.01,
       legend=c("D. octopetala", "Hybrid", "D. alaskensis"), 
       text.font = c(3,1,3),
       col=c("red", "black", "blue"), pch = 16)



