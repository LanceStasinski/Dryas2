#PLS - predicting ancestry from spectra
################################################################################
#Set up
################################################################################

library(tidyverse)
library(caret)
library(pls)
library(spectrolab)
library(matrixStats)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#Data setup 
################################################################################
spectra= readRDS("Clean-up/Clean_spectra/clean_all.rds")

spectra = spectra[!meta(spectra)$DA == "NaN",]

spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)
spec_df = as.data.frame(spectra.m)
spec_df = cbind(spec_df, spectra.df$DA)
colnames(spec_df)[colnames(spec_df) == 'spectra.df$DA'] <- 'DA'

#These objects are used for plotting later
spec_df2 = cbind(spec_df, spectra.df$Species_ID)
colnames(spec_df2)[colnames(spec_df2) == "spectra.df$Species_ID"]<- "Species_ID"
spec_df2 = cbind(spec_df2, spectra.df$Name)
colnames(spec_df2)[colnames(spec_df2) == "spectra.df$Name"] <- 'Name'

################################################################################
#Training and testing sets - all
################################################################################
#set number of components
ncomp = 23 

#create vectors and matrix to store metrics
rmse.fit <- matrix(nrow = ncomp)
error.list <- c()
r2 <- c()

#start of PLS code
for( i in 1:100){

#partition the data: 70% for training, 30% for testing
training.samples <- 
  createDataPartition(spec_df$DA, p = 0.7, list = FALSE)

train.data  <- spec_df[training.samples, ]
test.data <- spec_df[-training.samples, ]


#Fit model
plsFit <- train(
  DA~., data = train.data, method = "pls",
  trControl = trainControl("repeatedcv", number = 10, repeats = 3),
  tuneLength = ncomp
)

#objects for determining optimal number of components
rmse = assign(paste0('rmse', i), as.matrix(plsFit$results$RMSE))
rmse.fit <- cbind(rmse.fit, get('rmse'))

#Make predictions
predictions <- plsFit %>% predict(test.data)

#Model performance metrics
error = assign(paste0('error', i), caret::RMSE(predictions, test.data[,"DA"]))
error.list <- append(error.list, get('error'))

Rsquare = assign(paste0('Rsquare', i), caret::R2(predictions, test.data[,"DA"]))
r2 <- append(r2, get('Rsquare'))
}

################################################################################
#Determine optimal number of components to use via RMSE
################################################################################
#calculate RMSE
rmse.total = rmse.fit[,-1]
r.avg = as.matrix(rowMeans(rmse.total))
r.sd = as.matrix(rowSds(rmse.total))

r.lower = r.avg - r.sd
r.higher = r.avg + r.sd

#plot RMSE vs components to determine optimal number of components to use
x = 1:40
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, r.avg, type = 'p', pch = 16, cex = .75, ylab = 'RMSE',
     xlab = 'Component', xlim = c(1,60), main = 'RMSE by component')
arrows(x, r.lower, x, r.higher,length=0.05, angle=90, code=3)
abline(v = 23, col = 'blue')
abline(h = min(r.higher), col = "Red")
legend('topright', legend = c('Mean', 'Minimum RMSE','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

################################################################################
#Calculate error and r-squared
################################################################################

avg.rmse = mean(error.list)
sd.rmse = sd(error.list)
avg.rmse
sd.rmse 

avg.r2 = mean(r2)
sd.r2 = sd(r2)
avg.r2
sd.r2

################################################################################
#Plot one iteration of PLS for visualization 
################################################################################

training.samples <- 
  createDataPartition(spec_df$DA, p = 0.7, list = FALSE)

train.data  <- spec_df[training.samples, ]
test.data <- spec_df[-training.samples, ]
train.plot <-spec_df2[training.samples, ]
test.plot <- spec_df2[-training.samples, ]
test.plot$color = "black"
test.plot$color[test.plot$Species_ID == "DO"]="red"
test.plot$color[test.plot$Species_ID == "DA"]="blue"

test.hyp = test.plot[test.plot$Species_ID == 'DX',]
remove.in.hyb = setdiff(test.data, test.hyp[,1:2003])
test.hyb2 = test.data[!rownames(test.data) %in% remove.in.hyb,]



plsFit <- train(
  DA~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("repeatedcv", number = 10, repeats = 3),
  tuneLength = 23
)


# Make predictions
predictions <- plsFit %>% predict(test.data)

#Error and r-squred
rmse = caret::RMSE(predictions, test.data[,"DA"])
Rsquare = caret::R2(predictions, test.data[,"DA"])
rmse
Rsquare

#plot actual ancestry vs predicted ancestry 
dev.new(width = 6, height = 6, unit = 'in')
par(mfrow = c(1,1))
plot(test.plot$DA, predictions, cex.lab = 1.5, xlab = "Actual DA ancestry", 
     ylab = "Predicted DA ancestry", 
     xlim = c(-.3,1.3), ylim = c(-.3, 1.3))
lines(x = c(-2, 2), y = c(-2, 2))
points(test.plot$DA, predictions, col = test.plot$color, pch = 16)
legend("bottomright", inset = 0.01,
       legend=c("D. octopetala", "Hybrid", "D. alaskensis", "1:1 line"), 
       text.font = c(3,1,3,1),
       col=c("red", "black", "blue", "black"), pch = c(16,16,16,NA),
       lty = c(NA,NA,NA,1))
