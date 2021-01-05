#Caret PLSDA
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)
library(rlist)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Fit PLS_DA model all dry
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
spec_all = spec_all[!meta(spec_all)$Location == "NaN",]


spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Location)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Location"] <- "Location"

#Set number of components to be used
ncomp = 32

#Partition Data
accuracy <- c()
kappa <- c()
k.fit <- matrix(nrow = ncomp)
cm.list <- list()
loadings.c1 <- matrix(nrow = 201)
loadings.c2 <- matrix(nrow = 201)
loadings.c3 <- matrix(nrow = 201)

for(i in 1:100){

inTrain <- caret::createDataPartition(
  y = spec_df$Location,
  p = .7,
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
  Location ~ .,
  data = training,
  maxit = 10000,
  method = "pls",
  trControl = ctrl,
  tuneLength = ncomp)



#objects for determining n components
k = assign(paste0('k', i), as.matrix(plsFit$results$Kappa))
k.fit <- cbind(k.fit, get('k'))

#loadings
c1 = assign(paste0('c1',i), as.matrix(plsFit$finalModel$loadings[,1]))
loadings.c1 <- cbind(loadings.c1, get('c1'))

c2 = assign(paste0('c2',i), as.matrix(plsFit$finalModel$loadings[,2]))
loadings.c2 <- cbind(loadings.c2, get('c2'))

c3 = assign(paste0('c3',i), as.matrix(plsFit$finalModel$loadings[,3]))
loadings.c3 <- cbind(loadings.c3, get('c3'))

#test model
plsClasses <- predict(plsFit, newdata = testing)

#objects to assess accuracy 
cm = confusionMatrix(data = plsClasses, as.factor(testing$Location))
cm.m = assign(paste0("cm", i), as.matrix(cm))
cm.list <- list.append(cm.list, get('cm.m'))

ac <- assign(paste0('acc',i), cm$overall[1])
accuracy <- append(accuracy, get('ac'))

kap = assign(paste0("kap",i), cm$overall[2])
kappa <- append(kappa, get('kap'))

}

#accuracy and kappa
mean.acc = mean(accuracy)
sd.acc = sd(accuracy)

mean.kap = mean(kappa)
sd.kap = sd(kappa)

mean.acc
sd.acc

mean.kap
sd.kap

#kappa (for choosing components)

k.total = k.fit[,-1]
kavg = as.matrix(rowMeans(k.total))
ksd = as.matrix(rowSds(k.total))

klower = kavg - ksd
khigher = kavg + ksd
x = 1:60
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, kavg, type = 'p', pch = 16, cex = .75, ylab = 'Kappa', xlab = 'Component', 
     xlim = c(1,60), main = 'Kappa for Location')
arrows(x, klower, x, khigher,length=0.05, angle=90, code=3)
abline(v = 38, col = 'blue')
abline(h = max(klower), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum kappa','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))


#plot confusion matrix
cm.total = Reduce('+', cm.list)/100
cm.total = t(cm.total)
cm.total = cm.total/rowSums(cm.total)

cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('DA', 'DO', 'DX')
colnames(cm.total) <- c('DA', 'DO', 'DX')


write.csv(cm.total, "Figures/cm_final/cm_Location_mean.csv")

#sp loc special code
cm.total = read.csv("Figures/cm_final/cm_sp_loc_mean.csv", stringsAsFactors = T)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- cm.total[,1]
cm.total = cm.total[,-1]
cm.total = mapply(cm.total, FUN = as.numeric)
cm.total = matrix(data = cm.total, ncol = 12, nrow = 12)
rownames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 'MD', 'WDA', 'WDB', 'ES', 'TM', 'WDB')

colnames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 'MD', 'WDA', 'WDB', 'ES', 'TM', 'WDB')

#plot
pdf(file= "Figures/cm_final/dry/test.pdf", width = 6, height = 6)

dev.new(width = 6, height = 8, unit = 'in')

cols = colorRampPalette(c('#f5f5f5', '#b35806'))

par(mar = c(1,2,4,1), oma = c(1,1,3,1))
corrplot(cm.total, is.corr = T, method = 'square', col = cols(10), addCoef.col = '#542788',
         tl.srt = 0, tl.offset = 1, number.digits = 2, tl.cex = 1.2, 
         cl.cex = 1, number.cex = 1.5,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Reference", side = 2, line = 2, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 3.5, line = 6)


#loadings
comp_to_spec = function(x){
  t.comp = t(x)
  colnames(t.comp) <- seq(400,2400, by = 10)
  s.comp = as_spectra(t.comp)
}


#plot loadings

component1 = comp_to_spec(loadings.c1)
component2 = comp_to_spec(loadings.c2)
component3 = comp_to_spec(loadings.c3)


par(mar = c(5,4,1,1), oma = c(1,1,1,1), mfrow = c(1,1))
plot(mean(component1), lwd = 2, lty = 1, col = rgb(1,0,0,1), 
     cex.lab = 1.5, ylim = c(-.2, .15), ylab = "Loading Values", 
     xlab = "Wavelength (nm)")
plot_quantile(component1, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(component1, regions = default_spec_regions(), add = TRUE)
plot(mean(component2), lwd = 1.5, lty = 1, col = rgb(0,0,1,1), add = TRUE)
plot_quantile(component2, total_prob = 0.95, col = rgb(0, 0, 1, 0.25), border = FALSE, add = TRUE)

abline(h = 0, lty = 2, lwd = 1.5)
legend('bottomright',inset = .02, legend=c("Component 1", "Component 2", 'Component 3'),
              col=c(rgb(1,0,0,1), rgb(0,0,1,1), rgb(0,1,0,1)), lty=1, cex=0.8, bg ='white')


plot(mean(component3), lwd = 1.5, lty = 1, col = "darkgreen", add = TRUE)
plot_quantile(component3, total_prob = 0.95, col = rgb(0, .5, 0, 0.25), border = FALSE, add = TRUE)

#hybrid stuff
hybrid.spec = spec_all[meta(spec_all)$Location == "DX",]
hybrids = as.matrix(hybrid.spec)
hyb.meta = meta(hybrid.spec)

set.seed(7)
pred_hyb = predict(plsFit, newdata = hybrids)
hyb.df = as.data.frame(pred_hyb)
View(hyb.df)
hyb.df = cbind(hyb.df, hyb.meta)
write.csv(hyb.df, "Figures/hybrid_cms/hybrid_pop_predictions_all.csv")
hyb_cm = read.csv('Figures/hybrid_cms/hybrid matrices/hybrid_cm_all.csv', 
                  stringsAsFactors = F)
rownames(hyb_cm) <- hyb_cm$X
hyb_cm = hyb_cm[,-1]
hyb_cm = hyb_cm %>% replace_with_na_all(condition = ~.x == 0)
hyb_cm = as.matrix(hyb_cm)
rownames(hyb_cm) <- c('es', 'tm', 'wdb')

par(mar = c(0,0,0,0), oma = c(0,0,3,0))
corrplot(hyb_cm, is.corr = F, method = 'color', addCoef.col = 'darkorange2',
         tl.srt = 0, tl.offset = 1.5, number.digits = 2, tl.cex = 1,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Collection Location", side = 2, line = -10, cex = 1.5)
mtext("Predicted Population", side = 3, cex = 1.5, at = 2, line = 1)






