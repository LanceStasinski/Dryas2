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


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

oct = spec_all[meta(spec_all)$Species_ID == "DO"]
ala = spec_all[meta(spec_all)$Species_ID == "DA"]
hby = spec_all[meta(spec_all)$Species_ID == "DX"]
################################################################################
# Fit PLS_DA model all dry
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
hybrid.spec = spec_all[meta(spec_all)$Species_ID == "DX",]
hybrids = as.matrix(hybrid.spec)
hyb.meta = meta(hybrid.spec)

spec_all = spec_all[!meta(spec_all)$Species_ID == "NaN",]
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

for(i in 1:10){

set.seed(i)

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
cm = confusionMatrix(data = plsClasses, testing$Species_ID)



cm.m = as.matrix(cm)

assign(paste0("cm", i), cm.m)
}

#hybrid stuff
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

#plot confusion matrix

cm.total = (cm1 + cm2 + cm3 + cm4 + cm5 + cm6 + cm7 + cm8+ cm9 + cm10)/10
cm.total = t(cm.total)
write.csv(cm.total, "Figures/raw confusion matrices/Species_ID_10it_20")
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('DA', 'DO', 'DX')

write.csv(cm.total, "Figures/raw confusion matrices/Species_ID_test.csv")

cm.total = read.csv("Figures/raw confusion matrices/Species_ID_test.csv", stringsAsFactors = F)
row.names(cm.total) <- cm.total[,1]
cm.total = cm.total[,-1]
cm.total = as.matrix(cm.total)

par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
corrplot(cm.total, is.corr = F, method = 'color', addCoef.col = 'darkorange2',
         tl.srt = 0, tl.offset = 1.5, number.digits = 2, tl.cex = .75,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Reference", side = 2, line = -5, cex = 1.5)
mtext("Prediction", side = 3, cex = 1.5, at = 2, line = 2)



#loadings

component1 = Reduce(cbind, 
                    list(comp1_1, comp1_2, comp1_3, comp1_4, comp1_5, comp1_6, 
                         comp1_7, comp1_8, comp1_9, comp1_10))

component2 = Reduce(cbind, 
                    list(comp2_1, comp2_2, comp2_3, comp2_4, comp2_5, comp2_6, 
                         comp2_7, comp2_8, comp2_9, comp2_10))

component3 = Reduce(cbind, 
                    list(comp3_1, comp3_2, comp3_3, comp3_4, comp3_5, comp3_6, 
                         comp3_7, comp3_8, comp3_9, comp3_10))


#plot loadings

comp_to_spec = function(x){
  t.comp = t(x)
  colnames(t.comp) <- seq(400,2400, by = 10)
  s.comp = as.spectra(t.comp)
}

component1 = comp_to_spec(component1)
component2 = comp_to_spec(component2)
component3 = comp_to_spec(component3)

par(mar = c(4,4,1,1), oma = c(1,1,1,1))
plot(mean(component1), lwd = 2, lty = 1, col = rgb(1,0,0,1), 
     cex.lab = 1.2, ylim = c(-.2, .15), ylab = "Loading Values", 
     xlab = "Wavelength")
plot_quantile(component1, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(component1, regions = default_spec_regions(), add = TRUE)
plot(mean(component2), lwd = 1.5, lty = 1, col = rgb(0,0,1,1), add = TRUE)
plot_quantile(component2, total_prob = 0.95, col = rgb(0, 0, 1, 0.25), border = FALSE, add = TRUE)
plot(mean(component3), lwd = 1.5, lty = 1, col = "darkgreen", add = TRUE)
plot_quantile(component3, total_prob = 0.95, col = rgb(0, .5, 0, 0.25), border = FALSE, add = TRUE)
abline(h = 0, lty = 2, lwd = 1.5)
legend('bottomright',inset = .02, legend=c("Component 1", "Component 2", "Component 3"),
              col=c(rgb(1,0,0,1), rgb(0,0,1,1), "darkgreen"), lty=1, cex=0.8)










