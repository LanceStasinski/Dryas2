#Plot ancestry vs hybrid population prediction
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
# Fit PLS_DA model all dry
################################################################################

#data
hyb = read.csv("Figures/hybrid_cms/hybrid_pop_predictions_all.csv", stringsAsFactors = F)


par(mar = c(5,5,4,4))
boxplot(DA~Location, data = hyb, ylab = 'Proportion of D. alaskensis ancestry',
        xlab = 'Location', ylim = c(.3,.7), notch = F,
        main = 'D. alaskensis ancestry of hybrids versus location')
abline(h = 0.5, lty = 2, lwd = 1)
