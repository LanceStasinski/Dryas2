
library(corrplot)
library(plyr)
library(vegan)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

leaf = read.csv('morphology.csv', stringsAsFactors = F)



#data sets
traits = leaf[, -c(1,2,3,4,10,11,12)]
parent = traits[!traits$Species == 'DX',]
parent = parent[,-5]
hyb = traits[traits$Species == 'DX',]
hyb = hyb[,-5]


#parent trait correlations
p.c = cor(parent, method = 'pearson')
colnames(p.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')
rownames(p.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')


#hybrid trait correlations
hyb.c = cor(hyb, method = 'pearson')
colnames(hyb.c) <-c('Glands', 'Scales', 'Length', 'Tomentum')
rownames(hyb.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')

#corrplots

cols = colorRampPalette(c('#2166ac', '#d73027'))

par(mfrow = c(1,2))
corrplot::corrplot(p.c, method = "number", type = "lower", diag = F, cl.cex = 1,
                   addCoef.col = 'black', tl.col = 'black', cl.length = 5,
                   col = cols(15), number.cex = 1.25)
mtext('Parent Leaves', side = 2, line = -1, at = 2, cex = 1.5)
corrplot::corrplot(hyb.c, method = "number", type = 'lower', diag = F,
                   cl.cex = 1, addCoef.col = 'black',
                   tl.col = 'black', cl.length = 5, col = cols(15),
                   number.cex = 1.25)
mtext('Hybrid Leaves', side = 2, line = -1, at = 2, cex = 1.5)