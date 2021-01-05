name = c('tm45','esa3','esa6','wdb6','tm44','wdb41','wdb45','wdb27','wdb46')
anc = c(.52,.69,.53,.61,.47,.48,.58,.43,.36)
rank = c(1,2,2,2,3,3,4,4,5)

df = cbind(name, anc)
df = cbind(df, rank)

cor.test(anc, rank, data = df, method = 'spearman')
