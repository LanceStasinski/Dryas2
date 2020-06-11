library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
f = function(spec, specnumber){
  spec[specnumber,]
}

################################################################################
#Dry leaves
################################################################################

#low
l1 = f(spec_all, 9) #ala, tm
l2 = f(spec_all, 12) #ala, tm
l3 = f(spec_all, 10) #ala, tm
l4 = f(spec_all, 32)
l5 = f(spec_all, 51)
l6 = f(spec_all, 79)
l7 = f(spec_all, 82)
l8 = f(spec_all, 70)
l9 = f(spec_all, 39)
l10 = f(spec_all, 81)
l11 = f(spec_all, 663)
l12 = f(spec_all, 750)
l13 = f(spec_all, 743)
l14 = f(spec_all, 683)
l15 = f(spec_all, 651)
l16 = f(spec_all, 673)
l17 = f(spec_all, 707)
l18 = f(spec_all, 784)
l19 = f(spec_all, 744)
l20 = f(spec_all, 603)

low = meta(Reduce(combine, list(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l12,l13,l14,l15,
                                l16,l17,l18,l19,l20)))
#species: ala, oct, and hyb
#locations: tm, bg, wda, mdb

#medium

m1 = f(spec_all, 583)
m2 = f(spec_all, 582)
m3 = f(spec_all, 89)
m4 = f(spec_all, 572)
m5 = f(spec_all, 575)
m6 = f(spec_all, 560)
m7 = f(spec_all, 97)
m8 = f(spec_all, 143)
m9 = f(spec_all, 95)
m10 = f(spec_all, 17)
m11 = f(spec_all, 497)
m12 = f(spec_all, 722)
m13 = f(spec_all, 682)
m14 = f(spec_all, 619)
m15 = f(spec_all, 733)
m16 = f(spec_all, 662)
m17 = f(spec_all, 766)
m18 = f(spec_all, 489)
m19 = f(spec_all, 644)
m20 = f(spec_all, 367)

  
medium = meta(Reduce(combine, list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,
                                   m14,m15,m16,m17,m18,m19,m20)))
#species: oct and hyb
#locations: es, wdb, bg, tm, mdb

h1 = f(spec_all, 194)
h2 = f(spec_all, 195)
h3 = f(spec_all, 327)
h4 = f(spec_all, 297)
h5 = f(spec_all, 386)
h6 = f(spec_all, 494)
h7 = f(spec_all, 344)
h8 = f(spec_all, 205)
h9 = f(spec_all, 272)
h10 = f(spec_all, 319)
h11 = f(spec_all, 371)
h12 = f(spec_all, 372)
h13 = f(spec_all, 236)
h14 = f(spec_all, 384)
h15 = f(spec_all, 221)
h16 = f(spec_all, 230)
h17 = f(spec_all, 314)
h18 = f(spec_all, 302)
h19 = f(spec_all, 345)
h20 = f(spec_all, 262)

high = meta(Reduce(combine, list(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, 
                                 h12, h13, h14, h15, h16, h17, h18, h19, h20)))
#species: ala, oct, and hyb
#locations: 1 wdb, rest es

plot_interactive(spec_all[200:783,])

################################################################################
#Dry leaves
################################################################################
#low
l1 = f(spec_all, 343) #ala, tm
l2 = f(spec_all, 280) #ala, tm
l3 = f(spec_all, 298) #ala, tm
l4 = f(spec_all, 259)
l5 = f(spec_all, 36)
l6 = f(spec_all, 276)
l7 = f(spec_all, 185)
l8 = f(spec_all, 266)
l9 = f(spec_all, 269)
l10 = f(spec_all, 188)
l11 = f(spec_all, 155)
l12 = f(spec_all, 435)
l13 = f(spec_all, 20)
l14 = f(spec_all, 431)
l15 = f(spec_all, 223)
l16 = f(spec_all, 293)
l17 = f(spec_all, 213)
l18 = f(spec_all, 473)
l19 = f(spec_all, 390)
l20 = f(spec_all, 411)

loww = meta(Reduce(combine, list(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l12,l13,l14,l15,
                                l16,l17,l18,l19,l20)))
#species: ala (mostly), oct
#locations: es (mostly), wdb, tm

#medium

m1 = f(spec_all, 369)
m2 = f(spec_all, 461)
m3 = f(spec_all, 53)
m4 = f(spec_all, 76)
m5 = f(spec_all, 198)
m6 = f(spec_all, 283)
m7 = f(spec_all, 401)
m8 = f(spec_all, 96)
m9 = f(spec_all, 91)
m10 = f(spec_all, 93)
m11 = f(spec_all, 284)
m12 = f(spec_all, 90)
m13 = f(spec_all, 427)
m14 = f(spec_all, 140)
m15 = f(spec_all, 92)
m16 = f(spec_all, 385)
m17 = f(spec_all, 419)
m18 = f(spec_all, 421)
m19 = f(spec_all, 99)
m20 = f(spec_all, 105)


mediumw = meta(Reduce(combine, list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,
                                   m14,m15,m16,m17,m18,m19,m20)))
#species: oct, hyb, ala
#locations: es, tm, wdb

h1 = f(spec_all, 85)
h2 = f(spec_all, 145)
h3 = f(spec_all, 146)
h4 = f(spec_all, 380)
h5 = f(spec_all, 82)
h6 = f(spec_all, 113)
h7 = f(spec_all, 426)
h8 = f(spec_all, 464)
h9 = f(spec_all, 285)
h10 = f(spec_all, 114)
h11 = f(spec_all, 78)
h12 = f(spec_all, 97)
h13 = f(spec_all, 105)
h14 = f(spec_all, 359)
h15 = f(spec_all, 400)
h16 = f(spec_all, 99)
h17 = f(spec_all, 127)
h18 = f(spec_all, 12)
h19 = f(spec_all, 102)
h20 = f(spec_all, 141)

highw = meta(Reduce(combine, list(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, 
                                 h12, h13, h14, h15, h16, h17, h18, h19, h20)))
#species: ala, oct, and hyb
#locations: tm (mostly), es, wdb

plot_interactive(spec_all[200:783,])
