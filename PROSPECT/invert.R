detach("package:spectrolab", unload = TRUE)
library("PEcAnRTM")

################################################################################
# Read Data
################################################################################

obs = as.data.frame(readRDS("mat_all_spec_dry.rds"))

################################################################################
################################################################################

f = function(){
    c(N = 1.6, Cab = 30, Car = 8, Cw = 0.007, Cm = 0.005)
}

################################################################################
# Set up gibbs
################################################################################

invert.options              = default.settings.prospect
invert.options$n.tries      = 1
invert.options$nchains      = 4
invert.options$ngibbs       = 50000
invert.options$burnin       = 30000
invert.options$ngibbs.max   = 100000
invert.options$do.lsq.first = FALSE
invert.options$save.samples = FALSE
invert.options$quiet        = FALSE
invert.options$param.maxs   = c(3, 100, 25, 0.05, 0.02)
invert.options$param.mins   = c(1, 0, 0, 0, 0)

invert.options$model        =  function(params, seed = NULL){
    prospect(params, "5")[ 1:2001, 1]
}
invert.options$inits.function = f

################################################################################
# Run
################################################################################

########################################
# First spec
########################################

r = PEcAnRTM::invert.auto(observed        = obs[ , 1],
                          invert.options  = invert.options,
                          return.samples  = FALSE,
                          quiet           = FALSE)$results

r = sapply(r, c)
s = paste(c(colnames(obs)[1], r), collapse = ", ")

file.remove("log.csv")
file.create("log.csv")
cat(paste0(c("id", names(r)), collapse = ", "), sep = "\n",
    file = "log.csv", append = FALSE)
cat(s, sep = "\n", file = "log.csv", append = TRUE)


########################################
# Remainder in parallel
########################################

for(i in 2:length(obs)){
    r = PEcAnRTM::invert.auto(observed        = obs[ , i],
                              invert.options  = invert.options,
                              return.samples  = FALSE,
                              quiet           = FALSE)$results
    r = sapply(r, c)
    s = paste(c(colnames(obs)[i], r), collapse = ", ")

    cat(s, sep = "\n", file = "log.csv", append = TRUE)
}
