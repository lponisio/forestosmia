## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(vegan)
source('src/misc.R')
source('src/pcoa.R')

load("../data/parasite.Rdata")

infections <- c("Apicystis", "Crithidia", "Ascophaera")

# Does parasite community differ between individuals? between stands?
Individual <- parasite$UniqueID
Stand <- parasite$Stand

parasite.comms <- calcPcoa(parasite, infections, nperm=1000, Individual,
                           Stand)
parasite.comms$tests


#$ plotting
plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             parasite.comms$dist$genus, "parasite")









