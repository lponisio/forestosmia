## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(vegan)
source('src/misc.R')
source('src/pcoa.R')

load("../data/parasite.Rdata")

hist(parasite$standintensity, breaks=30)

parasite$StandRank <- "low"
parasite$StandRank[parasite$standintensity >= 7 &
                   parasite$standintensity <= 12] <- "medium"

parasite$StandRank[parasite$standintensity >= 13 &
                   parasite$standintensity <= 18] <- "high"

parasite$StandRank[parasite$standintensity >= 19] <- "very high"


infections <- c("Apicystis", "Crithidia", "Ascophaera")

# Does parasite community differ between stands? between stands of different stand intensities?
Stand <- parasite$Stand


## by stand intensity
parasite.comms <- calcPcoa(parasite, infections, nperm=1000,
                           parasite$StandRank)

parasite.comms$tests
# plotting
plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             "parasite_stand_intensity")


## by site
parasite.comms.site <- calcPcoa(parasite, infections, nperm=1000,
                           parasite$Stand)

parasite.comms.site$tests
#$ plotting
plotCommDist(parasite.comms.site$dist$dist, parasite.comms.site$dist$sites,
             "parasite_site")



