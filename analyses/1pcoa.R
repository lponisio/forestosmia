setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')

setwd("analyses")
rm(list=ls())
library(vegan)
source('src/misc.R')
source('src/pcoa.R')

# *************************************************************************
# Make a PCOA Based on Broadleaf cover, a measure of management intensity
# *************************************************************************

load("../data/indivdata.Rdata")

hist(indiv.data$BLcover, breaks=30)

infections <- c("Apicystis", "Crithidia", "Ascophaera")
Stand <- indiv.data$Stand

## PCOA by stand intensity
parasite.comms <- calcPcoa(indiv.data, infections, nperm=1000,
                           indiv.data$StandIntensity)
parasite.comms$tests

# plotting
plotCommDist(parasite.comms$dist$dist, parasite.comms$dist$sites,
             "parasite_stand_intensity2")


# *************************************************************************
# Make a PCOA Based on Site
# *************************************************************************
## by site
parasite.comms.site <- calcPcoa(indiv.data, infections, nperm=1000,
                                indiv.data$Stand)

parasite.comms.site$tests
#$ plotting
plotCommDist(parasite.comms.site$dist$dist, parasite.comms.site$dist$sites,
             "parasite_site2")


# *************************************************************************
# Make a PCOA Based on DBH
# *************************************************************************

load("../data/indivdata.Rdata")

hist(indiv.data$MeanDBH, breaks=30)

indiv.data$StandDBH <- "low"
indiv.data$StandDBH[indiv.data$MeanDBH >= 6 &
                      indiv.data$MeanDBH <= 10] <- "medium"

indiv.data$StandDBH[indiv.data$MeanDBH >= 11 &
                      indiv.data$MeanDBH <= 15] <- "high"

indiv.data$StandDBH[indiv.data$MeanDBH >= 16] <- "very high"


infections <- c("Apicystis", "Crithidia", "Ascophaera")

Stand <- indiv.data$Stand

## PCOA by stand dbh
parasite.comms.dbh <- calcPcoa(indiv.data, infections, nperm=1000,
                               indiv.data$StandDBH)

parasite.comms.dbh$tests

#$ plotting
plotCommDist(parasite.comms.dbh$dist$dist, parasite.comms.dbh$dist$sites,
             "parasiteDBH")

