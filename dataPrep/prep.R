rm(list=ls())
library(dplyr)
library(lubridate)
## prepares raw data and creates dataset for analyses

setwd("~/Dropbox/forestosmia_saved")
dir <- "~/Dropbox/forestosmia_saved"
parasite <- read.csv(file.path(dir,
                               "data/parasite.csv"),
                     stringsAsFactors=FALSE)

dbh <- read.csv(file.path(dir,
                          "data/dbh.csv"),
                stringsAsFactors=FALSE, na.strings=c("","NA"))

floral <- read.csv(file.path(dir,"data/floral2.csv"), stringsAsFactors=FALSE)


repro <- read.csv(file.path(dir,
                            "data/reproductive.csv"),
                  stringsAsFactors=FALSE)

standinfo <- read.csv(file.path(dir,
                                "data/standinfo.csv"),
                      stringsAsFactors=FALSE)

vegcover <- read.csv(file.path(dir,
                               "data/vegcover.csv"),
                     stringsAsFactors=FALSE)

bee <- read.csv(file.path(dir,
                          "data/bee-all.csv"),
                stringsAsFactors=FALSE)

foraging <- read.csv(file.path(dir,
                          "data/osmiaforaging.csv"),
                stringsAsFactors=FALSE)


## ************CLEAN PARASITE DATA AND ADD SUMMARY METRICS**************************************
source("../forestosmia/dataPrep/src/siteNames.R")

dim(parasite)

## drop the NOTES column
parasite$NOTES <- NULL

## drop all the no template controls under UniqueID because they were
## all neg for parasites
parasite <- parasite[!parasite$UniqueID == "C",]
dim(parasite)

## drop any sample where Apidae isnt 1
parasite <- parasite[!parasite$ApidaeCtrl == 0,]
dim(parasite)

## We may want to compare parasites in unemerged and foraging bees, in
## which case use "parasitecontrol" object. However, for most
## analyses, we want to lose the unemerged, control bees, use
## "parasite"

parasitecontrol <- parasite
parasite <- parasitecontrol[!parasitecontrol$Stand == "control",]
dim(parasite)

## community health metrics
parasitenames <- c("Crithidia", "Ascophaera", "Apicystis")

parasite$ParasiteRichness <- rowSums(parasite[, parasitenames])

#### trues/falses times 1 (true is a 1 and false is 0)
parasite$AnyParasite <- (parasite$ParasiteRichness > 0)*1

## calculate avg sick individuals for each site

sick.totals <- parasite %>%
  group_by(Stand) %>%
  summarise(TestedTotals = length(UniqueID),
            ParasitismRate=mean(AnyParasite, na.rm=TRUE),
            InfectedIndividuals=sum(AnyParasite, na.rm=TRUE))

## calculate crithidia, apicystis, and asospheara for each site

sick.parasites <- parasite %>%
  group_by(Stand) %>%
  summarise(InfectedCrith=sum(Crithidia, na.rm=TRUE),
            InfectedApicys=sum(Apicystis, na.rm=TRUE),
            InfectedAsco=sum(Ascophaera, na.rm=TRUE),
            CrithRate=mean(Crithidia, na.rm=TRUE),
            AscoRate=mean(Ascophaera, na.rm=TRUE),
            ApicysRate=mean(Apicystis, na.rm=TRUE))

# create a stand-level dataset and add parasite averages to it
standinfo <- merge(standinfo, sick.parasites, all.x=TRUE)
standinfo <- merge(standinfo, sick.totals, all.x=TRUE)


# create boxplot of parasitism at females at a site



## ***********************************************************
## CLEAN VEG/STAND COVER DATA
## ***********************************************************

## calculate stand-level broadleaf cover (stand=site).  keep only data
## for "broadleaf" trees, srubs, and forbs.  when dropping a variable,
## check the class with class(vegcover$broadleaf) and if its an
## integar dont put quotes around 0.

vegcover <- vegcover[!vegcover$Broadleaf == 0,]

## replace "tr" (tracce) with 1.00
vegcover$perCover[vegcover$perCover == "tr" |
                    vegcover$perCover == "TR"|
                    vegcover$perCover ==  "Tr"]  <- "1.0"

## average abundance of broadleaf cover among plots at a site

vegcover$perCover <- as.numeric(vegcover$perCover)

BLcover.stand <- tapply(vegcover$perCover,
                        vegcover$Stand,
                        mean, na.rm=TRUE)

## add BL Cover as a new column on our stand-level database
standinfo$BLcover <- BLcover.stand[match(standinfo$Stand,
                                         names(BLcover.stand))]

## ***********************************************************
##  CLEAN DBH DATA
## ***********************************************************

## let's remove the columns that are not correctly filled out
dbh$LiveorDead <- NULL
dbh$SlopePerc <- NULL
dbh$Aspect <- NULL

## Right now blanks are assigned as NA, but a blank TreeNum means no
## trees, which we actually want to account for.  If there's an NA
## under TreeNum and DBH, make it zeroes

dbh$TreeNum[is.na(dbh$TreeNum)] <- 0
dbh$DBH[is.na(dbh$DBH)] <- 0

## we want tree species richness, tree species abundance, tree DBH
## averaged across trees at a site

dbh$DBH <- as.numeric(dbh$DBH)
dbh$TreeNum <- as.numeric(dbh$TreeNum)

stand.sum <- dbh %>%
  group_by(Stand) %>%
  summarise(TreeRichness = length(unique(Species)),
            TreeAbund=mean(TreeNum, na.rm=TRUE),
            MeanDBH=mean(DBH, na.rm=TRUE))

## then merge to dataset
standinfo <- merge(standinfo, stand.sum)

## ***********************************************************
## CLEAN FLORAL DATA, MERGE
## ***********************************************************

## This dataset includes floral data from first two sampling rounds, when we had osmia out\

## we want abundance flowers (abundflw), richness flowering species
## (richnessflwingsp), and abundfloweringsp

## let's drop the column "vials" since it's not correctly filled out
floral$Vials <- NULL
floral$Other_flowers <- NULL
floral$Notes <- NULL

# Floral Richness
# we don't want zeroes counted as 1 though
floral$Flower_sci[floral$Flower_sci == 0] <- NA

floral.sum <- floral %>%
  group_by(Stand, Transect) %>%
  summarise(FlowerRichness =
              length(unique(Flower_sci[!is.na(Flower_sci)])),
            BloomAbund=sum(Blooms, na.rm=TRUE),
            StemAbund=sum(Stems, na.rm=TRUE),
            CanopyCover=mean(Canopy_cent,  na.rm=TRUE))

mean.floral <- floral.sum %>%
  group_by(Stand) %>%
  summarise(FlowerRichness = mean(FlowerRichness, na.rm=TRUE),
            MeanBloomAbund=mean(BloomAbund, na.rm=TRUE),
            MeanStemAbund=mean(StemAbund, na.rm=TRUE),
            MeanCanopy=mean(CanopyCover,  na.rm=TRUE))

standinfo <- merge(standinfo, mean.floral, all.x=TRUE)

## ***********************************************************
## BEE DIVERSITY DATA, MERGE
## ***********************************************************

## I decided to calculate bee diversity over 2 sample rounds

# drop any sample from round 3
bee <- bee[!bee$Round == 3,]

bee$Caste <- NULL

bee$GenSp[bee$GenSp=="0 0"] <- NA

bee$GenSp[bee$GenSp=="NA NA"] <- NA

## generate abundance of bees in each sample round
## calculate bee abundance with just net data
## calculate bee richness with net + pan data
bee.sum <- bee %>%
  group_by(Stand, Round, Trap.type, Year) %>%
  summarise(BeeRichness = length(unique(GenSp[!is.na(GenSp)])),
            BeeAbund=sum(Count, na.rm=TRUE))

mean.bee.rich <- bee.sum %>%
  group_by(Stand, Year) %>%
  summarise(MeanBeeRichness = mean(BeeRichness, na.rm=TRUE))


mean.bee.abund <- bee.sum %>%
 group_by(Stand, Trap.type, Year) %>%
    summarise(MeanBeeAbund=mean(BeeAbund, na.rm=TRUE))

# We only want bees that were netting and bees from 2019
mean.bee.abund <- mean.bee.abund[mean.bee.abund$Trap.type == "Net" &
                        mean.bee.abund$Year == "2019",]


mean.bee.rich <- mean.bee.rich[mean.bee.rich$Year == "2019",]

# remove columsn we don't want
mean.bee.abund$Trap.type <- NULL
mean.bee.abund$Year <- NULL
mean.bee.rich$Year <- NULL

# merge to site-level data
standinfo <- merge(standinfo, mean.bee.abund, all.x=TRUE)
standinfo <- merge(standinfo, mean.bee.rich, all.x=TRUE)


## ***********************************************************
## honey bees
## ***********************************************************

#honey bees from net data

hb.bee <- bee[bee$Genus == "Apis",]

hb.sum <- hb.bee %>%
  group_by(Stand, Round, Trap.type, Year) %>%
  summarise(BeeAbund=sum(Count, na.rm=TRUE))

mean.hb <- hb.sum %>%
  group_by(Stand, Trap.type, Year) %>%
  summarise(MeanHBAbund=mean(BeeAbund, na.rm=TRUE))


#We only want hb that were netting and bees from 2019
mean.hb <- mean.hb[mean.hb$Trap.type == "Net" &
                                  mean.hb$Year == "2019",]

mean.hb$Trap.type <- NULL
mean.hb$Year <- NULL

standinfo <- merge(standinfo, mean.hb, all.x=TRUE)

# we are done with the stand-level dataset, let's rename it
site.data <- standinfo

## ***********************************************************
## CLEAN REPRODUCTIVE DATA, MERGE
## ***********************************************************

## we are just using females from Osmia nests in blocks
## where samples were not pulled for parasite tests

## TUBE LEVEL DATA ********

## these columns not filled out
repro$NestNum_Bri <- NULL
repro$Notes <- NULL

## not correctly calculated
repro$Total <- NULL

## drop any data from the Pathogen study
repro <- repro[!repro$box_type == "Path",]

## confirm that only data is from reproductive study, "FPH" blocks
unique(repro$box_type)

## rename "nest number" as column
colnames(repro)[colnames(repro) == "NestNum_XRAY"] <- "Column"

## to get average number of Females/Males at a stand, by nest number, 
## we want to give each tube a unique ID, called "BlockNestTube"
repro$BlockNestTube  <- paste0(repro$Block, repro$Row, repro$Column)
parasite$BlockNestTube  <- paste0(parasite$Block,
                                  parasite$Nest)

repro$SumOffspring <- repro$Females + repro$Males

## BLOCK LEVEL DATA ********

## for each stand, average the number of females and females at each block
## end up with Block A averages and Block B averages for each site

repro.block <- aggregate(list(Females=repro$Females,
                              Males=repro$Males),
                         list(Stand=repro$Stand,
                              Block=repro$Block),
                         sum, na.rm=TRUE)

## add more offspring summary data to block-level averages
## above 1 F > M, below 1 M>F, = 1 M=F
repro.block$FM_ratio <-  log(repro.block$Females +1)/log(repro.block$Males +1)

repro.block$SumOffspring <-  repro.block$Females + repro.block$Males


## create a new block level dataset to describe reprodata, with stand info merged
repro.block <- merge(repro.block, standinfo)

## merge site level sick totals to our reproductive, block-level dataset
repro.block <- merge(repro.block, sick.totals)
repro.block <- merge(repro.block, sick.parasites)

## INDIVIDUAL LEVEL DATASET ********

## the parasite dataset is our individual level data,
## let's merge tube level data to this dataset and rename it as our indiv dataset

parasite$Females <- repro$Females[match(
  paste(parasite$BlockNestTube,
        parasite$Stand),
  paste(repro$BlockNestTube,
        repro$Stand))]

parasite$Males <- repro$Males[match(
  paste(parasite$BlockNestTube,
        parasite$Stand),
  paste(repro$BlockNestTube,
        repro$Stand))]

parasite$SumOffspring <- repro$SumOffspring[match(
  paste(parasite$BlockNestTube,
        parasite$Stand),
  paste(repro$BlockNestTube,
        repro$Stand))]

## merge stand info onto parasite (individual data) and rename it as indivi. data
dim(parasite)
indiv.data <- merge(parasite, standinfo)
dim(indiv.data)

## ***********************************************************
## CLEAN FORAGING DATA, MERGE
## ***********************************************************

## we are just looking at foraging distance of females from Osmia nests in blocks
## where samples were not pulled for parasite tests

## these columns not filled out entirely, have blanks
foraging$videosanalyzedby <- NULL
foraging$collectorinitials <- NULL
foraging$tempF <- NULL
foraging$NumCappedCells <- NULL
foraging$videostart.h.min <- NULL
foraging$Activity <- NULL
foraging$Notes <- NULL
foraging$BeeLeaves.h.mm.ss <- NULL
foraging$BeeReturns.h.mm.ss <- NULL

## drop any data where trip length could not be calculated due to lack of footage

foraging$TripLength.h.mm.ss[foraging$TripLength.h.mm.ss==""] <- NA
foraging$TripLength.h.mm.ss[foraging$TripLength.h.mm.ss==" "] <- NA
foraging <- na.omit(foraging)


## BLOCK LEVEL DATA ********

## for each stand, average the trip length of females at each block
## end up with Block A averages and Block B averages for each site

## conver triplenth from character to time with hms in lubridate package
foraging$TripLength.h.mm.ss <- hms(foraging$TripLength.h.mm.ss)

foraging.block <- foraging %>%
  group_by(Stand, NestBlock) %>%
  summarise(MeanTripLengthBlock=mean(TripLength.h.mm.ss, na.rm = TRUE))


## merge to our block-level dataset
repro.block <- merge(repro.block, foraging.block)


## SITE/INDIVIDUAL LEVEL DATASET ********

## let's create site-level data and merge onto site and individual level datasets

foraging.stand <- foraging %>%
  group_by(Stand) %>%
  summarise(MeanTripLengthStand=mean(TripLength.h.mm.ss, na.rm = TRUE))

## merge stand-level
standinfo <- merge(standinfo, foraging.stand, all.x=TRUE)

## merge individual data
dim(indiv.data)
indiv.data <- merge(indiv.data, foraging.stand)
dim(indiv.data)

## *************************************************************
## write out final datasets: stand-level, block-level, indiv-level
## *************************************************************

dir <- "~/Dropbox/forestosmia_saved/cleaneddata"
#site.data <- merge(standinfo, sick.totals, all.x=TRUE)

write.csv(indiv.data, file=file.path(dir,
                                     "indivData.csv"), row.names=FALSE)

write.csv(standinfo, file=file.path(dir,
                                    "standinfo.csv"), row.names=FALSE)

write.csv(repro.block, file=file.path(dir,
                                      "NestRepro.csv"), row.names=FALSE)


save.dir <- "~/Dropbox/forestosmia/data"
save(indiv.data,
     file=file.path(save.dir, "indiv.data.Rdata"))

save(repro.block,
     file=file.path(save.dir, "repro.block.Rdata"))

save(site.data,
     file=file.path(save.dir, "site.data.Rdata"))

