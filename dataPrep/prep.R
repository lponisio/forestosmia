rm(list=ls())
library(dplyr)
## prepares raw data and creates dataset for analyses

setwd("~/Dropbox/forestosmia_saved")
dir <- "~/Dropbox/forestosmia_saved"
parasite <- read.csv(file.path(dir,
                           "data/parasite.csv"),
                 stringsAsFactors=FALSE)

dbh <- read.csv(file.path(dir,
                               "data/dbh.csv"),
                     stringsAsFactors=FALSE, na.strings=c("","NA"))

floral <- read.csv(file.path(dir,"data/floral.csv"), stringsAsFactors=FALSE)


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

standinfo <- merge(standinfo, sick.parasites)


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

## merge  stand and stand intensity
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

## TBA WHICH ESTIMATE SHOULD BE USE FOR THE NESTS/PARASITES, THE MEAN
## OVER THE SEASON OR THE CLOSEST SAMPLE ROUND, OR THE AVERAGE OVER
## THE FRIST 2

## drop any sample where Apidae isnt 1
bee <- bee[!bee$Round == 3,]

bee$Caste <- NULL

bee$GenSp[bee$GenSp=="0 0"] <- NA

bee$GenSp[bee$GenSp=="NA NA"] <- NA

## generate abundance of bees in each sample round
bee.sum <- bee %>%
  group_by(Stand, Round, Trap.type, Year) %>%
  summarise(BeeRichness = length(unique(GenSp[!is.na(GenSp)])),
            BeeAbund=sum(Count, na.rm=TRUE))

mean.bee <- bee.sum %>%
  group_by(Stand, Trap.type, Year) %>%
  summarise(MeanBeeRichness = mean(BeeRichness, na.rm=TRUE),
            MeanBeeAbund=mean(BeeAbund, na.rm=TRUE))

mean.bee.net <- mean.bee[mean.bee$Trap.type == "Net" &
                        mean.bee$Year == "2019",]

mean.bee.net$Trap.type <- NULL
mean.bee.net$Year <- NULL

standinfo <- merge(standinfo, mean.bee.net, all.x=TRUE)

## ***********************************************************
## CLEAN REPRODUCTIVE DATA, MERGE
## ***********************************************************

## average number of Females/Males at a stand, by nest number
colnames(repro)[colnames(repro) == "NestNum"] <- "Column"

repro$BlockNestTube  <- paste0(repro$Block, repro$Row, repro$Column)
parasite$BlockNestTube  <- paste0(parasite$Block,
                                     parasite$Nest)

repro$SumOffspring <- repro$Females + repro$Males

## BLOCK LEVEL DATA ********

repro.nest <- aggregate(list(Females=repro$Females,
                             Males=repro$Males),
                       list(Stand=repro$Stand,
                            Block=repro$Block),
                       sum, na.rm=TRUE)

## above 1 F > M, below 1 M>F, = 1 M=F
repro.nest$FM_ratio <-  log(repro.nest$Females +1)/log(repro.nest$Males +1)

repro.nest$SumOffspring <-  repro.nest$Females + repro.nest$Males


## block level data with stand info merged

repro.nest <- merge(repro.nest, standinfo)

## merge site level sick totals
repro.nest <- merge(repro.nest, sick.totals)
repro.nest <- merge(repro.nest, sick.parasites)

## INDIVIDUAL LEVEL DATA ********
## an merge tube level data to the patasite data mom

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

## merge parasite and stand info
dim(parasite)
indiv.data <- merge(parasite, standinfo)
dim(indiv.data)




## *************************************************************
## write out final data
dir <- "~/Dropbox/forestosmia_saved/cleaneddata"
site.data <- merge(standinfo, sick.totals, all.x=TRUE)

write.csv(indiv.data, file=file.path(dir,
                     "indivData.csv"), row.names=FALSE)

write.csv(standinfo, file=file.path(dir,
                     "standinfo.csv"), row.names=FALSE)

write.csv(repro.nest, file=file.path(dir,
                       "NestRepro.csv"), row.names=FALSE)


save.dir <- "~/Dropbox/forestosmia/data"
save(indiv.data,
     file=file.path(save.dir, "indivData.Rdata"))

save(repro.nest,
     file=file.path(save.dir, "NestRepro.Rdata"))

save(site.data,
     file=file.path(save.dir, "siteData.Rdata"))



#for analysis, do path analysis
#where is bee richness?
# stand itensity (broadlaef), dbh, flowers, bee richness -> parasitsm (any parasite, dont use parasite richness) -> reprodcutive output (fecundity of males and females, just try with females)
# stand itensity
# dont try to put the path analyses yet, do one glmm for parasitism, and do one glmm for reproduction, with parasitism in the model for reproduction
