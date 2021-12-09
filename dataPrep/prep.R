rm(list=ls())
library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)
library(Taxonstand)
## prepares raw data and creates dataset for analyses

setwd('/Volumes/bombus/Dropbox (University of Oregon)/')

setwd("forestosmia_saved")

parasite <- read.csv("data/parasite.csv",
                     stringsAsFactors=FALSE)

floral <- read.csv("data/floral2.csv", stringsAsFactors=FALSE)

repro <- read.csv("data/reproductive.csv",
                  stringsAsFactors=FALSE)

standinfo <- read.csv("data/standinfo.csv",
                      stringsAsFactors=FALSE)

insect <- read.csv("data/insect-all.csv",
                   stringsAsFactors=FALSE)

source("../forestosmia/dataPrep/src/siteNames.R")
source("../forestosmia/dataPrep/src/misc.R")

dir <- "cleaneddata"

## ***********************************************************
## parasite data
## ***********************************************************

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

## dropped about 10 specimens

## We may want to compare parasites in unemerged and foraging bees, in
## which case use "parasitecontrol" object. However, for most
## analyses, we want to lose the unemerged, control bees, use
## "parasite"

parasitecontrol <- parasite
parasite <- parasitecontrol[!parasitecontrol$Stand == "control",]
dim(parasite)

## community health metrics
parasitenames <- c("Crithidia", "Ascophaera", "Apicystis")

## sum the 0,1 for each parasite
parasite$ParasiteRichness <- rowSums(parasite[, parasitenames],
                                     na.rm=TRUE)

## ranges from 0 to 3 (the max)
hist(parasite$ParasiteRichness)

## does an indivual have any 1+ parasites? trues/falses times 1 (true
## is a 1 and false is 0)
parasite$AnyParasite <- (parasite$ParasiteRichness > 0)*1
hist(parasite$AnyParasite)

## calculate avg sick individuals for each stand
## calculate crithidia, apicystis, and asospheara for each site
sick.totals <- parasite %>%
    group_by(Stand) %>%
    summarise(TestedTotals = length(UniqueID[!is.na(AnyParasite)]),
              ParasitismRate=mean(AnyParasite, na.rm=TRUE),
              InfectedIndividuals=sum(AnyParasite, na.rm=TRUE),
              InfectedCrith=sum(Crithidia, na.rm=TRUE),
              InfectedApicys=sum(Apicystis, na.rm=TRUE),
              InfectedAsco=sum(Ascophaera, na.rm=TRUE),
              CrithRate=mean(Crithidia, na.rm=TRUE),
              AscoRate=mean(Ascophaera, na.rm=TRUE),
              ApicysRate=mean(Apicystis, na.rm=TRUE))

## merge stand data and stand-level parasite data
standinfo <- merge(standinfo, sick.totals, all.x=TRUE)

f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(standinfo$ParasitismRate)
    hist(standinfo$CrithRate)
    hist(standinfo$AscoRate)
    hist(standinfo$ApicysRate)

}
pdf.f(f.hist, file="figures/parasitism.pdf", height=12, width=4)

## ***********************************************************
## floral data
## ***********************************************************

## This dataset includes floral data from first two sampling rounds,
## when we had osmia out\

## we want abundance flowers (abundflw), richness flowering species
## (richnessflwingsp), and abundfloweringsp

floral$Vials <- NULL
floral$Other_flowers <- NULL
floral$Notes <- NULL

## Floral Richness
## we don't want zeroes counted as 1 though

floral$Flower_sci[floral$Flower_sci == 0] <- NA
floral$Flower_sci <- fix.white.space(floral$Flower_sci)
id(floral$Flower_sci)

floral$Flower_sci <- as.character(floral$Flower_sci)

## clean up names
floral$Flower_sci[floral$Flower_sci == "Claytonia perfoliata"]  <-
    "Claytonia parviflora"
floral$Flower_sci[floral$Flower_sci == "Lupinus rivularis "]  <-
    "Lupinus rivularis"
floral$Flower_sci[floral$Flower_sci == "Unknown forb"]  <- NA
floral <- floral[!floral$Flower_sci == "NA",]
id(floral$Flower_sci)

checked.plant.names <- TPL(id(floral$Flower_sci))

floral$Family  <- NA
floral$Family <- checked.plant.names$Family[match(floral$Flower_sci,
                                                  checked.plant.names$Taxon)]

floral <- floral[!is.na(floral$Year)]
## all unIDed plants removed

floral.sum <- floral %>%
    group_by(Stand, Transect) %>%
    summarise(FlowerRichness =
                  length(unique(Flower_sci)),
              FlowerDiversity=
                  vegan:::diversity(table(Flower_sci),
                                    index="shannon"),
              BloomAbund=sum(Blooms, na.rm=TRUE),
              StemAbund=sum(Stems, na.rm=TRUE),
              CanopyCover=mean(Canopy_cent,  na.rm=TRUE))

mean.floral <- floral.sum %>%
    group_by(Stand) %>%
    summarise(FlowerRichness = mean(FlowerRichness, na.rm=TRUE),
              FlowerDiversity = mean(FlowerDiversity, na.rm=TRUE),
              MeanBloomAbund=mean(BloomAbund, na.rm=TRUE),
              MeanStemAbund=mean(StemAbund, na.rm=TRUE),
              MeanCanopy=mean(CanopyCover,  na.rm=TRUE))

standinfo <- merge(standinfo, mean.floral, all.x=TRUE)

standinfo$MeanBloomAbund[is.na(standinfo$MeanBloomAbund)] <- 0
standinfo$FlowerRichness[is.na(standinfo$FlowerRichness)] <- 0
standinfo$FlowerDiversity[is.na(standinfo$FlowerDiversity)] <- 0


f.hist <- function(){
    layout(matrix(1:3, ncol=1))
    hist(standinfo$MeanBloomAbund)
    hist(standinfo$FlowerRichness)
    hist(standinfo$FlowerDiversity)

}
pdf.f(f.hist, file="figures/floral.pdf",
      height=8, width=4)


## ***********************************************************
## BEE DIVERSITY DATA, MERGE
## ***********************************************************

insect$Caste <- NULL

bee.fams <- c("Halictidae", "Apidae", "Colletidae", "Megachilidae",
              "Andrenidae")
bee <- insect[insect$Family %in% bee.fams,]

id(bee$Genus)
## bee <- bee[!is.na(bee$Genus),]
## bee <- bee[bee$Genus != "0",]

bee$GenSp[bee$GenSp=="0 0"] <- NA
bee$GenSp[bee$GenSp=="NA NA"] <- NA


bee$GenSp <- fix.white.space(bee$GenSp)
id(bee$GenSp)

bee$GenSp[bee$GenSp  == "Andrena sp."]  <- NA
bee$GenSp[bee$GenSp  == "Ceratina damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Bombus damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Andrena damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Halictus sp."]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum dialictus"]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum dialictus sp."]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum dialictus damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum evylaeus damaged"]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum sp."]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum evylaeus"]  <- NA
bee$GenSp[bee$GenSp  == "Lasioglossum egregium"]  <- NA
bee$GenSp[bee$GenSp  == "Melissodes sp."]  <- NA

id(bee$GenSp)

## generate abundance of bees in each sample round
## calculate bee abundance with just net data
## calculate bee richness with net + pan data
bee.sum <- bee %>%
    group_by(Stand, Round, Trap.type, Year) %>%
    summarise(BeeRichness = length(unique(GenSp[!is.na(GenSp)])),
              BeeDiversity =
                  vegan::diversity(table(GenSp[!is.na(GenSp)]),
                                   index="shannon"),
              BeeAbund=sum(Count))

## mean richness across all collection types
mean.bee.rich <- bee.sum %>%
    group_by(Stand, Year) %>%
    summarise(MeanBeeRichness = mean(BeeRichness, na.rm=TRUE))

## mean for each collection type
mean.bee.abund <- bee.sum %>%
    group_by(Stand, Trap.type, Year) %>%
    summarise(MeanBeeAbund=mean(BeeAbund, na.rm=TRUE),
              MeanBeeDiversity=mean(BeeDiversity, na.rm=TRUE))

## mean.bee.abund <- mean.bee.abund[mean.bee.abund$Trap.type == "Net" &
##                                  mean.bee.abund$Year == "2019",]

## if we drop the non-netted data, most sites don't have data

mean.bee.abund <- bee.sum %>%
    group_by(Stand, Year) %>%
    summarise(MeanBeeAbund=mean(BeeAbund, na.rm=TRUE),
              MeanBeeDiversity=mean(BeeDiversity, na.rm=TRUE))

mean.bee.abund$MeanBeeAbund[is.na(mean.bee.abund$MeanBeeAbund)] <- 0
mean.bee.abund$MeanBeeDiversity[is.na(mean.bee.abund$MeanBeeDiversity)] <- 0

## We only want bees from 2019 when the parasite screenings took place
mean.bee.abund <- mean.bee.abund[mean.bee.abund$Year == "2019",]
mean.bee.rich <- mean.bee.rich[mean.bee.rich$Year == "2019",]

## remove columsn we don't want
mean.bee.abund$Trap.type <- NULL
mean.bee.abund$Year <- NULL
mean.bee.rich$Year <- NULL

                                        # merge to site-level data
standinfo <- merge(standinfo, mean.bee.abund, all.x=TRUE)
standinfo <- merge(standinfo, mean.bee.rich, all.x=TRUE)

standinfo$MeanBeeRichness[is.na(standinfo$MeanBeeRichness)] <- 0
standinfo$MeanBeeDiversity[is.na(standinfo$MeanBeeDiversity)] <- 0
standinfo$MeanBeeAbund[is.na(standinfo$MeanBeeAbund)] <- 0

## ***********************************************************
## honey bees
## ***********************************************************

## honey bees from net data

hb.bee <- bee[bee$Genus == "Apis",]

hb.sum <- hb.bee %>%
    group_by(Stand, Round, Trap.type, Year) %>%
    summarise(BeeAbund=sum(Count, na.rm=TRUE))

mean.hb <- hb.sum %>%
    group_by(Stand, Trap.type, Year) %>%
    summarise(MeanHBAbund=mean(BeeAbund, na.rm=TRUE))


## We only want hb that were netting and bees from 2019
mean.hb <- mean.hb[mean.hb$Trap.type == "Net" &
                   mean.hb$Year == "2019",]

mean.hb$Trap.type <- NULL
mean.hb$Year <- NULL

standinfo <- merge(standinfo, mean.hb, all.x=TRUE)


f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(standinfo$MeanBeeRichness)
    hist(standinfo$MeanBeeDiversity)
    hist(standinfo$MeanBeeAbund)
    hist(standinfo$MeanHBAbund)

}
pdf.f(f.hist, file="figures/bees.pdf",
      height=12, width=4)


## ***********************************************************
## nest block data
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


## create a new block level dataset to describe reprodata, with stand
## info merged
repro.block <- merge(repro.block, standinfo)

## merge site level sick totals to our reproductive, block-level dataset
repro.block <- merge(repro.block, sick.totals)

## INDIVIDUAL LEVEL DATASET ********

## the parasite dataset is our individual level data, let's merge tube
## level data to this dataset and rename it as our indiv dataset

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

## merge stand info onto parasite (individual data) and rename it as
## indivi. data
dim(parasite)
indiv.data <- merge(parasite, standinfo, all=TRUE)
dim(indiv.data)

## *************************************************************
## write out final datasets: stand-level, block-level, indiv-level
## *************************************************************
## we are done with the stand-level dataset, let's rename it
site.data <- standinfo

write.csv(indiv.data, file=file.path(dir,
                                     "indivData.csv"), row.names=FALSE)

write.csv(standinfo, file=file.path(dir,
                                    "standinfo.csv"), row.names=FALSE)

write.csv(repro.block, file=file.path(dir,
                                      "NestRepro.csv"), row.names=FALSE)


save.dir <- "../forestosmia/data"
save(indiv.data,
     file=file.path(save.dir, "indivdata.Rdata"))


save(repro.block,
     file=file.path(save.dir, "reproblock.Rdata"))

save(site.data,
     file=file.path(save.dir, "sitedata.Rdata"))



## summary stats for ms

only.screened <- indiv.data[!is.na(indiv.data$TestedTotals),]

sites.screened <- unique(only.screened$Stand)

bee.screened  <- bee[bee$Stand %in% sites.screened,]
bee.screened <- bee.screened[bee.screened$Year == "2019",]

bee.sum.screened <- bee.screened %>%
    summarise(BeeRichness = length(unique(GenSp[!is.na(GenSp)])),
              BeeDiversity =
                  vegan::diversity(table(GenSp[!is.na(GenSp)]),
                                   index="shannon"),
              BeeAbund=sum(Count),
              BeeGenera=length(unique(Genus[!is.na(Genus)])),)



floral.screened  <- floral[floral$Stand %in% sites.screened,]
floral.screened <- floral.screened[floral.screened$Year == "2019",]

floral.sum.screened <- floral.screened%>%
    summarise(FlowerRichness =
                  length(unique(Flower_sci)),
              FlowerDiversity=
                  vegan:::diversity(table(Flower_sci),
                                    index="shannon"),
              BloomAbund=sum(Blooms, na.rm=TRUE))
