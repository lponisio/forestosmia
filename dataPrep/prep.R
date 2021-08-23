rm(list=ls())
library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)
## prepares raw data and creates dataset for analyses

setwd('/Volumes/bombus/Dropbox (University of Oregon)/')

setwd("forestosmia_saved")
dir <- "forestosmia_saved"

parasite <- read.csv("data/parasite.csv",
                     stringsAsFactors=FALSE)

dbh <- read.csv("data/dbh.csv",
                stringsAsFactors=FALSE, na.strings=c("","NA"))

floral <- read.csv("data/floral2.csv", stringsAsFactors=FALSE)


repro <- read.csv("data/reproductive.csv",
                  stringsAsFactors=FALSE)

standinfo <- read.csv("data/standinfo.csv",
                      stringsAsFactors=FALSE)

vegcover <- read.csv("data/vegcover.csv",
                     stringsAsFactors=FALSE)

insect <- read.csv("data/insect-all.csv",
                   stringsAsFactors=FALSE)

foraging <- read.csv("data/osmiaforaging.csv",
                     stringsAsFactors=FALSE)

source("../forestosmia/dataPrep/src/siteNames.R")
source("../forestosmia/dataPrep/src/misc.R")

## **********************************************************
## stand level calculations
## **********************************************************

standinfo$AgePoly1 <-  poly(standinfo$Age,3)[,1]
standinfo$AgePoly2 <-  poly(standinfo$Age,3)[,2]
standinfo$AgePoly3 <-  poly(standinfo$Age,3)[,3]
standinfo$AgePoly4 <-  poly(standinfo$Age,4)[,4]
standinfo$AgePoly5 <-  poly(standinfo$Age,5)[,5]

## ***********************************************************
## veg and stand data
## ***********************************************************

## calculate stand-level broadleaf cover (stand=site).  keep only data
## for "broadleaf" trees, srubs, and forbs.  when dropping a variable,
## check the class with class(vegcover$broadleaf) and if its an
## integar dont put quotes around 0.

## not sure why this was there before...
## vegcover <- vegcover[!vegcover$Broadleaf == 0,]

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
standinfo$BLcover <- c(BLcover.stand[match(standinfo$Stand,
                                           names(BLcover.stand))])

standinfo$BLcover[is.na(standinfo$BLcover)] <- 0

## *************************************************************
## create stand intensity
## *************************************************************

hist(standinfo$BLcover, breaks=30,
     xlim = c(0, max(standinfo$BLcover, na.rm=TRUE)))

standinfo$StandIntensity <- "3high"
standinfo$StandIntensity[standinfo$BLcover > 5 &
                         standinfo$BLcover <= 12] <- "2medium"

standinfo$StandIntensity[standinfo$BLcover > 12] <- "1low"

standinfo$StandIntensity[is.na(standinfo$BLcover)] <- NA

table(standinfo$StandIntensity, standinfo$Age)

## ***********************************************************
##  dbh
## ***********************************************************

## let's remove the columns that are not correctly filled out
dbh$LiveorDead <- NULL
dbh$SlopePerc <- NULL
dbh$Aspect <- NULL

## Right now blanks are assigned as NA, but a blank TreeNum means no
## trees, which we actually want to account for.  If there's an NA
## under TreeNum and DBH, make it zeroes

dbh$TreeNum[is.na(dbh$TreeNum)] <- 999 # place holder
dbh$DBH[is.na(dbh$DBH)] <- 0

## we want tree species richness, tree species abundance, tree DBH
## averaged across trees at a site

dbh$DBH <- as.numeric(dbh$DBH)
dbh$TreeNum <- as.numeric(dbh$TreeNum)

stand.sum <- dbh %>%
    group_by(Stand) %>%
    summarise(TreeRichness = length(unique(Species)),
              TreeDiversity=vegan:::diversity(table(Species),
                                              index="shannon"),
              TreeAbund=length(TreeNum),
              MeanDBH=mean(DBH, na.rm=TRUE))

## then merge to dataset
standinfo <- merge(standinfo, stand.sum)



f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(standinfo$TreeRichness)
    hist(standinfo$TreeDiversity)
    hist(standinfo$MeanDBH)
    hist(standinfo$BLcover)

}
pdf.f(f.hist, file="figures/standInfo.pdf", height=12, width=4)

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

## let's drop the column "vials" since it's not correctly filled out
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

## drop round 3 since it is after the osmia were put out?
## bee <- bee[!bee$Round == 3,]

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

## drop count of 0 which are placeholders for sites where no bees
## where caught. Add back in later
## bee <- bee[bee$Count == 1,]

## drop unIDed bees
## bee <- bee[!is.na(bee$GenSp),]

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
indiv.data <- merge(parasite, standinfo, all.x=TRUE)
dim(indiv.data)

## ***********************************************************
## CLEAN FORAGING DATA, MERGE
## ***********************************************************

## we are just looking at foraging distance of females from Osmia
## nests in blocks where samples were not pulled for parasite tests

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

## drop any data where trip length could not be calculated due to lack
## of footage

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

## let's create site-level data and merge onto site and individual
## level datasets

foraging.stand <- foraging %>%
    group_by(Stand) %>%
    summarise(MeanTripLengthStand=mean(TripLength.h.mm.ss, na.rm = TRUE))

## merge stand-level
standinfo <- merge(standinfo, foraging.stand, all.x=TRUE)

## merge individual data
dim(indiv.data)
indiv.data <- merge(indiv.data, foraging.stand, all.x=TRUE)
dim(indiv.data)


## *************************************************************
## write out final datasets: stand-level, block-level, indiv-level
## *************************************************************
## we are done with the stand-level dataset, let's rename it
site.data <- standinfo


dir <- "cleaneddata"
parasite <- merge(parasite, standinfo, all.x=TRUE)

write.csv(indiv.data, file=file.path(dir,
                                     "indivData.csv"), row.names=FALSE)



## write.csv(parasite, file=file.path(dir,
##                                    "parasite.csv"), row.names=FALSE)


write.csv(standinfo, file=file.path(dir,
                                    "standinfo.csv"), row.names=FALSE)

write.csv(repro.block, file=file.path(dir,
                                      "NestRepro.csv"), row.names=FALSE)


save.dir <- "../forestosmia/data"
save(indiv.data,
     file=file.path(save.dir, "indivdata.Rdata"))


## save(parasite,
##      file=file.path(save.dir, "parasite.Rdata"))

save(repro.block,
     file=file.path(save.dir, "reproblock.Rdata"))

save(site.data,
     file=file.path(save.dir, "sitedata.Rdata"))


## **********************************************************
## data exploration
## **********************************************************

p <- ggplot(data=site.data, aes(x=Stand, y=ParasitismRate, fill=MeanBeeAbund)) +
    geom_bar(stat="identity") +
    labs(y="Proportion infected") +
     scale_color_viridis(option = "D")

# Horizontal bar plot
p + coord_flip()

ggsave(file="figures/bysite_beeAbund.pdf")

p <- ggplot(data=site.data, aes(x=Stand, y=ParasitismRate, fill=MeanBeeAbund)) +
    geom_bar(stat="identity") +
    labs(y="Proportion infected") +
     scale_color_viridis(option = "D")

# Horizontal bar plot
p + coord_flip()
ggsave(file="figures/bysite_beeAbund.pdf")


f2 <- function(){
    layout(matrix(1:5, ncol=1))
    par(oma=c(2,2,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))

    plot(log(MeanBloomAbund + 1) ~Age, data=standinfo, pch=16,
         ylab ="",
         xlab="", cex=1.2, cex.lab=1.2)
    mtext("Mean Flower Abundance (log)", 2, line=2)


    plot(FlowerDiversity~Age, data=standinfo, pch=16,
         ylab ="",
         xlab="", cex=1.2, cex.lab=1.2)
    mtext("Mean Flower Diversity", 2, line=2)

    plot(MeanBeeAbund ~Age, data=standinfo, pch=16,
         ylab ="",
         xlab="", cex=1.2, cex.lab=1.2)
    mtext("Mean Bee Abundance", 2, line=2)


    plot(MeanBeeRichness~Age, data=standinfo, pch=16,
         ylab ="",
         xlab="", cex=1.2, cex.lab=1.2)
    mtext("Mean Bee Richness", 2, line=2)

    plot(MeanBeeDiversity~Age, data=standinfo, pch=16,
         ylab ="",
         xlab="", cex=1.2, cex.lab=1.2)

    mtext("Mean Bee Diversity", 2, line=2)
    mtext("Years post harvest", 1, line=2)
}


pdf.f(f2, file="figures/yph.pdf",
      height=12, width=4)


## *************************************************************
##  map
## *************************************************************

## library(sp)
## library(rgdal)

## standinfo <- standinfo[standinfo$Age <= 9,]
## table(standinfo$StandIntensity, standinfo$Age)

## pts.standinfo <- standinfo
## coordinates(pts.standinfo)=~CenterLon+CenterLat
## proj4string(pts.standinfo)=CRS("+init=epsg:4326") # set it to lat-long

## sites <- SpatialPointsDataFrame(coordinates(pts.standinfo),
##                                    data=standinfo)

## rgdal::writeOGR(sites, dsn="data/spatial", "sites",
## driver="ESRI Shapefile")

## *************************************************************
##  plant lists
## *************************************************************
## floral <- merge(floral, standinfo)

## floral$AgeCat <- "yph0-10"
## ## floral$AgeCat[floral$Age > 3] <- "yph4-10"
## floral$AgeCat[floral$Age > 10] <-  "yph11+"


## ## by age
## flowers.age <- table(floral$Flower_sci, floral$Age)

## flowers.mean <-  floral %>%
##   group_by(Flower_sci, Age) %>%
##   summarise(FloralMeanAbund=mean(Blooms, na.rm=TRUE))

## write.csv(flowers.age, file="data/floral.csv")
## write.csv(flowers.mean, file="data/floralMeanAge.csv")

## flowers <- read.csv(file="data/floral.csv")
## colnames(flowers) <- c("PlantGenusSpecies", colnames(flowers)[-1])

## flowers$sumSites <- apply(flowers[,-1], 1, sum)
## write.csv(flowers, file="data/floral.csv", row.names=FALSE)

## ## by age group

## flowers.age.cat <-  floral %>%
##   group_by(Flower_sci, AgeCat) %>%
##   summarise(FloralMeanAbund=mean(Blooms, na.rm=TRUE))

## write.csv(flowers.age.cat, file="data/floralMeanCatAge.csv",
##           row.names=FALSE)

## flowers.age.cat <- flowers.age.cat[flowers.age.cat$AgeCat == "yph0-10",]

## ## *****

## mixes <- read.csv("../beebettertimber/seeds/mixes.csv")
## mixes$GenusSpecies <- fix.white.space(mixes$GenusSpecies)

## id(mixes$GenusSpecies)

## mixes$Genus <- sapply(strsplit(mixes$GenusSpecies, " "),
##                       function(x)x[1])

## flowers.age.cat$Genus <- sapply(strsplit(flowers.age.cat$Flower_sci, " "),
##                       function(x)x[1])


## matchMixes <- function(flower.data, flower.name.col, flower.mix.col,
##                        mixes){
##     flower.data <- data.frame(flower.data)
##     tt <- mixes[, flower.mix.col][mixes$Mix == "Tough and tenacious"]
##     dp <- mixes[, flower.mix.col][mixes$Mix == "Diverse prarie"]
##     us <- mixes[, flower.mix.col][mixes$Mix == "Understory"]
##     bp <- mixes[, flower.mix.col][mixes$Mix == "Burn pile"]
##     dg <- mixes[, flower.mix.col][mixes$Mix == "Disturbed ground"]
##     cm <- mixes[, flower.mix.col][mixes$Mix == "Custom"]
##     cc <- mixes[, flower.mix.col][mixes$Mix == "Clearcut_herbicide"]
##     ccn <- mixes[, flower.mix.col][mixes$Mix == "Clearcut_natural"]
##     es <- mixes[, flower.mix.col][mixes$Mix == "Seral Ref"]
##     usfs <- mixes[, flower.mix.col][mixes$Mix == "USFS"]

##     flower.data$Tenacious <- 0
##     flower.data$Tenacious[flower.data[, flower.name.col] %in% tt]  <- 1

##     flower.data$DiversePrarie <- 0
##     flower.data$DiversePrarie[flower.data[, flower.name.col] %in% dp]  <- 1

##     flower.data$Understory <- 0
##     flower.data$Understory[flower.data[, flower.name.col] %in% us]  <- 1

##     flower.data$BurnPile <- 0
##     flower.data$BurnPile[flower.data[, flower.name.col] %in% bp]  <- 1

##     flower.data$Disturbed <- 0
##     flower.data$Disturbed[flower.data[, flower.name.col] %in% dg] <- 1

##     flower.data$Custom <- 0
##     flower.data$Custom[flower.data[, flower.name.col] %in% cm] <- 1

##     flower.data$EarlySeralRef <- 0
##     flower.data$EarlySeralRef[flower.data[, flower.name.col] %in% es] <- 1

##     flower.data$ClearcutHerbicide <- 0
##     flower.data$ClearcutHerbicide[flower.data[, flower.name.col] %in% cc] <- 1

##     flower.data$Clearcut21year <- 0
##     flower.data$Clearcut21year[flower.data[, flower.name.col] %in%
##                                ccn] <- 1

##     flower.data$USFS <- 0
##     flower.data$USFS[flower.data[, flower.name.col] %in% usfs] <- 1


##     return(flower.data)
## }


## flowers.age.cat.sp <- matchMixes(flowers.age.cat, "Flower_sci",
##                               "GenusSpecies", mixes)

## write.csv(flowers.age.cat.sp, file="data/ClearCutFlowerMean.csv",
##           row.names=FALSE)

## ## by genus

## flowers.age.cat.genus <- matchMixes(flowers.age.cat, "Genus", "Genus",
##                                     mixes)

## flowers.age.cat.genus$Flower_sci  <- NULL
## flowers.age.cat.genus$FloralMeanAbund <- NULL

## flowers.age.cat.genus <-
##     flowers.age.cat.genus[!duplicated(flowers.age.cat.genus),]


## write.csv(flowers.age.cat.genus, file="data/GenusClearCutFlowerMean.csv",
##           row.names=FALSE)

## ## *****

## seral <- read.csv("../beebettertimber/seeds/earlySeral.csv")

## seral$Group <- fix.white.space(seral$Group)
## seral$GenusSpecies <- fix.white.space(seral$GenusSpecies)
## seral$Family <- fix.white.space(seral$Family)

## cats <- c("early seral herb",
##           "late seral herb",
##           "forest generalist shrub",
##           "forest generalist herb",
##           "unclassified herb",
##           "unclassified shrub")

## seral <- seral[seral$Group %in% cats,]
## seral <- seral[seral$Family != "Poaceae",]

## seral.mixes <- matchMixes(seral, "GenusSpecies",
##                           "GenusSpecies", mixes)

## write.csv(seral.mixes, file="data/SeralRefMixes.csv",
##           row.names=FALSE)



