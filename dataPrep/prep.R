rm(list=ls())
library(dplyr)
library(lubridate)
library(Taxonstand)
## prepares raw data and creates dataset for analyses
setwd('~/Dropbox (University of Oregon)/')

setwd("forestosmia_saved")

parasite <- read.csv("data/parasite.csv",
                     stringsAsFactors=FALSE)
floral2 <- read.csv("data/final_Floral_2018_2019_01192023.csv",
                   stringsAsFactors=FALSE)
repro <- read.csv("data/final_repro_cleaned_SG_20230110.csv",
                  stringsAsFactors=FALSE)
stand.info <- read.csv("data/standinfo_Nov9_2022.csv",
                       stringsAsFactors=FALSE)
insect2 <- read.csv("data/final_AllSpecimens_2018_2019_12312022.csv",
                   stringsAsFactors=FALSE)
source("../forestosmia/dataPrep/src/siteNames.R")
source("../forestosmia/dataPrep/src/misc.R")

dir <- "cleaneddata"

## ## ***********************************************************
## ## prep for merging stand data with multiple years
## ## ***********************************************************
floral.stands <- unique(floral2$Stand[floral2$Year == 2019])
insect.stands <- unique(insect2$Stand[insect2$Year == 2019])

## check that all stands in the 2019 insect data are in the floral,
## and vice versa
insect.stands[!insect.stands %in% floral.stands]
floral.stands[!floral.stands %in% insect.stands]

## check to see that all stands in stand info are in the insect data
stand.info$Stand[!stand.info$Stand %in% insect.stands]
stand.info$Stand[stand.info$Stand == "Alexander Rd"]  <-
    "Alexander Rd."
repro$Stand[repro$Stand == "Alexander Rd"]  <- "Alexander Rd."

stand.info$Stand[!stand.info$Stand %in% insect.stands]
## Honeygrove was sampled in 2018 but access was blocked in 2019
stand.info$Year <- 2018
stand.info$Year[stand.info$Stand %in% floral.stands]  <- 2019

sum(stand.info$Year == 2018)

## drop the 2018 only stands, which is only one stand Honeygrove
stand.info <- stand.info[stand.info$Year == 2019,]

stand.info$Stand[!stand.info$Stand %in% unique(repro$Stand)]

## ***********************************************************
## parasite data
## ***********************************************************
dim(parasite)

## drop the NOTES column
parasite$NOTES <- NULL

parasite$Date <- as.Date(parasite$Date, "%m/%d/%y")
parasite$Year <- format(parasite$Date, "%Y")

## drop all the no template controls under UniqueID because they were
## all neg for parasites
parasite <- parasite[!parasite$UniqueID == "C",]
dim(parasite)

## specimens where bee DNA did not extract properly
sum(parasite$ApidaeCtrl == 0)

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

## calculate avg sick individuals for each stand
## calculate crithidia, apicystis, and asospheara for each site
sick.totals <- parasite %>%
    group_by(Stand,Year) %>%
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
stand.info <- merge(stand.info, sick.totals, all.x=TRUE)

f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(stand.info$ParasitismRate, main="Any parasite", xlab="")
    hist(stand.info$CrithRate, main= "Crithidia prevalence", xlab="")
    hist(stand.info$AscoRate, main= "Ascophaera prevalence", xlab="")
    hist(stand.info$ApicysRate, main="Apicystis prevalence", xlab="")

}
pdf.f(f.hist, file="figures/parasitism.pdf", height=12, width=4)

## ***********************************************************
## floral data
## ***********************************************************

## This dataset includes floral data from first two sampling rounds,
## when we had osmia out

## we want abundance flowers (abundflw), richness flowering species
## (richnessflwingsp), and abundfloweringsp

cleanFloral <- function(floral){
    floral$Other_flowers <- NULL
    floral$Notes <- NULL
    floral$Flower_common <- NULL

    ## Floral Richness
    ## we don't want zeroes counted as 1 though

    floral$Flower_sci[floral$Flower_sci == 0] <- "No flowers"
    floral$Flower_sci <- fix.white.space(floral$Flower_sci)
    floral$Flower_sci <- as.character(floral$Flower_sci)

    ## clean up names
    floral$Flower_sci[floral$Flower_sci == "Claytonia perfoliata"]  <-
        "Claytonia parviflora"
    floral$Flower_sci[floral$Flower_sci == "Lupinus rivularis "]  <-
        "Lupinus rivularis"
    floral$Flower_sci[floral$Flower_sci == "Unknown forb"]  <- NA
    floral$Flower_sci[floral$Flower_sci == "Iris sp."]  <-
        "Iris tenax"
    floral$Flower_sci[floral$Flower_sci == "Corylus sp."]  <- "Corylus cornuta"
    print(id(floral$Flower_sci))

    print(paste("number of duplicated rows to drop", sum(duplicated(floral))))
    ## remove duplicates
    dup.floral <- floral[duplicated(floral),]
    write.csv(dup.floral, file="data_cleaning/duplicated_veg_rows.csv")

    floral <- floral[!duplicated(floral),]

    ## checked.plant.names <- TPL(id(floral$Flower_sci))
    ## write.csv(checked.plant.names, file="checkedPlants.csv")

    ## floral$Family  <- NA
    ## floral$Family <- checked.plant.names$Family[match(floral$Flower_sci,
    ##                                                   checked.plant.names$Taxon)]

    floral$Date <- as.Date(floral$Date, "%m/%d/%y")
    floral$Year <- format(floral$Date, "%Y")

    floral <- floral[!is.na(floral$Stand),]
    return(floral)
}


floral <- cleanFloral(floral2)
nrow(floral)

## only use 2019 data
floral <- floral[floral$Year == 2019,]
nrow(floral)

split.floral <- split(floral, floral$Stand)

## rounds are all over the please, fill in based on date
split.floral <- lapply(split.floral, function(x){
    x <- x[order(x$Date),]
    rounds.count <- c(table(x$Date))
    nround <- length(names(rounds.count))
    rounds <- 1:nround
    rs <- lapply(rounds, function(r)  rep(r, rounds.count[r]))
    rs <- do.call(c, rs)
    x$Round <- rs
    return(x)
})

floral <- do.call(rbind, split.floral)

## drop round 4 which does not overlap when Osmia where flying
floral <- floral[floral$Round != 4,]
nrow(floral)

floral.check <- floral %>%
    group_by(Stand,  Year, Round) %>%
    summarise(length(unique(Transect)))

write.csv(floral.check, file = "floral_samples.csv",
          row.names=FALSE)

## all unIDed plants removed

## calculate richness and diversity for each transect, dropping the
## "no flowers" rows
floral.sum <- floral %>%
    group_by(Stand, Transect, Year, Round) %>%
    summarise(FlowerRichness =
                  length(unique(Flower_sci[Flower_sci != "No flowers"])),
              FlowerDiversity=
                  vegan:::diversity(table(Flower_sci[Flower_sci !=
                                                     "No flowers"]),
                                    index="shannon"),
              BloomAbund=sum(Blooms, na.rm=TRUE),
              StemAbund=sum(Stems, na.rm=TRUE),
              CanopyCover=mean(Canopy_cent,  na.rm=TRUE))

## take the mean across the year, since that would be have the average
## conditions the osmia were exposed to


mean.floral <- floral.sum %>%
    group_by(Stand, Year) %>%
    summarise(FlowerRichness = mean(FlowerRichness, na.rm=TRUE),
              FlowerDiversity = mean(FlowerDiversity, na.rm=TRUE),
              MeanBloomAbund=mean(BloomAbund, na.rm=TRUE),
              MeanStemAbund=mean(StemAbund, na.rm=TRUE),
              MeanCanopy=mean(CanopyCover,  na.rm=TRUE))

stand.info <- merge(stand.info, mean.floral, all.x=TRUE)

only.na <- stand.info$Stand[is.na(stand.info$MeanBloomAbund)]

## check stands with with NA for site summaries
floral[floral$Stand %in% only.na,]

## all should be zero because no flowers ever found. No longer
## necessary as of Dec 8 2022 because all stands have values

f.hist <- function(){
    layout(matrix(1:3, ncol=1))
    hist(stand.info$MeanBloomAbund, main= "Flower abundance", xlab="")
    hist(stand.info$FlowerRichness, main = "Flower Richness", xlab="")
    hist(stand.info$FlowerDiversity, main= "Flower diversity", xlab="")

}
pdf.f(f.hist, file="figures/floral_new_data.pdf",
      height=8, width=4)

checkMax <- function(floral){

    ## investigate stand with ~max blooms
    etreme.stand <- stand.info[stand.info$MeanBloomAbund ==
                           max(stand.info$MeanBloomAbund),]
    ## given by lots and lots of foxglove
    print("Stand with max blooms")
    floral[floral$Stand == unique(etreme.stand$Stand),]
}

checkMax(floral)

## ***********************************************************
## BEE DIVERSITY DATA, MERGE
## ***********************************************************

cleanBees <- function(insect){
    insect$Caste <- NULL
    insect$Notes <- NULL
    bee.fams <- c("Halictidae", "Apidae", "Colletidae", "Megachilidae",
                  "Andrenidae")
    bee <- insect[insect$Family %in% bee.fams,]
    print(id(bee$Genus))

    bee$GenusSpecies <- fix.white.space(paste(bee$Genus, bee$Species))
    bee$GenusSpecies[bee$GenusSpecies=="NA NA"] <- NA
    bee$GenusSpecies[grep("damaged", bee$GenusSpecies)] <- NA

    ## print(id(bee$GenusSpecies))
    ## maybe leave Andrena sp. if it is only at one place?
    bee$GenusSpecies[bee$GenusSpecies  == "Andrena sp."]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Halictus sp."]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Lasioglossum dialictus"]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Lasioglossum dialictus sp."]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Lasioglossum sp."]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Lasioglossum evylaeus"]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Lasioglossum egregium"]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == "Melissodes sp."]  <- NA
    bee$GenusSpecies[bee$GenusSpecies  == ""]  <- NA

    ## lump together based on Koch et al. 2018 molecular evidence and
    ## Linc Best confirmation
    bee$GenusSpecies[bee$GenusSpecies  == "Bombus californicus"]  <- "Bombus fervidus"

    print(id(bee$GenusSpecies))

    bee$Date <- as.Date(bee$Date, "%Y-%m-%d")
    bee$Year <- format(bee$Date, "%Y")

    ## print(paste("rows before dropping NA species ids", nrow(bee)))
    ## bee <- bee[!is.na(bee$GenusSpecies),]

    ## print(paste("rows after dropping NA species ids", nrow(bee)))

    print(paste("number of duplicated rows to be dropped", sum(duplicated(bee))))

    bee <- bee[!duplicated(bee),]
    return(bee)
}

bee <- cleanBees(insect2)

## drop round 4
bee <- bee[bee$Round != 4,]

## drop years before 2019
bee <- bee[bee$Year == 2019,]

## create a method column to distingish pans and van
bee$Method <- "Net"
bee$Method[bee$Trap.type %in% c("BPT", "WPT", "YPT")]  <- "Pan"
bee$Method[bee$Trap.type %in% c("BVT")]  <- "Vane"

print("number of rows with duplicated uniqueIDs with the same species ID")
print(sum(duplicated(cbind(bee$ID, bee$GenusSpecies))))

print("number of rows with duplicated uniqueIDs")
print(sum(duplicated(bee$ID)))

bee.check <- bee %>%
    group_by(Stand,  Year, Round) %>%
    summarise(methods=length(unique(Trap.type)))

write.csv(bee.check, file = "samples.csv", row.names=FALSE)

sample.type <- unique(bee$Trap.type)

## create a master sheet of all the possible vane, pan and net samples
combos <- expand.grid(Stand= unique(bee$Stand),
                      Round= 1:3,
                      Trap.type  =sample.type)

combos$Year  <- 2019
write.csv(combos, file = "master_samples.csv", row.names=FALSE)


## generate abundance of bees in each sample round, combine each method (pan, net, vane)
bee.sum <- bee %>%
    group_by(Stand, Round,  Year, Trap.type) %>%
    summarise(BeeRichness = length(unique(GenusSpecies[!is.na(GenusSpecies)])),
              BeeDiversity =
                  vegan::diversity(table(GenusSpecies[!is.na(GenusSpecies)]),
                                   index="shannon"),
              BeeAbund=sum(Count), ## 1 if its a real bee
              HBAbund = sum(GenusSpecies == "Apis mellifera", na.rm=TRUE))

bee.sum.combo <- merge(combos, bee.sum, all.x=TRUE)

write.csv(bee.sum.combo, file = "all_combos.csv", row.names=FALSE)

## need to fix
bee.sum.combo[is.na(bee.sum.combo$BeeRichness), c("BeeRichness",
                                                  "BeeDiversity",
                                                  "BeeAbund",
                                                  "HBAbund")]  <- 0

## take the mean across sample types across the year at a stand
mean.bee <- bee.sum.combo %>%
    group_by(Stand, Year) %>%
    summarise(MeanBeeAbund=mean(BeeAbund, na.rm=TRUE),
              MeanBeeDiversity=mean(BeeDiversity, na.rm=TRUE),
              MeanBeeRichness=mean(BeeRichness, na.rm=TRUE),
              MeanHBAbund=mean(HBAbund, na.rm=TRUE))

## merge to site-level data
stand.info <- merge(stand.info, mean.bee, all.x=TRUE)

## stands with no bees netted or in pans
no.bees <- unique(stand.info$Stand[is.na(stand.info$MeanBeeAbund)])
bee$Stand[bee$Stand %in% no.bees]

stand.info$MeanBeeRichness[is.na(stand.info$MeanBeeRichness)] <- 0
stand.info$MeanBeeDiversity[is.na(stand.info$MeanBeeDiversity)] <- 0
stand.info$MeanBeeAbund[is.na(stand.info$MeanBeeAbund)] <- 0
stand.info$MeanHBAbund[is.na(stand.info$MeanHBAbund)] <- 0

## investigate stand with ~40 bees
etreme.stand.bee <- stand.info[stand.info$MeanBeeAbund ==
                           max(stand.info$MeanBeeAbund),]

print("Stand with max bees")
unique(etreme.stand.bee$Stand)
print("number of bees")
print(max(stand.info$MeanBeeAbund))

## bee[bee$Stand == etreme.stand.bee$Stand,]

f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(stand.info$MeanBeeRichness, main="Bee richness", xlab="")
    hist(stand.info$MeanBeeDiversity, main= "Bee diversity", xlab="")
    hist(stand.info$MeanBeeAbund, main= "Bee abundance", xlab="")
    hist(stand.info$MeanHBAbund, main= "Honey bee abundance", xlab="")

}
pdf.f(f.hist, file="figures/bees_newdata.pdf",
      height=12, width=4)


## ***********************************************************
## nest block data
## ***********************************************************
## we are just using females from Osmia nests in blocks
## where samples were not pulled for parasite tests

## TUBE LEVEL DATA ********

## drop any data from the Pathogen study
repro <- repro[repro$box_type != "Path",]

## confirm that only data is from reproductive study, "FPH" blocks
unique(repro$box_type)

## ## rename "nest number" as column
## colnames(repro)[colnames(repro) == "NestNum_XRAY"] <- "Column"

## to get average number of Females/Males at a stand, by nest number,
## we want to give each tube a unique ID, called "BlockNestTube"
## repro$BlockNestTube  <- paste0(repro$Block, repro$Row,
## repro$Column)

parasite$BlockNestTube  <- paste0(parasite$Block,
                                  parasite$Nest)

## repro$SumOffspring <- repro$Females + repro$Males

## BLOCK LEVEL DATA ********

## for each stand, average the number of females and females at each block
## end up with Block A averages and Block B averages for each site

repro <- repro[repro$Year == "2019",]

## dropping due to an incomplete xray
repro <- repro[! (repro$Stand == "Moonshine" & repro$Block == "A"),]

repro.block <- aggregate(list(Females=repro$Females,
                              Males=repro$Males),
                         list(Stand=repro$Stand,
                              Block=repro$Block,
                              Year=repro$Year),
                         sum, na.rm=TRUE)

## check over data
repro.check <- aggregate(repro$Block,
                         list(Stand=repro$Stand,
                              Block=repro$Block,
                              Row=repro$Row,
                              Year=repro$Year),
                         length)

## all 8 as expected!!! Jan 12-2023
write.csv(repro.check, file = "block_row_counts.csv", row.names=FALSE)

repro.check2 <- aggregate(repro$Block,
                         list(Stand=repro$Stand,
                              Block=repro$Block,
                              Year=repro$Year),
                         length)

## all 32 as expected!! Jan 12-2023
write.csv(repro.check2, file = "block_counts.csv", row.names=FALSE)


## add more offspring summary data to block-level averages
## above 1 F > M, below 1 M>F, = 1 M=F
repro.block$FM_ratio <-  log(repro.block$Females +1)/log(repro.block$Males +1)

repro.block$SumOffspring <-  repro.block$Females + repro.block$Males

## create a new block level dataset to describe reprodata, with stand
## info merged
repro.block <- merge(repro.block, stand.info)

f.hist <- function(){
    layout(matrix(1:4, ncol=1))
    hist(repro.block$SumOffspring, main="Any parasite", xlab="")
    hist(repro.block$Females, main= "Crithidia prevalence", xlab="")
    hist(repro.block$AscoRate, main= "Ascophaera prevalence", xlab="")
    hist(repro.block$ApicysRate, main="Apicystis prevalence", xlab="")

}
pdf.f(f.hist, file="figures/offspring.pdf", height=12, width=4)


## merge site level sick totals to our reproductive, block-level dataset
repro.block <- merge(repro.block, sick.totals, all.x=TRUE)

## merge stand info onto parasite (individual data) and rename it as
## indivi. data
dim(parasite)
indiv.data <- merge(parasite, stand.info, all=TRUE)
dim(indiv.data)

## *************************************************************
## write out final datasets: stand-level, block-level, indiv-level
## *************************************************************
## we are done with the stand-level dataset, let's rename it
site.data <- stand.info

drop.cols <- c("CenterLat", "CenterLon", "MeanStemAbund", "Age_LandOwner")

indiv.data <- indiv.data[!colnames(indiv.data) %in% drop.cols]
stand.info <- stand.info[!colnames(stand.info) %in% drop.cols]

repro.block <- repro.block[, c("Stand", "ParasitismRate",
                               "SumOffspring",
                               "Owner", "Hectares",
                               "Age_LandTrendr", "Elev",
                               "CenterLat", "CenterLon",
                               "FlowerDiversity",
                               "MeanBloomAbund",
                               "MeanBeeAbund", "MeanBeeDiversity",
                               "MeanHBAbund", "Block", "Year")]

write.csv(indiv.data, file=file.path(dir,
                                     "indivData.csv"), row.names=FALSE)
write.csv(stand.info, file=file.path(dir,
                                    "standInfo.csv"), row.names=FALSE)
write.csv(repro.block, file=file.path(dir,
                                      "NestRepro.csv"), row.names=FALSE)


save.dir <- "../forestosmia/data"
save(indiv.data,
     file=file.path(save.dir, "indivdata.Rdata"))


save(repro.block,
     file=file.path(save.dir, "reproblock.Rdata"))

save(site.data,
     file=file.path(save.dir, "sitedata.Rdata"))


bee.sum.screened <- bee %>%
    summarise(BeeRichness = length(unique(GenusSpecies[!is.na(GenusSpecies)])),
              BeeDiversity =
                  vegan::diversity(table(GenusSpecies[!is.na(GenusSpecies)]),
                                   index="shannon"),
              BeeAbund=sum(Count),
              BeeGenera=length(unique(Genus[!is.na(Genus)])),)

bee.sum.screened

floral.sum.screened <- floral %>%
    summarise(FlowerRichness =
                  length(unique(Flower_sci)),
              FlowerDiversity=
                  vegan:::diversity(table(Flower_sci),
                                    index="shannon"),
              BloomAbund=sum(Blooms, na.rm=TRUE))
floral.sum.screened

parasite.sum.screened <- indiv.data[indiv.data$ApidaeCtrl == 1,] %>%
    summarise(Apicystis = mean(Apicystis, na.rm=TRUE),
              Ascophaera= mean(Ascophaera, na.rm=TRUE),
              Crithidia= mean(Crithidia, na.rm=TRUE),
              AnyParasite= mean(AnyParasite, na.rm=TRUE))


parasite.sum.screened


## Stands total
print("Total stands with insect/plant data")
print(length(unique(site.data$Stand)))

print("Total stands with parasite data")
print(length(unique(indiv.data$Stand[!is.na(indiv.data$TestedTotals)])))

print("Total stands with offspring data")
print(length(unique(repro.block$Stand)))

print("Stands with nest boxes but no parasites")
unique(repro.block$Stand)[!unique(repro.block$Stand) %in%
                          unique(indiv.data$Stand[!is.na(indiv.data$TestedTotals)])]

print("Stands with parasite data but no nest box")
unique(indiv.data$Stand[!is.na(indiv.data$TestedTotals)])[!
    unique(indiv.data$Stand[!is.na(indiv.data$TestedTotals)])
    %in% unique(repro.block$Stand)
    ]
## the nest boxes were damaged by birds at these sites


print("Stands no nest box")
site.data$Stand[!site.data$Stand %in% unique(repro.block$Stand)]


bee.sum <- bee %>%
    group_by(Genus) %>%
    summarise(Species = length(unique(GenusSpecies[!is.na(GenusSpecies)])),
              Indiv = sum(Count))
