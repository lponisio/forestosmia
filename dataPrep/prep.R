rm(list=ls())
library(dplyr)
## prepares raw data and creates dataset for analyses

dir <- "~/Dropbox/forestosmia_saved"
parasite <- read.csv(file.path(dir,
                           "data/parasite.csv"),
                 stringsAsFactors=FALSE)

dbh <- read.csv(file.path(dir,
                               "data/dbh.csv"),
                     stringsAsFactors=FALSE, na.strings=c("","NA"))

floral <- read.csv(file.path(dir,
                               "data/floral.csv"),
                     stringsAsFactors=FALSE)

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
                          "data/beediversity.csv"),
                stringsAsFactors=FALSE)


## ************CLEAN PARASITE DATA AND ADD SUMMARY METRICS**************************************
source("../forestosmia/dataPrep/src/siteNames.R")

dim(parasite)

## drop the NOTES column
parasite$NOTES <- NULL

## drop all the no template controls under UniqueID because they were all neg for parasites
parasite <- parasite[!parasite$UniqueID == "C",]
dim(parasite)

## drop any sample where Apidae isnt 1
parasite <- parasite[!parasite$ApidaeCtrl == 0,]
dim(parasite)

## We may want to compare parasites in unemerged and foraging bees, in which case use "parasitecontrol" object
## However, for most analyses, we want to lose the unemerged, control bees, use "parasite"
parasitecontrol <- parasite
parasite <- parasitecontrol[!parasitecontrol$Stand == "control",]
dim(parasite)

## community health metrics
parasitenames <- c("Crithidia", "Ascophaera", "Apicystis")

parasite$ParasiteRichness <- rowSums(parasite[,parasitenames])

#### trues/falses times 1 (true is a 1 and false is 0)
parasite$AnyParasite <- (parasite$ParasiteRichness > 0)*1

## calculate avg sick individuals for each site
## do this with dplyr (as in ffar) instead?

sick.totals <- parasite %>%
  group_by(Stand) %>%
  summarise(TestedTotals = length(UniqueID),
            ParasitismRate=mean(AnyParasite, na.rm=TRUE),
            InfectedIndividuals=sum(AnyParasite, na.rm=TRUE))

dim(parasite)

## ************CLEAN VEG COVER DATA, ADD SUMMARY METRICS, MERGE TO PARASITE DATA*********************************

## calculate stand-level broadleaf cover (stand=site)

## keep only data for "broadleaf" trees, srubs, and forbs

## when dropping a variable, check the class with class(vegcover$broadleaf) and if its an integar dont put quotes around 0
vegcover <- vegcover[!vegcover$Broadleaf == 0,]

## replace "tr" (tracce) with 1.00
vegcover$perCover[vegcover$perCover ==
                 "tr"]  <- "1.0"

vegcover$perCover[vegcover$perCover ==
                    "TR"]  <- "1.0"

vegcover$perCover[vegcover$perCover ==
                    "Tr"]  <- "1.0"

## average abundance of broadleaf cover among plots at a site

vegcover$perCover<-as.numeric(vegcover$perCover)

BLcover.stand<-tapply(vegcover$perCover,vegcover$Stand,mean,na.rm=TRUE)

## add to parasite data as a new variable called stand intensity
## default is a left merge, merging the second object to the first object

## parasite$standintensity <- BLcover.stand[match(parasite$Stand,names(BLcover.stand))]

##look for any wonky things that didnt merge
## unique(parasite$Stand[is.na(parasite$standintensity)])

## #FIX THIS -- Luckiamute isnt in the vegcover data. What to do? They don't have data from that site, I asked Sara.

## ## drop Luckiamute
## parasite <- parasite[parasite$Stand != "Luckiamute", ]
## standinfo <- standinfo[standinfo$Stand != "Luckiamute", ]

## ## ************CLEAN STAND INFO DATA, MERGE *********************************

## merge  stand and stand intensity
standinfo$standintensity <- BLcover.stand[match(standinfo$Stand,
                                                names(BLcover.stand))]




## ## merge stand info data onto a new object called allspec
## ## default is a left merge, merging the second object to the first object

## allspecdata <- merge(parasite, standinfo)
## colnames(allspecdata)
## dim(allspecdata)
## sort(unique(allspecdata$Stand))

## ##look for any wonky things that didnt merge
## unique(allspecdata$Stand[is.na(allspecdata$CenterLat)])


## ************CLEAN DBH DATA, MERGE *********************************

## let's remove the columns that are not correctly filled out
dbh$LiveorDead <- NULL
dbh$SlopePerc <- NULL
dbh$Aspect <- NULL

## Right now blanks are assigned as NA, but a blank TreeNum means no trees, which we actually want to account for
## If there's an NA under TreeNum and DBH, make it zeroes
dbh$TreeNum[is.na(dbh$TreeNum)] <- 0
dbh$DBH[is.na(dbh$DBH)] <- 0

## we want tree species richness, tree species abundance, tree DBH averaged across trees at a site

# Tree DBH
dbh$DBH<-as.numeric(dbh$DBH)
DBH.stand<-tapply(dbh$DBH, dbh$Stand,mean,na.rm=TRUE)

## Tree Abundance:
dbh$TreeNum <- as.numeric(dbh$TreeNum)
TreeAbund.stand <- tapply(dbh$TreeNum,dbh$Stand,mean,na.rm=TRUE)

## Tree Richness: need to count distinct number of tree species at a site
## keep species as NA, potentially drop these rows or this will treat an NA as a tree species

dbh <- dbh[!is.na(dbh$Species), ]
TreeRichness.stand <- tapply(dbh$Species,dbh$Stand,function(x) length(unique(x)))

## then merge to dataset
standinfo$MeanDBH <- DBH.stand[match(standinfo$Stand,
                                       names(DBH.stand))]
standinfo$TreeAbund <- TreeAbund.stand[match(standinfo$Stand,
                                       names(TreeAbund.stand))]
standinfo$TreeRichness <- TreeRichness.stand[match(standinfo$Stand,
                                       names(TreeRichness.stand))]



## apply is by every row or every column.
## tapply is aggregated by one category.
## lapply is list apply.

## ## then merge to dataset
## allspecdata$DBH.stand <- DBH.stand[match(allspecdata$Stand,names(DBH.stand))]
## allspecdata$TreeAbund.stand <- TreeAbund.stand[match(allspecdata$Stand,names(TreeAbund.stand))]
## allspecdata$TreeRichness.stand <- TreeRichness.stand[match(allspecdata$Stand,names(TreeRichness.stand))]

## ## look for any wonky things that didnt merge for some reason(missing or NA data)
## unique(allspecdata$Stand[is.na(allspecdata$DBH.stand)])

## ************CLEAN FLORAL DATA, MERGE *********************************

## we want abundance flowers (abundflw), richness flowering species (richnessflwingsp), and abundfloweringsp

## let's drop the column "vials" since it's not correctly filled out
floral$Vials <- NULL
floral$Other_flowers <- NULL
floral$Notes <- NULL

# Floral Richness
# we don't want zeroes counted as 1 though
floral$Flower_sci[floral$Flower_sci==0] <- NA
FlowerRichness.stand<-tapply(floral$Flower_sci,floral$Stand,
                             function(x) length(unique(na.omit(x))))

#merge floral richness
standinfo$FlowerRichness <- FlowerRichness.stand[match(standinfo$Stand,
                                          names(FlowerRichness.stand))]

#take mean of many variabes together for each stand, this code doesnt take transect into account
#floral$Blooms<-as.numeric(floral$Blooms)
#floral$Stems<-as.numeric(floral$Stems)
#floral$Canopy_cent<-as.numeric(floral$Canopy_cent)

#floral.site <- aggregate(list(BloomAbund.stand=floral$Blooms,
#                           StemAbund.stand=floral$Stems,
#                           Canopy.stand=floral$Canopy_cent),
#                      list(Stand=floral$Stand),
#                      mean)

# take mean of floral variables across each stand, this code takes transect into account
## I DIDNT DO THIS FOR FLORAL RICHNESS THOUGH BUT I THINK THATS OK

floral.SR <- aggregate(list(BloomAbund.stand=floral$Blooms,
                              StemAbund.stand=floral$Stems,
                              Canopy.stand=floral$Canopy_cent),
                         list(Stand=floral$Stand,
                              SampleRound=floral$Transect),
                         sum, na.rm=TRUE)

## average the abundnce of floral features a sample round at a site
BloomAbund.SR.stand <- tapply(floral.SR$BloomAbund.stand,
                          floral.SR$Stand, mean,  na.rm=TRUE)

StemAbund.SR.stand <- tapply(floral.SR$StemAbund.stand,
                              floral.SR$Stand, mean,  na.rm=TRUE)

Canopy.stand.SR.stand <- tapply(floral.SR$Canopy.stand,
                              floral.SR$Stand, mean,  na.rm=TRUE)

#merge abundance per stand data (USE THIS IF we arent averaging across transect)
#dim(allspecdata)
#allspecdata <- merge(allspecdata, floral.site, all.x=TRUE)
#dim(allspecdata)

#merge abundance per stand data (USE THIS IF we average across transect)

standinfo$MeanBloomAbund <- BloomAbund.SR.stand[match(standinfo$Stand,
                                          names(BloomAbund.SR.stand))]

standinfo$MeanStemAbund <- StemAbund.SR.stand[match(standinfo$Stand,
                                         names(StemAbund.SR.stand))]

standinfo$MeanCanopy <- Canopy.stand.SR.stand[match(standinfo$Stand,
                                      names(Canopy.stand.SR.stand))]

dim(standinfo)

## look for any wonky things that didnt merge for some reason(missing or NA data)
unique(standinfo$Stand[is.na(standinfo$MeanBloomAbund)])
unique(standinfo$Stand[is.na(standinfo$FlowerRichness)])

## ************BEE DIVERSITY DATA, MERGE *********************************
## TBA WHICH ESTIMATE SHOULD BE USE FOR THE NESTS/PARASITES, THE MEAN OVER THE SEASON OR THE CLOSEST SAMPLE ROUND, OR THE AVERAGE OVER THE FRIST 2

bee$Caste <- NULL

## generate abundance of bees in each sample round
AbundBee.SR <- aggregate(list(Abund=bee$Count),
                           list(Stand=bee$Stand,
                                Round=bee$Round),
                           sum, na.rm=TRUE)

## average the abundnce of bees in a sample round at a site
AbundBee.SR.stand <- tapply(AbundBee.SR$Abund,
                        AbundBee.SR$Stand, mean,  na.rm=TRUE)

## richness of bees, get rid of zeroes
bee$GenSp[bee$GenSp=="0 0"] <- NA
RichBees.stand <- aggregate(list(BeeRichness=bee$GenSp),
                       list(Stand=bee$Stand),
                       function(x) length(unique(na.omit(x))))

# merge to allspecdata
standinfo$MeanAbundBee <- AbundBee.SR.stand[
    match(standinfo$Stand, names(AbundBee.SR.stand))]

#GET AN ERROR??! NOT SURE WHY
standinfo$RichnessBees <- RichBees.stand$BeeRichness[
                                match(standinfo$Stand,
                                RichBees.stand$Stand)]




## ************CLEAN REPRODUCTIVE DATA, MERGE *********************************

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

## *************************************************************
## HC mve to analysis

## standardize varaibles in our master dataset
## right now this doesnt work cause "RichBees.SR.stand" isnt merged on yet

## path.variables <- c("ParasiteRichness","AnyParasite","standintensity","Acres",
##                     "Age","Elev","CenterLat","CenterLon","DBH.stand","TreeAbund.stand",
##                     "TreeRichness.stand","FlowerRichness.stand","BloomAbund.SR.stand",
##                     "Stemabund.SR.stand","Canopy.SR.stand","AbundBee.SR.stand",
##                     "RichBees.SR.stand","femalesLaid.nest.stand","femalesLaidratio.nest.stand" )

## standardize <- function(x)
##   (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

## allspecdata[,path.variables] <- apply(allspecdata[,path.variables], 2,
##                                       standardize)

## *************************************************************
## merge parasite and stand info
dim(parasite)
parasite <- merge(parasite, standinfo)
dim(parasite)

## *************************************************************
## write out final data

write.csv(parasite, file=file.path(dir,
                     "cleanedData/parasite.csv"), row.names=FALSE)

write.csv(standinfo, file=file.path(dir,
                     "cleanedData/standinfo.csv"), row.names=FALSE)


write.csv(repro.nest, file=file.path(dir,
                     "cleanedData/NestRepro.csv"), row.names=FALSE)

save.dir <- "~/Dropbox/forestosmia"
save(parasite,
     file=file.path(save.dir, "data/parasite.Rdata"))

save(repro.nest,
     file=file.path(save.dir, "data/NestRepro.Rdata"))


#for analysis, do path analysis
#where is bee richness?
# stand itensity (broadlaef), dbh, flowers, bee richness -> parasitsm (any parasite, dont use parasite richness) -> reprodcutive output (fecundity of males and females, just try with females)
# stand itensity
# dont try to put the path analyses yet, do one glmm for parasitism, and do one glmm for reproduction, with parasitism in the model for reproduction
