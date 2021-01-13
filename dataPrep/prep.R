rm(list=ls())
library(dplyr)
## prepares raw data and creates dataset for analyses

dir <- "~/Dropbox/forestosmia_saved"
parasite <- read.csv(file.path(dir,
                           "data/parasite.csv"),
                 stringsAsFactors=FALSE)

dbh <- read.csv(file.path(dir,
                               "data/dbh.csv"),
                     stringsAsFactors=FALSE)

floral <- read.csv(file.path(dir,
                               "data/floral.csv"),
                     stringsAsFactors=FALSE)

reproductive <- read.csv(file.path(dir,
                               "data/reproductive.csv"),
                     stringsAsFactors=FALSE)

standinfo <- read.csv(file.path(dir,
                               "data/standinfo.csv"),
                     stringsAsFactors=FALSE)

vegcover <- read.csv(file.path(dir,
                                "data/vegcover.csv"),
                      stringsAsFactors=FALSE)


## ************CLEAN PARASITE DATA AND ADD SUMMARY METRICS**************************************

## fix site names to match between datasets
parasite$Stand <- as.character(parasite$Stand)
parasite$Stand[parasite$Stand ==
              "norton hill"]  <- "Norton Hill"
parasite$Stand[parasite$Stand ==
                 "Norton hill"]  <- "Norton Hill"
parasite$Stand[parasite$Stand ==
                 "buttermilk creek"]  <- "Buttermilk Creek"
parasite$Stand[parasite$Stand ==
                 "lake lyons"]  <- "Lake Lyons"
parasite$Stand[parasite$Stand ==
                 "bethel vernon"]  <- "Bethel Vernon"
parasite$Stand[parasite$Stand ==
                 "Bethel vernon"]  <- "Bethel Vernon"
parasite$Stand[parasite$Stand ==
                 "luckiamute"]  <- "Luckiamute"
parasite$Stand[parasite$Stand ==
                 "luckiamate"]  <- "Luckiamute"
parasite$Stand[parasite$Stand ==
                 "Luckiamate"]  <- "Luckiamute"
parasite$Stand[parasite$Stand ==
                 "broken horn"]  <- "Broken Horn"
parasite$Stand[parasite$Stand ==
                 "black rock"]  <- "Black Rock"
parasite$Stand[parasite$Stand ==
                 "haybarn"]  <- "Haybarn"
parasite$Stand[parasite$Stand ==
                 "cougar ridge"]  <- "Cougar Ridge"
parasite$Stand[parasite$Stand ==
                 "miller creek"]  <- "Miller Creek"
parasite$Stand[parasite$Stand ==
                 "wild cab"]  <- "Wild Cab"
parasite$Stand[parasite$Stand ==
                 "gage"]  <- "Gage"
parasite$Stand[parasite$Stand ==
                 "elk city"]  <- "Elk City"
parasite$Stand[parasite$Stand ==
                 "south fork"]  <- "South Fork"
parasite$Stand[parasite$Stand ==
                 "hull oakes"]  <- "Hull Oakes"
parasite$Stand[parasite$Stand ==
                 "walker camp"]  <- "Walker Camp"
parasite$Stand[parasite$Stand ==
                 "hatchery falls"]  <- "Hatchery Falls"
parasite$Stand[parasite$Stand ==
                 "bow season"]  <- "Bow Season"
parasite$Stand[parasite$Stand ==
                 "sheythe"]  <- "Sheythe"
parasite$Stand[parasite$Stand ==
                 "Reese creek"]  <- "Reese Creek"
parasite$Stand[parasite$Stand ==
                 "back grove"]  <- "Backgrove"
parasite$Stand[parasite$Stand ==
                 "bald mountain"]  <- "Bald Mountain"
parasite$Stand[parasite$Stand ==
                 "snag"]  <- "Snag"
parasite$Stand[parasite$Stand ==
                 "logsden"]  <- "Logsden"
parasite$Stand[parasite$Stand ==
                 "Paul mort."]  <- "Paul Mortenson"
parasite$Stand[parasite$Stand ==
                 "Cougar ridge"]  <- "Cougar Ridge"
parasite$Stand[parasite$Stand ==
                 "Bigrock"]  <- "Big Rock Creek"
parasite$Stand[parasite$Stand ==
                 "Buttermilk"]  <- "Buttermilk Creek"
parasite$Stand[parasite$Stand ==
                 "Broken horn"]  <- "Broken Horn"
parasite$Stand[parasite$Stand ==
                 "Wolf cab"]  <- "Wolf Cabin"
parasite$Stand[parasite$Stand ==
                 "wolf cab"]  <- "Wolf Cabin"
parasite$Stand[parasite$Stand ==
                 "Wolf cabin"]  <- "Wolf Cabin"
parasite$Stand[parasite$Stand ==
                 "Miller creek"]  <- "Miller Creek"
parasite$Stand[parasite$Stand ==
                 "Walker camp"]  <- "Walker Camp"
parasite$Stand[parasite$Stand ==
                 "Lake lyons"]  <- "Lake Lyons"
parasite$Stand[parasite$Stand ==
                 "Gage "]  <- "Gage"

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
## do this with dplyr instead, check out ffar code (summary) - TO DO

sick.totals <- parasite %>%
  group_by(Stand) %>%
  summarise(TestedTotals = length(UniqueID),
            ParasitismRate=mean(AnyParasite, na.rm=TRUE),
            InfectedIndividuals=sum(AnyParasite, na.rm=TRUE))

dim(parasite)

## ************CLEAN VEG COVER DATA, ADD SUMMARY METRICS, MERGE TO PARASITE DATA*********************************

## fix site names to match between datasets
vegcover$Stand <- as.character(vegcover$Stand)
vegcover$Stand[vegcover$Stand ==
                 "Alexander Rd."]  <- "Alexander Rd"
vegcover$Stand[vegcover$Stand ==
                 "luckiamate"]  <- "Luckiamute"
vegcover$Stand[vegcover$Stand ==
                 "Luckiamate"]  <- "Luckiamute"
vegcover$Stand[vegcover$Stand ==
                 "Lake lyons"]  <- "Lake Lyons"
vegcover$Stand[vegcover$Stand ==
                 "Hull Oaks"]  <- "Hull Oakes"
vegcover$Stand[vegcover$Stand ==
                 "Walker camp"]  <- "Walker Camp"


## calculate stand-level broadleaf cover (stand=site)
## 
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


## merge with parasite data. default is a left merge, merging the second obect to the first object 

parasite$standintensity <- BLcover.stand[match(parasite$Stand,names(BLcover.stand))]

##look for any wonky things that didnt merge
unique(parasite$Stand[is.na(parasite$standintensity)])


#FIX THIS -- Luckiamute isnt in the vegcover data. What to do? I emailed Sara


## ************CLEAN STAND INFO DATA, MERGE *********************************
## merge stand info data onto a new object called allspec

allspecdata <- merge(parasite, standinfo)
colnames(allspecdata)
dim(allspecdata)

## ************CLEAN DBH DATA, MERGE *********************************

## let's remove the column "live or dead" since it's not correctly filled out
dbh$LiveorDead <- NULL

## TO DO? assign zeoes and NAs under blanks??

## we want tree species richness, tree species abundance, tree DBH averaged across trees at a site

dbh$DBH<-as.numeric(dbh$DBH)
DBH.stand<-tapply(dbh$DBH,dbh$Stand,mean,na.rm=TRUE)

dbh$TreeNum<-as.numeric(dbh$TreeNum)
TreeAbund.stand<-tapply(dbh$TreeNum,dbh$Stand,mean,na.rm=TRUE)

## THIS DIDNT WORK, MAYBE BECAUSE BLANKS? need to count distinct number of tree species at a site
TreeRichness.stand<-apply(dbh$Species,dbh$Stand,function(x) length(unique(x)), na.rm=TRUE)

## then merge to dataset
allspecdata$DBH.stand <- DBH.stand[match(allspecdata$Stand,names(DBH.stand))]
allspecdata$TreeAbund.stand <- TreeAbund.stand[match(allspecdata$Stand,names(TreeAbund.stand))]
allspecdata$TreeRichness.stand <- TreeRichness.stand[match(allspecdata$Stand,names(TreeRichness.stand))]

## look for any wonky things that didnt merge for some reason(missing or NA data)
unique(allspecdata$Stand[is.na(allspecdata$DBH.stand)])

## ************CLEAN FLORAL DATA, MERGE *********************************

## we want abundance flowers (abundflw), richness flowering species (richnessflwingsp), and abundfloweringsp 

## let's drop the column "vials" since it's not correctly filled out
floral$Vials <- NULL
floral$Other_flowers <- NULL
floral$Notes <- NULL

#Floral bloom at a stand across transects

floral$Blooms<-as.numeric(floral$Blooms)
Blooms.stand<-tapply(floral$Blooms,floral$Stand,mean,na.rm=TRUE)



## ************CLEAN REPRODUCTIVE DATA, MERGE *********************************

## want sex ratio, offspring production, foraging trip length

## ************BEE DIVERSE DATA, MERGE *********************************

## I don't have this data yet from them? will calculate wild bee richness and wild bee abundance



## *************************************************************
## standardize varaibles in our master dataset, 
path.variables <- c("BLcover","Acres")

standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

allspecdata[,path.variables] <- apply(allspecdata[,path.variables], 2,
                                      standardize)

## *************************************************************
## write out final data
write.csv(par.path, file=file.path(save.dir,
                                   "specimens-complete.csv"), row.names=FALSE)

save.dir.git <- "~/Dropbox/forestosmia/data"
save(allspecdata,
     file=file.path(save.dir.git, "specimens-complete.Rdata"))


#for analysis, do path analysis
#where is bee richness?
# stand itensity (broadlaef), dbh, flowers, bee richness -> parasitsm (any parasite, dont use parasite richness) -> reprodcutive output (fecundity of males and females, just try with females)
# stand itensity 
# dont try to put the path analyses yet, do one glmm for parasitism, and do one glmm for reproduction, with parasitism in the model for reproduction