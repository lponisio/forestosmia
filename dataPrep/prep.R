rm(list=ls())
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
parasite$Stand <- as.character(par.path$Site)
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
                 "back grove"]  <- "Back Grove"
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

## drop all the no template controls under UniqueID because they were all neg for parasites
parasite <- parasite[!parasite$UniqueID == "C",]

## drop any sample where Apidae isnt 1
parasite <- parasite[!parasite$ApidaeCtrl == "0",]

## We may want to compare parasites in unemerged and foraging bees, in which case use "parasitecontrol" object
## However, for most analyses, we want to lose the unemerged, control bees, use "parasite"
parasitecontrol <- parasite
parasite <- parasitecontrol[!parasitecontrol$Stand == "control",]

## community health metrics
parasitenames <- c("Crithidia", "Ascophaera", "Apicystis")

parasite$ParasiteRichness <- rowSums(parasite[,parasitenames])
parasite$AnyParasite <- (parasite$ParasiteRichness > 0)*1

## calculate total sick individuals for each site

sick.totals <- aggregate(parasite[c(parasitenames)],
                         list(Stand=parasite$Stand,
                              UniqueID=parasite$UniqueID),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(parasite[c(parasitenames)],
                           list(Stand=parasite$Stand,
                                UniqueID=parasite$UniqueID),
                           function(x) sum(!is.na(x)))

sick.totals$ScreenedPar <- tested.totals$Crithidia

sick.totals[,parasitenames] <-
  sick.totals[,parasitenames]/sick.totals$ScreenedPar

parasite<- merge(sick.totals, parasite)

## ************CLEAN VEG COVER DATA, ADD SUMMARY METRICS, MERGE TO PARASITE DATA*********************************

## fix site names to match between datasets
vegcover$Stand <- as.character(vegcover$Stand)
vegcover$Stand[vegcover$Stand ==
                 "Alexander Rd."]  <- "Alexander Rd"

## calculate stand-level broadleaf cover (stand=site)
## 
## keep only data for "broadleaf" trees, srubs, and forbs
vegcover <- vegcover[!vegcover$Broadleaf == "0",]

## replace "tr" (tracce) with .01
vegcover$perCover[vegcover$perCover ==
                 "tr"]  <- "0.01"

vegcover$perCover[vegcover$perCover ==
                    "TR"]  <- "0.01"

## average abundance of broadleaf cover among plots at a site

BLcover.site <- aggregate(list(BLcover=vegcover$perCover), 
                          list(Stand=vegcover$Stand),
                          mean)
                          

## merge with parasite data. 

allspecdata <- merge(BLcover.site, parasite)


## ************CLEAN STAND INFO DATA, MERGE *********************************
## merge stand info data onto allspec

##WHY DID I LOSE SO MANY UNIQUE IDs? I should have like 366?

allspecdata <- merge(allspecdata, standinfo)
colnames(allspecdata)
allspecdata$NOTES <- NULL
colnames(allspecdata)

## ************CLEAN DBH DATA, MERGE *********************************

## let's remove the column "live or dead" since it's not correctly filled out
dbh$LiveorDead <- NULL

## we want tree species richness, tree species abundance, tree DBH averaged across plots at a site

## then merge to dataset

## ************CLEAN FLORAL DATA, MERGE *********************************

## we want abundance flowers (abundflw), richness flowering species (richnessflwingsp), and abundfloweringsp 

## let's drop the column "vials" since it's not correctly filled out
floral$Vials <- NULL

## some ID data missing, I emailed Jim and Sara about how to proceed. 
## ...Can each observation be assumed to be a distinct plant species?

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
