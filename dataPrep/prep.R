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


## ************CLEAN PARASITE DATA AND ADD PARASITE SUMMARY METRICS**************************************

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


## ************CLEAN VEG COVER DATA**************************************

## fix site names to match between datasets
vegcover$Stand <- as.character(vegcover$Stand)
vegcover$Stand[vegcover$Stand ==
                 "Alexander Rd."]  <- "Alexander Rd"

## *************************************************************
## calculate site level characteristics for bees
## abundance (average between sample rounds)
abund.SR.bees <- aggregate(list(Abund=bees$no_individuals),
                           list(Site=bees$site,
                                SampleRound=bees$sample_pd),
                           sum)

abund.SR.bees <- tapply(abund.SR.bees$Abund,
                        abund.SR.bees$Site, mean)

## richness (total for a site)
site.bees <- aggregate(list(BeeRichness=bees$genus_sub_sp),
                       list(Site=bees$site),
                       function(x) length(unique(x)))

site.bees$BeeAbund <- abund.SR.bees[match(names(abund.SR.bees),
                                          site.bees$Site)]

## *************************************************************
## calculate site level characteristics for veg and merge with bee data

veg.site <- aggregate(list(AbundWoodyFlowers=veg$NoTreeShrubsFlower,
                           AbundAnnualFlowers=veg$NumTotalFlowers,
                           PlantRichness=veg$NumHerbPlantSpp,
                           PercentBareSoil=veg$PercentBareSoil),
                      list(Site=veg$Site),
                      mean)

site.char <- merge(veg.site, site.bees)

site.char <- merge(site.char,
                   unique(veg[, c("Site", "Size", "natural1000m",
                                  "natural2000m")]))

## *************************************************************
## calculate densities to control for garden size

site.char$BeeDensity <- site.char$BeeAbund/site.char$Size
site.char$BeeRichnessArea <- site.char$BeeRichness/site.char$Size

site.char$WoodyFlowerDensity <- site.char$AbundWoodyFlowers/site.char$Size
site.char$AnnualFlowerDensity <-
    site.char$AbundAnnualFlowers/site.char$Size
site.char$PlantRichnessArea <- site.char$PlantRichness/site.char$Size

## *************************************************************
## calculate total sick individuals for each site, bad thing

sick.totals <- aggregate(par.path[c(parasites, pathogens)],
                         list(Site=par.path$Site,
                              Genus=par.path$Genus),
                         sum, na.rm=TRUE)

tested.totals <- aggregate(par.path[c(parasites, pathogens)],
                           list(Site=par.path$Site,
                                Genus=par.path$Genus),
                           function(x) sum(!is.na(x)))

sick.totals$ScreenedPath <- tested.totals$CBPV
sick.totals$ScreenedPar <- tested.totals$Phorid

sick.totals[,parasites] <-
    sick.totals[,parasites]/sick.totals$ScreenedPar
sick.totals[,pathogens] <-
    sick.totals[,pathogens]/sick.totals$ScreenedPath

## *************************************************************
## standardize varaibles
path.variables <- c("WoodyFlowerDensity", "AnnualFlowerDensity",
                    "PlantRichnessArea", "natural1000m",
                    "natural2000m", "PercentBareSoil", "BeeDensity")

standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

site.char[,path.variables] <- apply(site.char[,path.variables], 2,
                                    standardize)

## *************************************************************
## merge pathogen and site data

dim(par.path)
## okay to drop controls
par.path$Site[!par.path$Site %in% site.char$Site]

par.path <- merge(par.path, site.char)
par.path$Date <- NULL
dim(par.path)

sick.totals <- merge(sick.totals, site.char)

## write out final data
write.csv(par.path, file=file.path(save.dir,
                                   "specimens-complete.csv"), row.names=FALSE)

save.dir.git <- "~/Dropbox/urbanbeeparasites/data"
save(par.path,
     site.char, sick.totals,
     file=file.path(save.dir.git, "specimens-complete.Rdata"))

