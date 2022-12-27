setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')

## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created. For the
## models for Osmia reproduction, missing parasite prevalence rate
## data is imputed before model fitting.

setwd("analyses")
rm(list=ls())

source("src/init.R")
load("../data/indivdata.Rdata")
load("../data/reproblock.Rdata")

source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/misc.R")


## set to the number of cores you would like the models to run on
ncores <- 10


## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## order data by stand for ease of viewing
indiv.data <- indiv.data[order(indiv.data$Stand),]
repro.block <- repro.block[order(repro.block$Stand),]

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("FlowerDiversity",
          "MeanBloomAbund",
          "MeanBeeDiversity",
          "MeanBeeAbund",
          "MeanHBAbund",
          "Age_LandTrendr",
          "Acres",
          "Elev")

## log variables for age and acres
indiv.data$Age_LandTrendr <- log(indiv.data$Age_LandTrendr)
indiv.data$Acres <- log(indiv.data$Acres)
repro.block$Age_LandTrendr <- log(repro.block$Age_LandTrendr)
repro.block$Acres <- log(repro.block$Acres)

##  center all of the x variables across the datasets
indiv.data[, vars] <- apply(indiv.data[, vars], 2, standardize)
repro.block[, vars] <- apply(repro.block[, vars], 2, standardize)

repro.block$ParasitismRate <- standardize(repro.block$ParasitismRate)

## create a dummy varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms

indiv.data <- makeDataMultiLevel(indiv.data, "Stand")
repro.block <- makeDataMultiLevel(repro.block, "Stand")
indiv.data$Owner <- factor(indiv.data$Owner,
                           levels= c("OwnerA", "OwnerB", "OwnerC", "ODF"))

## create a dummy varaible "WeightPar" for the parasite data. The
## intention is to keep stan from dropping data for site-level models,
## but weight is 0 for parasite models.
indiv.data$WeightsPar <- 1
indiv.data$WeightsPar[is.na(indiv.data$AnyParasite)] <- 0

## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
## to keep it in the models
indiv.data$AnyParasite[is.na(indiv.data$AnyParasite)] <- 0

## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## define all the formulas for the different parts of the models

## flower diversity
formula.flower.div <- formula(FlowerDiversity | weights(Weights) ~
                                  + Age_LandTrendr + I(Age_LandTrendr^2) +
                                      Elev + Owner
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund | weights(Weights) ~
                                        Age_LandTrendr  + Elev + Owner
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************

## bee diversity
formula.bee.div <- formula(MeanBeeDiversity | weights(Weights)~
                               FlowerDiversity +
                                   Age_LandTrendr + Acres)

## bee abund
formula.bee.abund <- formula(MeanBeeAbund | weights(Weights)~
                                 MeanBloomAbund +
                                     Age_LandTrendr + Acres)

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************

formula.parasite <- formula(AnyParasite | weights(WeightsPar) ~
                                MeanBeeAbund +
                                    FlowerDiversity +
                                    MeanBeeDiversity +
                                    MeanBloomAbund +
                                    (1|Stand)
                            )

formula.parasite.site <- formula(scale(ParasitismRate) | weights(Weights) ~
                                     MeanBeeAbund +
                                     MeanBeeDiversity +
                                      MeanBloomAbund +
                                      FlowerDiversity
                                 )

## convert to brms format
bf.par.site <- bf(formula.parasite.site)
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)

## **********************************************************
## Model 1 community effects on bee parasitism
## **********************************************************

bf.par <- bf(formula.parasite, family="bernoulli")

## full model
bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.par +
    set_rescor(FALSE)

## run model
fit <- brm(bform, indiv.data,
           cores=ncores,
           iter = 10^4,
           chains =3,
           thin=1,
           init=0,
           control = list(adapt_delta = 0.99))

write.ms.table(fit, "parasitism_poly")
save(fit, indiv.data,
     file="saved/parasiteFitMod.Rdata")

## dignostic figures
plot.res(fit, "parasite")

## ************************************************************
## Model 2.1 formulas for the site effects on nest reproduction
## ************************************************************

## can use total offspring or total female offspring as the response
## variable

## data  imputation
imp.dat.pre <- mice(repro.block, m = 2,
                    print = FALSE)

pred <- imp.dat.pre$pred
pred[, c("Weights", "Block", "SiteIDs", "Year")] <- 0
pred[, c("Stand", "Owner")] <- 1

meth <- imp.dat.pre$meth
meth[c(2)] <- "rf"

imp.dat <- mice(repro.block, m = 100, predictorMatrix=pred,
                method=meth,
                print = FALSE)


densityplot(imp.dat)
plot(imp.dat, c("ParasitismRate"))

ys <- c("SumOffspring")

xvar.NestRepro <- c("scale(ParasitismRate)",
    "MeanBloomAbund",
    "FlowerDiversity",
    "MeanBeeAbund")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)",
                         sep="+")))
})

bf.offspring <- bf(formulas.NestRepro[[1]], family="negbinomial")

## full model
bform2 <-    bf.par.site +
    bf.offspring +
    set_rescor(FALSE)

## run model
fit2 <- brm_multiple(bform2, data=imp.dat,
           cores=ncores,
           iter = 5*10^4,
           chains = 2,
           inits=0,
           thin=2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15))

write.ms.table(fit2, "offspring")

save(fit2, repro.block,
     file="saved/offspringFitMod.Rdata")

## diagnostic plots
mcmc_trace(fit2)
ggsave("figures/diagnostics/offspring.pdf",
       height=11, width=8.5)
