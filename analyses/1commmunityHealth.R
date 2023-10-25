## set your working directory to the forestosmia git repo
setwd('forestosmia')
setwd("analyses")
## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created. For the
## models for Osmia reproduction, missing parasite prevalence rate
## data is imputed before model fitting.

rm(list=ls())

source("src/init.R")
load("../data/indivdata.Rdata")
load("../data/reproblock.Rdata")

source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/misc.R")
source("src/runParasiteModels.R")

## set to the number of cores you would like the models to run on
ncores <- 3

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## order data by stand for ease of viewing
indiv.data <- indiv.data[order(indiv.data$Stand),]
repro.block <- repro.block[order(repro.block$Stand),]

## drop 1 outlier max MeanBeeAbund? Test based on reviewer
## suggestion.
max.abund.par.prep <- indiv.data[!is.na(indiv.data$AnyParasite),]
max.abund.par <-  unique(max.abund.par.prep$Stand[
  max.abund.par.prep$MeanBeeAbund ==
  max(max.abund.par.prep$MeanBeeAbund)])

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("FlowerDiversity",
          "MeanBloomAbund",
          "MeanBeeDiversity",
          "MeanBeeAbund",
          "MeanHBAbund",
          "Age_LandTrendr",
          "Hectares",
          "Elev")

## log variables for age and hectares
indiv.data$Age_LandTrendr <- log(indiv.data$Age_LandTrendr)
indiv.data$Hectares <- log(indiv.data$Hectares)
repro.block$Age_LandTrendr <- log(repro.block$Age_LandTrendr)
repro.block$Hectares <- log(repro.block$Hectares)

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
                           levels= c("ODF", "OwnerA", "OwnerB", "OwnerC"))

## create a dummy varaible "WeightPar" for the parasite data. The
## intention is to keep stan from dropping data for site-level models,
## but weight is 0 for parasite models.
indiv.data$WeightsPar <- 1
indiv.data$WeightsPar[is.na(indiv.data$AnyParasite)] <- 0

## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
## to keep it in the models
indiv.data$AnyParasite[is.na(indiv.data$AnyParasite)] <- 0
indiv.data$Crithidia[is.na(indiv.data$Crithidia)] <- 0
indiv.data$Apicystis[is.na(indiv.data$Apicystis)] <- 0
indiv.data$Ascophaera[is.na(indiv.data$Ascophaera)] <- 0

## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(FlowerDiversity | weights(Weights) ~
                                  + Age_LandTrendr + I(Age_LandTrendr^2) +
                                      Elev + Owner
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund | weights(Weights) ~
                                        Age_LandTrendr  + Elev
                                +  I(Age_LandTrendr^2) + Owner
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(MeanBeeDiversity | weights(Weights)~
                               FlowerDiversity +  I(Age_LandTrendr^2)+
                                   Age_LandTrendr + Hectares)

## bee abund
formula.bee.abund <- formula(MeanBeeAbund | weights(Weights)~
                                 MeanBloomAbund +  I(Age_LandTrendr^2)+
                                     Age_LandTrendr + Hectares)

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************
formula.parasite.site <- formula(scale(ParasitismRate) | weights(Weights) ~
                                     MeanBeeAbund +
                                     MeanBeeDiversity +
                                      MeanBloomAbund
                                 )

## convert to all models brms format
bf.par.site <- bf(formula.parasite.site)
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)

## **********************************************************
## Model 1 community effects on bee parasitism
## **********************************************************
## create a block-stand identified to be used as a random effect
indiv.data$StandBlock <- paste0(indiv.data$Stand, indiv.data$Block)

xvars.parasites <- c("MeanBeeAbund", "MeanBeeDiversity",
                     "MeanBloomAbund",
                     "(1|Stand)", "(1|StandBlock)")

## run models with all parasites simultaneously 
fit.combined.parasite <- runCombinedParasiteModels(indiv.data, "osmia_all_data",
                                                   parasites=c("Crithidia",
                                                               "Apicystis",
                                                               "Ascophaera"),
                                                   xvars.parasites,
                                                   iter = 10^4,
                                                   chains = 3)

## run model with "any parasite" i.e., did a bee test positive for any
## of the three parasites
fit.anyParasite <- runParasiteModels(indiv.data, "osmia_all_data",
                                     "AnyParasite",
                                     xvars.parasites,
                                     SEM=TRUE,
                                     iter = 10^4,
                                     chains = 3)

## ************************************************************
## Model 2.1 formulas for the site effects on nest reproduction
## ************************************************************

## can use total offspring or total female offspring as the response
## variable, both yield similar results. 
## data  imputation
imp.dat.pre <- mice(repro.block, m = 2,
                    print = FALSE)

pred <- imp.dat.pre$pred
## bookkeeping varaibles that should not be included
pred[, c("Weights", "Block", "SiteIDs", "Year")] <- 0
pred[, c("Stand", "Owner")] <- 1

meth <- imp.dat.pre$meth
meth[c(2)] <- "rf"

imp.dat <- mice(repro.block, m = 100, predictorMatrix=pred,
                method=meth,
                print = FALSE)

densityplot(imp.dat)
quartz()
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
                     init=0,
                     thin=2,
                     open_progress = FALSE,
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 15),
                     combine = FALSE)

fit2.raw <- fit2
fit2 <- combine_models(mlist=fit2, check_data=FALSE)

write.ms.table(fit2, "offspring")

## calculate a pooled R2 combining the models from the different
## imputed datasets
pool_R2 <- function(mlist, probs = c(0.025, 0.975), robust = FALSE, ...) {
  r2_post <- sapply(mlist, bayes_R2, summary = FALSE, ..., simplify = TRUE)
  posterior_summary(c(r2_post), probs = probs, robust = robust)
}

r2.offspring <- pool_R2(fit2.raw)
round(r2.offspring, 2)

save(fit2, repro.block, 
     fit2.raw, r2.offspring, 
     file="saved/offspringFitMod.Rdata")


## diagnostic plots
mcmc_trace(fit2)
ggsave("figures/diagnostics/offspring.pdf",
       height=11, width=8.5)

