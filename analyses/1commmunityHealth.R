setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')

## setwd('~/Dropbox (University of Oregon)/forestosmia')

setwd("analyses")
rm(list=ls())
library(ggplot2)
library(viridis)
library(brms)
library(bayesplot)
library(tidybayes)
library(tidyverse)

load("../data/indivdata.Rdata")
load("../data/sitedata.Rdata")
load("../data/reproblock.Rdata")

source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/misc.R")

ncores <- 10

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## order data by stand for ease of viewing
indiv.data <- indiv.data[order(indiv.data$Stand),]
repro.block <- repro.block[order(repro.block$Stand),]

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("TreeRichness",
          "FlowerDiversity",
          "MeanBloomAbund",
          "MeanBeeDiversity",
          "MeanBeeAbund",
          "Age",
          "Elev")

##  center all of the x variables across the datasets
site.data[, vars] <- apply(site.data[, vars], 2, standardize)
indiv.data[, vars] <- apply(indiv.data[, vars], 2, standardize)
repro.block[, vars] <- apply(repro.block[, vars], 2, standardize)

## standardize paristism rate, take the log first to help with skew
repro.block$ParasitismRate <- standardize(log(repro.block$ParasitismRate))

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms

repro.block <- makeDataMultiLevel(repro.block, "Stand")
indiv.data <- makeDataMultiLevel(indiv.data, "Stand")

indiv.data$Owner <- factor(indiv.data$Owner,
                    levels= c("Hancock", "Starker", "ODF", "Weyco"))


## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.

indiv.data$WeightsPar <- 1
indiv.data$WeightsPar[is.na(indiv.data$AnyParasite)] <- 0

## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
## to keep it in the models, this is commented out because we don't
## want to use the entire dataset in this publication

## indiv.data$AnyParasite[is.na(indiv.data$AnyParasite)] <- 0


## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## define all the formulas for the different parts of the models

## flower diversity
formula.flower.div <- formula(FlowerDiversity | weights(Weights) ~
                                  Owner + Age + Elev
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund | weights(Weights) ~
                                    Owner +  Age  + Elev
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************

## bee diversity
formula.bee.div <- formula(MeanBeeDiversity | weights(Weights)~
                               MeanBloomAbund +
                               FlowerDiversity +
                               Age)

                           )

## bee abund
formula.bee.abund <- formula(MeanBeeAbund | weights(Weights)~
                                 MeanBloomAbund +
                                 FlowerDiversity +
                                 Age)

                             )

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************

formula.parasite <- formula(AnyParasite | weights(WeightsPar) ~
                                MeanBeeAbund*
                                FlowerDiversity +
                                MeanBeeDiversity +
                                MeanBloomAbund +
                                (1|Stand)
                            )

formula.parasite.site <- formula(ParasitismRate | weights(Weights) ~
                                     MeanBeeAbund +
                                     MeanBeeDiversity +
                                     MeanBloomAbund +
                                     (FlowerDiversity)
                                 )


## **********************************************************
## Community models
## **********************************************************
## convert to brms format

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
           chains = 2,
           thin=1,
           inits=0,
           control = list(adapt_delta = 0.99))

write.ms.table(fit, "parasitism_interaction")

save(fit, site.data, indiv.data,
     file="saved/parasiteFitMod.Rdata")

## dignostic figures
mcmc_trace(fit)
ggsave("figures/diagnostics/parasite.pdf",
       height=11, width=8.5)

## ************************************************************
## Model 2.1 formulas for the site effects on nest reproduction
## ************************************************************

## can use total offspring or total female offspring as the response
## variable

ys <- c("SumOffspring",
        "Females")

bf.par.site <- bf(formula.parasite.site)

xvar.NestRepro <- c("ParasitismRate",
                   "(MeanBloomAbund)",
                    "(FlowerDiversity)",
                    "(MeanBeeAbund)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)",
                         sep="+")))
})

bf.offspring <- bf(formulas.NestRepro[[1]], family="poisson")

## full model
bform2 <-   bf.par.site +
    bf.offspring +
    set_rescor(FALSE)


## run model
fit2 <- brm(bform2, repro.block,
           cores=ncores,
           iter = 10^4,
           chains = 4,
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
