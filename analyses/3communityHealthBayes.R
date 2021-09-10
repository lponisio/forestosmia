setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')
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

ncores <- 4

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

indiv.data <- indiv.data[order(indiv.data$Stand),]
repro.block <- repro.block[order(repro.block$Stand),]

## create a dumby varaible to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms

## flower richness

vars <- c("TreeRichness",
          "FlowerDiversity",
          "MeanBloomAbund",
          "MeanBeeDiversity",
          "MeanBeeAbund",
          "Age"
          ## "RoseDiversity",
          ## "RoseMeanBloomAbund"
          )

site.data[, vars] <- apply(site.data[, vars], 2, scale)
indiv.data[, vars] <- apply(indiv.data[, vars], 2, scale)
repro.block[, vars] <- apply(repro.block[, vars], 2, scale)
repro.block$ParasitismRate <- scale(log(repro.block$ParasitismRate))

repro.block <- makeDataMultiLevel(repro.block, "Stand")
indiv.data <- makeDataMultiLevel(indiv.data, "Stand")


## flower diversity
formula.flower.div <- formula(FlowerDiversity | weights(Weights) ~
                                  ## (ManagementIntensity) +
                                  s(Age)
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund | weights(Weights) ~
                                    ## (ManagementIntensity) +
                                        s(Age)
                                )

## **********************************************************
## formula for forest effects on bee community
## **********************************************************

## bee diversity
formula.bee.div <- formula(MeanBeeDiversity | weights(Weights)~
                               ## (ManagementIntensity) +
                               (MeanBloomAbund) +
                               (FlowerDiversity) +
                                   s(Age)

                           )

## bee abund
formula.bee.abund <- formula(MeanBeeAbund | weights(Weights)~
                                 ## (ManagementIntensity) +
                                 (MeanBloomAbund) +
                                 (FlowerDiversity) +
                                     s(Age)

                             )

## **********************************************************
## formula for bee community effects on parasitism
## **********************************************************

formula.parasite <- formula(AnyParasite ~
                                MeanBeeAbund +
                                FlowerDiversity +
                                MeanBeeDiversity +
                                MeanBloomAbund +
                                (1|Stand)
                            )

formula.parasite.site <- formula(ParasitismRate | weights(Weights) ~
                                     (MeanBeeAbund) +
                                     (MeanBeeDiversity) +
                                     (MeanBloomAbund) +
                                     (FlowerDiversity)
                                 )


## **********************************************************
## community models
## **********************************************************
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)

prior <- c(set_prior("normal(0,1)", class="b"))

## **********************************************************
## SEM parasitism
## **********************************************************

bf.par <- bf(formula.parasite, family="bernoulli")

bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.par +
    set_rescor(FALSE)


fit <- brm(bform, indiv.data,
           cores=ncores,
           iter = 10^5,
           chains = 4,
           thin=2,
           inits=0,
           ## prior=prior,
           control = list(adapt_delta = 0.99))

write.ms.table(fit, "parasitism")

mcmc_trace(fit)
ggsave("figures/diagnostics/parasite.pdf",
       height=11, width=8.5)

save(fit, site.data, indiv.data,
     file="saved/parasiteFitMod.Rdata")

## **********************************************************
## formulas for the site effects on nest reproductiom
## **********************************************************

ys <- c("SumOffspring",
        "Females")

bf.par.site <- bf(formula.parasite.site)


xvar.NestRepro <- c("ParasitismRate",
                    ## "(ManagementIntensity)",
                    "(MeanBloomAbund)",
                    "(FlowerDiversity)",
                    "(MeanBeeAbund)",
                    "(MeanBeeDiversity)")

                    ## "RoseMeanBloomAbund", ## no effect
                    ## "RoseDiversity"


formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)",
                         sep="+")))
})

bf.offspring <- bf(formulas.NestRepro[[1]], family="poisson")

## bform2 <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv +  bf.par.site +
##     bf.offspring +
##     set_rescor(FALSE)

bform2 <-   bf.par.site +
    bf.offspring +
    set_rescor(FALSE)


fit2 <- brm(bform2, repro.block,
           cores=ncores,
           iter = 10^5,
           chains = 4,
           inits=0,
           ## prior=prior,
           thin=2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15))

write.ms.table(fit2, "offspring")


mcmc_trace(fit2)
ggsave("figures/diagnostics/offspring.pdf",
       height=11, width=8.5)

save(fit2, repro.block,
     file="saved/offspringFitMod.Rdata")
