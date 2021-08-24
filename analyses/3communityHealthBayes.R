setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(car)
library(ggplot2)
library(viridis)
library(brms)

load("../data/indivdata.Rdata")
load("../data/sitedata.Rdata")
load("../data/reproblock.Rdata")

ncores <- 10

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

indiv.data <- indiv.data[order(indiv.data$Stand),]
repro.block <- repro.block[order(repro.block$Stand),]

## create a dumby varaible to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms

site.ids <- unlist(tapply(indiv.data$Stand,
                          indiv.data$Stand,
                          function(x) 1:length(x)))
names(site.ids) <- NULL
indiv.data$SiteIDs <- site.ids
indiv.data$Weights <- indiv.data$SiteIDs
indiv.data$Weights[indiv.data$Weights > 1] <- 0


site.ids <- unlist(tapply(repro.block$Stand,
                          repro.block$Stand,
                          function(x) 1:length(x)))
names(site.ids) <- NULL
repro.block$SiteIDs <- site.ids
repro.block$Weights <- repro.block$SiteIDs
repro.block$Weights[repro.block$Weights > 1] <- 0

## flower richness

vars <- c("BLcover",
          "FlowerDiversity",
          "MeanBloomAbund",
          "MeanBeeDiversity",
          "MeanBeeAbund",
          "AgePoly1",
          "AgePoly2",
          "AgePoly3",
          "AgePoly4"
          )

site.data[, vars] <- apply(site.data[, vars], 2, scale)
indiv.data[, vars] <- apply(indiv.data[, vars], 2, scale)
repro.block[, vars] <- apply(repro.block[, vars], 2, scale)
repro.block$ParasitismRate <- scale(repro.block$ParasitismRate)


## flower diversity
formula.flower.div <- formula(FlowerDiversity | weights(Weights) ~
                                  (BLcover) +
                                  (AgePoly1) +
                                  (AgePoly2) +
                                  (AgePoly3) +
                                  (AgePoly4)
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund | weights(Weights) ~
                                    (BLcover) +
                                    (AgePoly1) +
                                    (AgePoly2) +
                                    (AgePoly3) +
                                    (AgePoly4)
                                )

## **********************************************************
## formula for forest effects on bee community
## **********************************************************

## bee diversity
formula.bee.div <- formula(MeanBeeDiversity | weights(Weights)~
                               (BLcover) +
                               (MeanBloomAbund) +
                               ## (MeanBeeAbund) +
                               (FlowerDiversity) +
                               (AgePoly1) +
                               (AgePoly2)

                           )

## bee abund
formula.bee.abund <- formula(MeanBeeAbund | weights(Weights)~
                                 (BLcover) +
                                 (MeanBloomAbund) +
                                 (FlowerDiversity) +
                                 (AgePoly1) +
                                 (AgePoly2)
                             )
## **********************************************************
## formula for bee community effects on parasitism
## **********************************************************

formula.parasite <- formula(AnyParasite ~
                                     (MeanBeeAbund) +
                                     (MeanBeeDiversity) +
                                     (MeanBloomAbund) +
                                (FlowerDiversity)+
                                (1|Stand)
                                 )

## **********************************************************
## SEM parasitism
## **********************************************************




bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)

bf.par <- bf(formula.parasite, family="bernoulli")

bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.par +
    set_rescor(FALSE)

prior <- c(set_prior("normal(0,1)", class="b"))

fit <- brm(bform, indiv.data,
           cores=ncores,
           iter = 10^5,
           chains = 4,
           prior=prior,
           control = list(adapt_delta = 0.99))


## **********************************************************
## formulas for the site effects on nest reproductiom
## **********************************************************

ys <- c("SumOffspring",
        "Females")

xvar.NestRepro <- c("(ParasitismRate)",
                    "(BLcover)",
                    "(MeanBloomAbund)",
                    "(FlowerDiversity)",
                    "(MeanBeeAbund)",
                    "(MeanBeeDiversity)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)",
                         sep="+")))
})

bf.offspring <- bf(formulas.NestRepro[[1]], family="poisson")

bform2 <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv +
    bf.offspring +
    set_rescor(FALSE)

fit2 <- brm(bform2, repro.block,
           cores=ncores,
           iter = 10^4,
           chains = 3,
           control = list(adapt_delta = 0.99),
           prior=prior)



## *************************************************************


