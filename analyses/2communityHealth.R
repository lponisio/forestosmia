setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(car)
library(ggplot2)
library(viridis)

load("../data/indivdata.Rdata")
load("../data/sitedata.Rdata")
load("../data/reproblock.Rdata")

## **********************************************************
## formula for site effects on the bee community
## **********************************************************
## tree richness, canopy and DBH are colinear

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
formula.flower.div <- formula(FlowerDiversity ~
                                  (BLcover) +
                                  (AgePoly1) +
                                  (AgePoly2) +
                                  (AgePoly3) +
                                  (AgePoly4)
                              )
## flower abund
formula.flower.abund <- formula(MeanBloomAbund ~
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
formula.bee.div <- formula(MeanBeeDiversity ~
                               (BLcover) +
                               (MeanBloomAbund) +
                               ## (MeanBeeAbund) +
                               (FlowerDiversity) +
                               (AgePoly1) +
                               (AgePoly2)
                               ## (AgePoly3)
                           )

## bee abund
formula.bee.abund <- formula(MeanBeeAbund ~
                                 (BLcover) +
                                 (MeanBloomAbund) +
                                 (FlowerDiversity) +
                                 (AgePoly1) +
                                 (AgePoly2)
                                 ## (AgePoly3)
                             )

## **********************************************************
## SEM bees and flowers
## **********************************************************

bees.flowers= psem(
    plants = lm(formula.flower.abund,
                      data = site.data),
    plants2 = lm(formula.flower.div,
                       data = site.data),
    bees = lm(formula.bee.abund,
                    data = site.data),
    bees2 = lm(formula.bee.div,
                     data = site.data)
)

summary(bees.flowers)
dTable <- dSep(bees.flowers)
dTable
fisherC(dTable)

## **********************************************************
## formula for bee community effects on parasitism
## **********************************************************

## Apicystis Ascophaera Crithidia AnyParasite

## bee abund and bloom abund colinear

site.data.par <- site.data[!is.na(site.data$InfectedIndividuals),]
## cannot use cbind in peicewise SEM

formula.parasite.site <- formula(AnyParasite ~
                                     (MeanBeeAbund) +
                                     (MeanBeeDiversity) +
                                     (MeanBloomAbund) +
                                     (FlowerDiversity) +
                                     (1|Stand)
                                 )


## **********************************************************
## SEM parasitism
## **********************************************************
indiv.data$Stand <- as.factor(indiv.data$Stand)

## trying with the quasi poison and offspet fro total screened,
## dispersion permaters looks silly
bees.flowers.par= psem(
    plants = lm(formula.flower.abund,
                data = site.data),
    plants2 = lm(formula.flower.div,
                 data = site.data),
    bees = lm(formula.bee.abund,
              data = site.data),
    bees2 = lm(formula.bee.div,
               data = site.data),
    par = glm(formula=formula.parasite.site,
              data = site.data.par,
              offset=site.data.par$TestedTotals,
              family="quasipoisson")
)

summary(glm(formula=formula.parasite.site,
                       data = site.data.par,
                       offset=site.data.par$TestedTotals,
                       family="quasipoisson"))


##
bees.flowers.par= psem(
    plants = lm(formula.flower.abund,
                data = site.data),
    plants2 = lm(formula.flower.div,
                 data = site.data),
    bees = lm(formula.bee.abund,
              data = site.data),
    bees2 = lm(formula.bee.div,
               data = site.data),
    par = glmer(formula=formula.parasite.site,
              data = indiv.data,
              family="binomial",
              glmerControl(optimizer="bobyqa"))
)

summary(bees.flowers.par)
dTable <- dSep(bees.flowers.par)
dTable
fisherC(dTable)

summary(glmer(formula=formula.parasite.site,
              data = indiv.data,
              family="binomial"))


## try to be at the site level with infected, total tested
## bees.flowers.par= psem(
##                   plants = lm(formula.flower.abund,
##                               data = site.data.par),
##                   plants2 = lm(formula.flower.div,
##                                data = site.data.par),
##                   bees = lm(formula.bee.abund,
##                             data = site.data.par),
##                   bees2 = lm(formula.bee.div,
##                              data = site.data.par),
##                   par = do.call(glm,
##                                 list(formula=formula.parasite.site,

##                                      data = site.data.par,
##                                      weights=site.data.par$TestedTotals,
##                                      family="binomial"))
##               )


## results are different than above, not as many stands that have
## parasite data

## **********************************************************
## formulas for the site effects on nest reproductiom
## **********************************************************
repro.block$Block <- as.factor(repro.block$Block)
repro.block$Stand <- as.factor(repro.block$Stand)

ys <- c("SumOffspring", "Females")

xvar.NestRepro <- c("(ParasitismRate)",
                    "(BLcover)",
                    "(MeanBloomAbund)",
                    "(FlowerDiversity)",
                    "(MeanBeeAbund)",
                    "(MeanBeeDiversity)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)", ## "(1|Block)",
                         sep="+")))
})

mod <- lmer(formulas.NestRepro[[1]], data=repro.block)

## it's count data, but numbers are big and few 0, resid look perfect!
plot(density(residuals(mod)))

flowers.repro= psem(
    plants = lm(formula.flower.abund,
                data = site.data),
    plants2 = lm(formula.flower.div,
                 data = site.data),
    bees = lm(formula.bee.abund,
              data = site.data),
    bees2 = lm(formula.bee.div,
               data = site.data),
    nest= lmer(formulas.NestRepro[[1]],
                data = repro.block)
)

summary(flowers.repro)

dTable <- dSep(flowers.repro)
dTable
fisherC(dTable)


## same results

flowers.repro.f= psem(
   plants = lm(formula.flower.abund,
                data = site.data),
    plants2 = lm(formula.flower.div,
                 data = site.data),
    bees = lm(formula.bee.abund,
              data = site.data),
    bees2 = lm(formula.bee.div,
               data = site.data),
    nest= lmer(formulas.NestRepro[[2]],
                data = repro.block)
)

summary(flowers.repro.f)

dTable <- dSep(flowers.repro.f)
dTable
fisherC(dTable)



## *************************************************************


