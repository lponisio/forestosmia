## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(car)

load("../data/indivData.Rdata")
load("../data/siteData.Rdata")
load("../data/NestRepro.Rdata")


## *************************************************************
## make formulas for path analyses
## *************************************************************


## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## site.data <- site.data[!is.na(site.data$InfectedIndividuals),]

#formula.bee <- formula(MeanBeeRichness ~
                           #scale(standintensity) +
                           #scale(MeanCanopy) + scale(TreeRichness)  +
                           #scale(MeanBloomAbund) +    scale(MeanDBH) +
                           #scale(FlowerRichness) + scale(Acres))

formula.bee <- formula(MeanBeeRichness ~
                         scale(standintensity) +
                         scale(MeanCanopy) + scale(TreeRichness) +    scale(MeanDBH) + 
                         scale(FlowerRichness) + scale(Acres))
summary(lm(formula.bee,
           data = site.data))

vif(lm(formula.bee,
       data = site.data))


## *****

#formula.bee <- formula(MeanBeeAbund ~
                    #scale(standintensity) +
                    #scale(MeanCanopy)    +     scale(TreeRichness)  +
                    #scale(MeanBloomAbund) +    scale(MeanDBH) +
                    #scale(FlowerRichness) +    scale(Acres))

formula.bee <- formula(MeanBeeAbund ~
                         scale(standintensity) +
                         scale(MeanCanopy) + scale(TreeRichness) + scale(MeanDBH) + 
                         scale(FlowerRichness) + scale(Acres))

summary(lm(formula.bee,
           data = site.data))

## **********************************************************
## formula for bee community effects on parasitism
## **********************************************************

formula.parasite <- formula(cbind(InfectedIndividuals,
                                  TestedTotals)~
                                scale(MeanBeeRichness))

parasite.mod <-  glm(formula.parasite,
                     data = site.data,
                     family="binomial")
summary(parasite.mod)

## **** with abundance

formula.parasite <- formula(cbind(InfectedIndividuals,
                                  TestedTotals)~
                              scale(MeanBeeAbund))

parasite.mod <-  glm(formula.parasite,
                     data = site.data,
                     family="binomial")
summary(parasite.mod)

## *****try it out with other variables

formula.parasite <- formula(cbind(InfectedIndividuals,
                                  TestedTotals)~
                              scale(MeanBeeAbund)
                              + scale(TreeRichness) + scale(standintensity)
                              + scale(MeanBloomAbund) + scale(MeanDBH)
                              + scale(FlowerRichness) + scale(Acres))


vif(parasite.mod)

## **********************************************************
## formulas for the site effects on nest reproductiom
## **********************************************************

ys <- c("FM_ratio", "SumOffspring", "Females")

#xvar.NestRepro <- c("scale(CrithRate)")

xvar.NestRepro <- c("scale(ParasitismRate)",
      "scale(MeanBeeAbund)",
      "scale(standintensity)",
      "scale(TreeRichness)",
      "scale(MeanBloomAbund)", "scale(MeanDBH)",
      "scale(FlowerRichness)", "scale(Acres)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)", "(1|Block)",
                         sep="+")))
})


## offspring
repro.mod <- glmer(formulas.NestRepro[[2]],
                   family="poisson",
                   data = repro.nest)

summary(repro.mod)
## plot(density(residuals(repro.mod)))

## *****
## fm ratio
fm.mod <- lmer(formulas.NestRepro[[1]],
                data = repro.nest)

summary(fm.mod)

## *****
## females
f.mod <- glmer(formulas.NestRepro[[3]],
               data = repro.nest,
               family="poisson")

summary(f.mod)

## *************************************************************

mod.offspring  <-  psem(
    Bee = lm(formula.bee,
                   data = site.data),
    Parasite = glm(formula.parasite,
                   data = site.data,
                   family="binomial"),
    NestRepro = lmer(formulas.NestRepro[[2]],
                      data = repro.nest))


summary(mod.offspring)

lapply(osmia.mods, rsquared)



## *************************************************************
## sanity check using a linear models of Any Parasite
## *************************************************************

library(car)

bee.par.mod <- lm(AnyParasite ~ standintensity + MeanDBH +
                     TreeRichness + FlowerRichness,
                       data=parasite)
AIC(bee.par.mod)
vif(bee.par.mod)
summary(bee.par.mod)
plot(density(bee.par.mod$resid))
