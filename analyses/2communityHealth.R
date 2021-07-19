## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(car)

load("../data/indivdata.Rdata")
load("../data/sitedata.Rdata")
load("../data/reproblock.Rdata")


## *************************************************************
## make formulas for path analyses
## *************************************************************


## **********************************************************
## formula for site effects on the bee community
## **********************************************************

site.data <- site.data[!is.na(site.data$InfectedIndividuals),]


formula.bee <- formula(MeanBeeRichness ~
                         scale(BLcover) +
                         scale(MeanCanopy) + scale(TreeRichness) +  scale(MeanDBH)
                       + scale(MeanBloomAbund))
#+ scale(FlowerRichness) 

summary(lm(formula.bee,
           data = site.data))

vif(lm(formula.bee,
       data = site.data))


## *****

formula.bee <- formula(MeanBeeAbund ~
                         scale(BLcover) +
                         scale(MeanCanopy) + scale(TreeRichness) + scale(MeanDBH) + 
                         #scale(FlowerRichness) + 
                         scale(MeanBloomAbund))

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

## tried above with InfectedCrith, InfectedAsco, InfectedApicys

## *****try it out with other variables

formula.parasite <- formula(cbind(InfectedCrith,
                                  TestedTotals)~
                              scale(MeanBeeRichness)
                            + scale(TreeRichness) 
                            + scale(BLcover)
                            + scale(MeanCanopy)
                            + scale(MeanBloomAbund) 
                            + scale(MeanDBH))
                           # + scale(FlowerRichness))

parasite.mod <-  glm(formula.parasite,
                     data = site.data,
                     family="binomial")
summary(parasite.mod)


vif(parasite.mod)

## **********************************************************
## formulas for the site effects on nest reproductiom
## **********************************************************

ys <- c("FM_ratio", "SumOffspring", "Females")

#xvar.NestRepro <- c("scale(ParasitismRate)")

xvar.NestRepro <- c("scale(ParasitismRate)",
                    "scale(MeanBeeRichness)",
                    "scale(BLcover)",
                    "scale(TreeRichness)",
                    "scale(MeanCanopy)",
                    "scale(MeanBloomAbund)", 
                    "scale(MeanDBH)")

                    #"scale(FlowerRichness)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)", "(1|Block)",
                         sep="+")))
})


## offspring
repro.mod <- glmer(formulas.NestRepro[[2]],
                   family="poisson",
                   data = repro.block)

summary(repro.mod)
## plot(density(residuals(repro.mod)))

## *****
## fm ratio
fm.mod <- lmer(formulas.NestRepro[[1]],
               data = repro.block)

summary(fm.mod)

vif(fm.mod)


## *****
## females
f.mod <- glmer(formulas.NestRepro[[3]],
               data = repro.block,
               family="poisson")

summary(f.mod)

vif(f.mod)

#********with bee abund

ys <- c("FM_ratio", "SumOffspring", "Females")

#xvar.NestRepro <- c("scale(ApicysRate)")

xvar.NestRepro <- c("scale(ParasitismRate)",
                    "scale(MeanBeeAbund)",
                    "scale(BLcover)",
                    "scale(TreeRichness)",
                    "scale(MeanCanopy)",
                    "scale(MeanBloomAbund)", 
                    "scale(MeanDBH)")

formulas.NestRepro <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar.NestRepro, collapse="+"),
                         "(1|Stand)", "(1|Block)",
                         sep="+")))
})


## offspring
repro.mod <- glmer(formulas.NestRepro[[2]],
                   family="poisson",
                   data = repro.block)

summary(repro.mod)
## plot(density(residuals(repro.mod)))

## *****
## fm ratio
fm.mod <- lmer(formulas.NestRepro[[1]],
               data = repro.block)

summary(fm.mod)

## *****
## females
f.mod <- glmer(formulas.NestRepro[[3]],
               data = repro.block,
               family="poisson")

summary(f.mod)


vif(f.mod)
## *************************************************************

mod.offspring  <-  psem(
  Bee = lm(formula.bee,
           data = site.data),
  Parasite = glm(formula.parasite,
                 data = site.data,
                 family="binomial"),
  NestRepro = lmer(formulas.NestRepro[[2]],
                   data = repro.block))


summary(mod.offspring)

lapply(osmia.mods, rsquared)



## *************************************************************
## sanity check using a linear models of Any Parasite
## *************************************************************


bee.par.mod <- lm(AnyParasite ~ BLcover + MeanDBH +
                    TreeRichness + FlowerRichness,
                  data=site.data)
AIC(bee.par.mod)
vif(bee.par.mod)
summary(bee.par.mod)
plot(density(bee.par.mod$resid))
