## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)

load("../data/parasite.Rdata")
load("../data/NestRepro.Rdata")


## *************************************************************
## make formulas for path analyses
## *************************************************************

## formula for site effects on the bee community

formula.bee <- formula(AnyParasite~ standintensity +
                                    MeanDBH + TreeRichness  + FlowerRichness)

## formulas for the site effects on nest reproductiom

ys <- c("FM_ratio", "SumOffspring")
xvar.NestRepro <- c("AnyParasite",
                   "standintensity",
                   "MeanDBH",
                   "TreeRichness",
                   "FlowerRichness")

formulas.NestRepro <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.NestRepro, collapse="+"),
                           "(1|Stand)",
                           sep="+")))
})

## *************************************************************

calcMods <- function(this.formula, formula.bee,
                     dats,
                     site.char){
    mod = psem(
        Parasite = lm(formula.bee,
                        data = parasite),
        NestRepro = lmer(this.formula,
                       data = dats))
}
## *************************************************************
osmia <- parasite

## osmia
osmia.mods <- lapply(formulas.NestRepro, calcMods,
                      formula.bee, osmia, parasite)


print("osmia")
lapply(osmia.mods, summary)
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
