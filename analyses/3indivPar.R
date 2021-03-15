## setwd("~/Dropbox/forestosmia")
setwd("analyses")
rm(list=ls())
library(piecewiseSEM)
library(lme4)
library(lmerTest)

load("../data/parasite.Rdata")
load("../data/NestRepro.Rdata")

## *************************************************************
## make formulas for path analyses
## *************************************************************

## formula for site effects on the bee community
formula.bee <- formula(AnyParasite~ standintensity +
                         MeanDBH + TreeRichness  + FlowerRichness)

## formulas for the site effects on parasites and pathogens

infections <- c("Apicystis", "Crithidia", "Ascophaera")

xvar.par.path <- c("standintensity",
                    "MeanDBH",
                    "TreeRichness",
                    "FlowerRichness")

formulas.par <-lapply(infections, function(x) {
    as.formula(paste(x, "~",
                     paste(xvar.par.path, collapse="+")))
})


## *************************************************************

calcMods <- function(this.formula,
                     formula.bee,
                     col.trials,
                     dats,
                     site.char){
    colnames(dats)[colnames(dats) == col.trials] <- "trials"
    mod = psem(
        BeeDensity = do.call(lm,
                             list(formula=formula.bee,
                                  data=site.char)),
        ParPath = do.call(glm,
                          list(formula=this.formula,
                               data = dats,
                               weights=dats$trials,
                               family="binomial"))
    )
    print(summary(mod))
    return(mod)
}

## *************************************************************
## osmia
## *************************************************************

osmia.mods.par <- lapply(formulas.par, calcMods,
                          formula.bee=formula.bee,
                          col.trials= "ScreenedPar",
                          dats=bombus,
                          site.char=site.char)

names(osmia.mods.par) <- parasites


print("bombus parasites")
lapply(bombus.mods.par, summary)
lapply(bombus.mods.par, rsquared)



