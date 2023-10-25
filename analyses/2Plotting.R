## set your working directory to the forestosmia git repo
setwd('forestosmia')
setwd("analyses")
## Script for plotting all of the important explanatory variables.

## Plotting code is based on this tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-average-predictions-with-multilevel-models

rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")

load("../data/indivdata.Rdata")
load("../data/reproblock.Rdata")

indiv.data.orig1 <- indiv.data
indiv.data.orig <- indiv.data
repro.block.orig <- repro.block

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************

indiv.data.orig$Age_LandTrendr <- log(indiv.data.orig$Age_LandTrendr)
labs.age.x <- (pretty(c(indiv.data.orig$Age_LandTrendr),
                      n=8))
axis.age.x <-  standardize.axis(labs.age.x,
                                indiv.data.orig$Age_LandTrendr)
## the max floral abundance makes is a single point, and is makes the
## rest of the points really squished so we dropped it for plotting
no.max.abund <- indiv.data.orig$MeanBloomAbund[indiv.data.orig$MeanBloomAbund !=
                                        max(indiv.data.orig$MeanBloomAbund)]
labs.bloom.abund <- (pretty(c(0, no.max.abund), n=8))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      no.max.abund)
indiv.data.orig$Hectares <- log(indiv.data.orig$Hectares)
labs.hectares <- (pretty(indiv.data.orig$Hectares, n=8))
axis.hectares <-  standardize.axis(labs.hectares,
                                indiv.data.orig$Hectares)
labs.bee.abund2 <- (pretty(c(0,
                             indiv.data.orig$MeanBeeAbund), n=8))
axis.bee.abund2 <-  standardize.axis(labs.bee.abund2,
                                     indiv.data.orig$MeanBeeAbund)
labs.flower.div <- (pretty(indiv.data.orig$FlowerDiversity, n=8))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     indiv.data.orig$FlowerDiversity)
labs.flower.div.repro <- (pretty(repro.block.orig$FlowerDiversity, n=8))
axis.flower.div.repro <-  standardize.axis(labs.flower.div.repro,
                                           repro.block.orig$FlowerDiversity)

## the max bee abundance makes is a single point, and is makes the
## rest of the points really squished so we dropped it for plotting
no.max.abund <- indiv.data.orig$MeanBeeAbund[indiv.data.orig$MeanBeeAbund !=
                                        max(indiv.data.orig$MeanBeeAbund)]
labs.bee.abund <- (pretty(c(0, no.max.abund), n=8))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                      no.max.abund)

labs.bee.div <- (pretty(c(0, indiv.data.orig$MeanBeeDiversity), n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  indiv.data.orig$MeanBeeDiversity)

## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************

makeBeeDivAbundFigs <- function(yvar, ylabs, yrate,
                                abund.sig="gt95",
                                div.sig="gt95",
                                bloomabund.sig="gt95",
                                tags=c("A", "B", "C")){
  ## function for plotting model results, curve color is determined by
  ## level of support imputed via the XX.sig arguments. 
  
  chooseColors <- function(sig){
    if(sig == "gt95"){
      cols <-  "Blues"
    } else if(sig == "gt90"){
      cols <- "Oranges"
    }else{
      cols <-  "Greys"
    }
    return(cols)
  }

  col.abund <- chooseColors(abund.sig)
  col.div <- chooseColors(div.sig)
  col.bloomabund <- chooseColors(bloomabund.sig)

    ## parasitism ~ bee diversity
  newdata.beediv <- crossing(MeanBeeDiversity =
                               seq(min(data.par$MeanBeeDiversity),
                                   max(data.par$MeanBeeDiversity),
                                   length.out=10),
                             Age_LandTrendr =mean(data.par$Age_LandTrendr),
                             Elev=mean(data.par$Elev),
                             Owner="OwnerC",
                             Stand="Backgrove",
                             StandBlock= "BackgroveA",
                             MeanBloomAbund=mean(data.par$MeanBloomAbund),
                             Hectares=mean(data.par$Hectares),
                             FlowerDiversity=mean(data.par$FlowerDiversity),
                             MeanBeeAbund=mean(data.par$MeanBeeAbund)
                             )

  ## predict values based on generated data and model parameters
  pred_beediv <- fit.parasite %>%
    epred_draws(newdata = newdata.beediv,
                resp = yvar)

  ## to see range of predicted values
  pred_beediv %>%
    group_by(MeanBeeDiversity) %>%
    summarise(mean(.epred))
  
  p1.parasite <- ggplot(pred_beediv, aes(x = MeanBeeDiversity, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = col.div) +
    labs(x = "Bee community diversity", y = ylabs, tag=tags[2],
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    ylim(0,1) +
    scale_x_continuous(
      breaks = axis.bee.div,
      labels =  labs.bee.div) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_ms() +
    ## theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=!! sym(yrate), x=MeanBeeDiversity),
               color="black", cex=2)

  ## parasitism ~ bee abundance
  newdata.beeabund <- crossing(MeanBeeAbund =
                                 seq(min(data.par$MeanBeeAbund),
                                     max(data.par$MeanBeeAbund),
                                     length.out=10),
                               Age_LandTrendr =mean(data.par$Age_LandTrendr),
                               Elev=mean(data.par$Elev),
                               Owner="OwnerC",
                               Stand="Backgrove",
                               StandBlock= "BackgroveA",
                               MeanBloomAbund=mean(data.par$MeanBloomAbund),
                               Hectares=mean(data.par$Hectares),
                               FlowerDiversity=mean(data.par$FlowerDiversity),
                               MeanBeeDiversity=mean(data.par$MeanBeeDiversity)
                               )

  pred_beeabund <- fit.parasite %>%
    epred_draws(newdata = newdata.beeabund ,
                resp = yvar)

  pred_beeabund %>%
    group_by(MeanBeeAbund) %>%
    summarise(mean(.epred))

  p2.parasite <- ggplot(pred_beeabund, aes(x = MeanBeeAbund, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = col.abund) +
    labs(x = "Bee abundance", y = ylabs, tag=tags[1],
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    ylim(0,1) +
    scale_x_continuous(
      breaks = axis.bee.abund,
      labels =  labs.bee.abund) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_ms() +
    ## theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=!! sym(yrate), x=MeanBeeAbund),
               color="black", cex=2)


  ## parasitism ~ floral abundance
  newdata.bloomabund <- crossing(MeanBloomAbund =
                                   seq(min(data.par$MeanBloomAbund),
                                       max(data.par$MeanBloomAbund),
                                       length.out=10),
                                 Age_LandTrendr =mean(data.par$Age_LandTrendr),
                                 Elev=mean(data.par$Elev),
                                 Owner="OwnerC",
                                 Stand="Backgrove",
                                 StandBlock= "BackgroveA",
                                 MeanBeeAbund=mean(data.par$MeanBeeAbund),
                                 Hectares=mean(data.par$Hectares),
                                 FlowerDiversity=mean(data.par$FlowerDiversity),
                                 MeanBeeDiversity=mean(data.par$MeanBeeDiversity)
                                 )

  pred_bloomabund <- fit.parasite %>%
    epred_draws(newdata = newdata.bloomabund ,
                resp = yvar)

  pred_bloomabund %>%
    group_by(MeanBloomAbund) %>%
    summarise(mean(.epred))

  p3.parasite <- ggplot(pred_bloomabund, aes(x = MeanBloomAbund, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = col.bloomabund) +
    labs(x = "Floral abundance", y = ylabs, tag=tags[3],
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    ylim(0,1) +
    scale_x_continuous(
      breaks = axis.bloom.abund,
      labels =  labs.bloom.abund) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_ms() +
    ## theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=!! sym(yrate), x=MeanBloomAbund),
               color="black", cex=2)

  
  ggsave(p1.parasite, file=sprintf("figures/%s_beeDiv.pdf", yvar),
         height=4, width=7.5)
  ggsave(p2.parasite, file=sprintf("figures/%s_beeAbund.pdf", yvar),
         height=4, width=7.5)
  ggsave(p3.parasite, file=sprintf("figures/%s_bloomAbund.pdf", yvar),
         height=4, width=7.5)

  parasite.all <- grid.arrange(p2.parasite, p1.parasite, p3.parasite, ncol=3)

  ggsave(parasite.all, file=sprintf("figures/all_%s.pdf", yvar),
         height=4, width=12)

}


load('saved/parasiteFit_osmia_all_data_CrithidiaApicystisAscophaera.Rdata')
round(r2, 2) 
data.par <- spec.data[spec.data$WeightsPar == 1, ]
indiv.data <- spec.data
makeBeeDivAbundFigs("Apicystis", "Apicystis Prevalence",
                    "ApicysRate",
                    abund.sig="gt95",
                    div.sig=FALSE,
                    bloomabund.sig="gt90",
                    tags=c("G", "H", "I"))

makeBeeDivAbundFigs("Ascophaera", "Ascophaera Prevalence",
                    "AscoRate",
                    abund.sig="gt95",
                    div.sig="gt95",
                    bloomabund.sig=FALSE,
                    tags=c("D", "E", "F"))

makeBeeDivAbundFigs("Crithidia", "Crithidia Prevalence",
                    "CrithRate",
                      abund.sig="gt90",
                    div.sig="gt90",
                    bloomabund.sig=FALSE,
                    tags=c("J", "K", "L"))


load('saved/parasiteFit_osmia_all_data_AnyParasite.Rdata')
round(r2, 2)
makeBeeDivAbundFigs("AnyParasite", "Any Parasite Prevalence",
                    "ParasitismRate",
                    abund.sig="gt95",
                    div.sig="gt95",
                    bloomabund.sig=FALSE)


## ***********************************************************************
## bee offspring
## ***********************************************************************

load(file="saved/offspringFitMod.Rdata")

newdata.floraldiv <- crossing(FlowerDiversity =
                                  seq(min(repro.block$FlowerDiversity),
                                      max(repro.block$FlowerDiversity),
                                      length.out=2),
                              Age_LandTrendr =mean(repro.block$Age_LandTrendr),
                              Elev=mean(repro.block$Elev),
                              Owner="OwnerC",
                              Stand="Backgrove",
                              MeanBloomAbund=mean(repro.block$MeanBloomAbund),
                              MeanBeeAbund=mean(repro.block$MeanBeeAbund),
                              MeanBeeDiversity=mean(repro.block$MeanBeeDiversity),
                              ParasitismRate=mean(repro.block$ParasitismRate, na.rm=TRUE)
                              )

pred_floraldiv <- fit2 %>%
    epred_draws(newdata = newdata.floraldiv ,
                resp = "SumOffspring")

pred_floraldiv %>%
    group_by(FlowerDiversity) %>%
     summarise(mean(.epred))

p1.offspring <- ggplot(pred_floraldiv, aes(x = FlowerDiversity, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Floral diversity", y = "Offspring",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.flower.div.repro,
        labels =  labs.flower.div.repro) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=repro.block,
               aes(y=SumOffspring, x=FlowerDiversity), color="grey40",
               cex=2)

ggsave(p1.offspring, file="figures/offspring_floraldiv.pdf",
       height=4, width=5)

## ***********************************************************************
## bee community- floral div, hectares
## ***********************************************************************

## bee abundance ~ hectares

indiv.data <- indiv.data[indiv.data$Weights == 1,]
fit <- fit.parasite

newdata.hectares <- crossing(Hectares =
                              seq(min(indiv.data$Hectares),
                                  max(indiv.data$Hectares),
                                  length.out=10),
                          Age_LandTrendr =mean(indiv.data$Age_LandTrendr),
                          Elev=mean(indiv.data$Elev),
                          Owner="ODF",
                          Stand="Backgrove",
                          MeanBloomAbund=mean(indiv.data$MeanBloomAbund),
                          FlowerDiversity=mean(indiv.data$FlowerDiversity),
                          MeanBeeAbund=mean(indiv.data$MeanBeeAbund),
                          MeanBeeDiversity=mean(indiv.data$MeanBeeDiversity)
                          )

pred_hectares <- fit %>%
    epred_draws(newdata = newdata.hectares ,
                resp = "MeanBeeAbund")

p1.bee <- ggplot(pred_hectares, aes(x = Hectares, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Oranges") +
    labs(x = "Stand area (log hectares)", y = "Bee abundance", tag="A",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.hectares,
        labels =  labs.hectares) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund) +
    coord_cartesian(ylim = range(axis.bee.abund)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$MeanBeeAbund != max(indiv.data$MeanBeeAbund),],
               aes(y=MeanBeeAbund, x=Hectares), color="grey40",
               cex=2)

## Bee diversity ~ floral diversity

newdata.floraldiv <- crossing(FlowerDiversity =
                                  seq(min(indiv.data$FlowerDiversity),
                                      max(indiv.data$FlowerDiversity),
                                      length.out=10),
                              Age_LandTrendr =mean(indiv.data$Age_LandTrendr),
                              Elev=mean(indiv.data$Elev),
                              Owner="ODF",
                              Stand="Backgrove",
                              MeanBloomAbund=mean(indiv.data$MeanBloomAbund),
                              Hectares=mean(indiv.data$Hectares),
                              MeanBeeAbund=mean(indiv.data$MeanBeeAbund)
                              )

pred_floraldiv <- fit %>%
    epred_draws(newdata = newdata.floraldiv ,
                resp = "MeanBeeDiversity")

p2.bee <- ggplot(pred_floraldiv, aes(x = FlowerDiversity, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Floral diversity", y = "Bee diversity", tag="B",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    coord_cartesian(ylim = range(axis.bee.div)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data,
               aes(y=MeanBeeDiversity, x=FlowerDiversity), color="grey40",
               cex=2) 
  

bee.plots <- grid.arrange(p1.bee, p2.bee, ncol=2)

ggsave(p1.bee, file="figures/bee_hectares.pdf",
       height=4, width=5)

ggsave(p2.bee, file="figures/bee_floralDiv.pdf",
       height=4, width=5)

ggsave(bee.plots, file="figures/beeComm.pdf",
       height=4, width=8)

## ***********************************************************************
## bee community- age
## ***********************************************************************
## floral abundance ~ year post harvest

newdata.age <- crossing(Age_LandTrendr=
                            seq(min(indiv.data$Age_LandTrendr),
                                max(indiv.data$Age_LandTrendr),
                                length.out=10),
                        Hectares =mean(indiv.data$Hectares),
                        Elev=mean(indiv.data$Elev),
                        Owner="ODF",
                        Stand="Backgrove",
                        MeanBloomAbund=mean(indiv.data$MeanBloomAbund),
                        FlowerDiversity=mean(indiv.data$FlowerDiversity),
                        MeanBeeAbund=mean(indiv.data$MeanBeeAbund),
                        MeanBeeDiversity=mean(indiv.data$MeanBeeDiversity)
                        )

## flower diversity ~ years post harvest
pred_fdiv_age <- fit %>%
    epred_draws(newdata = newdata.age,
                resp = "FlowerDiversity")

p1.flower.age <- ggplot(pred_fdiv_age, aes(x = Age_LandTrendr, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Year post-harvest (log)", y = "Floral diversity", tag="A",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.age.x,
        labels =  labs.age.x) +
    scale_y_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    coord_cartesian(ylim = range(axis.flower.div)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data,
               aes(y=FlowerDiversity, x=Age_LandTrendr), color="grey40",
               cex=2)


## flower abundance ~ years post harvest
pred_fabund_age <- fit %>%
    epred_draws(newdata = newdata.age,
                resp = "MeanBloomAbund")

p2.flower.age <- ggplot(pred_fabund_age, aes(x = Age_LandTrendr, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Year post-harvest (log)", y = "Floral abundance", tag="B",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.age.x,
        labels =  labs.age.x) +
    scale_y_continuous(
        breaks = axis.bloom.abund,
        labels =  labs.bloom.abund) +
    coord_cartesian(ylim = range(axis.bloom.abund)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data,
               aes(y=MeanBloomAbund, x=Age_LandTrendr), color="grey40",
               cex=2)

## bee diversity ~ years post harvest
pred_bdiv_age <- fit %>%
    epred_draws(newdata = newdata.age,
                resp = "MeanBeeDiversity")

p3.bee.age <- ggplot(pred_bdiv_age, aes(x = Age_LandTrendr, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Year post-harvest (log)", y = "Bee diversity", tag="C",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.age.x,
        labels =  labs.age.x) +
    scale_y_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    coord_cartesian(ylim = range(axis.bee.div)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data,
               aes(y=MeanBeeDiversity, x=Age_LandTrendr), color="grey40",
               cex=2)

## bee abundance ~ years post harvest

pred_babund_age <- fit %>%
    epred_draws(newdata = newdata.age,
                resp = "MeanBeeAbund")

p4.bee.age <- ggplot(pred_babund_age, aes(x = Age_LandTrendr, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Year post-harvest (log)", y = "Bee abundance", tag="D",
         fill = "Credible interval") +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.age.x,
        labels =  labs.age.x) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund) +
    coord_cartesian(ylim = range(axis.bee.abund)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
   theme_ms() +
 ##     theme_dark_black()+
    geom_point(data=indiv.data,
               aes(y=MeanBeeAbund, x=Age_LandTrendr), color="grey40",
               cex=2)

age.plots <- grid.arrange(p1.flower.age, p2.flower.age, p3.bee.age,
                          p4.bee.age,
                          ncol=2)

ggsave(age.plots, file="figures/standage.pdf",
       height=7, width=8)

## ***********************************************************************
## Landowner
## ***********************************************************************
load("../data/sitedata.Rdata")
site.data$AgeBin <- NA
site.data$AgeBin[site.data$Age_LandTrendr  <= 3]  <- "1-3"
site.data$AgeBin[site.data$Age_LandTrendr  > 3]  <- "4-6"
site.data$AgeBin[site.data$Age_LandTrendr  > 6]  <- "7-9"
site.data$AgeBin[site.data$Age_LandTrendr  >= 10]  <- "10+"

site.data$AgeBin <- factor(site.data$AgeBin,
                           levels=c("1-3", "4-6", "7-9", "10+"))

p1.owner <- ggplot(site.data , aes(x=AgeBin, y=FlowerDiversity)) +
    geom_bar(aes(fill = Owner),
             stat = "identity", position = position_dodge(0.8),
             width = 0.7) +
    coord_flip() +
    scale_fill_brewer(palette="Blues") +
   theme_ms() +
 ##     theme_dark_black()+
    theme(legend.position="none") +
    ylab("Floral community diversity") +
    xlab("")


ggsave(p1.owner, file="figures/flowerDiv_owner.pdf",
       height=4, width=5)


p2.owner <- ggplot(site.data , aes(x=AgeBin, y=MeanBloomAbund)) +
    geom_bar(aes(fill = Owner),
             stat = "identity", position = position_dodge(0.8),
             width = 0.7) +
    coord_flip() +
    scale_fill_brewer(palette="Blues") +
   theme_ms() +
 ##     theme_dark_black()+
    ylab("Floral abundance") +
    xlab("")


owner.all <- grid.arrange(p1.owner, p2.owner, ncol=2)

ggsave(owner.all, file="figures/all_owner.pdf",
       height=4, width=10)


