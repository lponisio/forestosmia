setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')

## setwd('~/Dropbox (University of Oregon)/forestosmia')
## Script for plotting all of the important explanatory variables.

setwd("analyses")
rm(list=ls())

source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")

load("../data/indivdata.Rdata")
load("../data/reproblock.Rdata")

indiv.data.orig1 <- indiv.data
indiv.data.orig <- indiv.data
repro.block.orig <- repro.block

load(file="saved/offspringFitMod.Rdata")
load(file="saved/parasiteFitMod.Rdata")

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************

indiv.data.orig$Age_LandTrendr <- log(indiv.data.orig$Age_LandTrendr)

labs.age.x <- (pretty(c(indiv.data.orig$Age_LandTrendr),
                      n=10))

axis.age.x <-  standardize.axis(labs.age.x, indiv.data.orig$Age_LandTrendr)


labs.bloom.abund <- (pretty(c(0, indiv.data.orig$MeanBloomAbund), n=10))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      indiv.data.orig$MeanBloomAbund)


indiv.data.orig$Acres <- log(indiv.data.orig$Acres)

labs.acres <- (pretty(indiv.data.orig$Acres, n=10))
axis.acres <-  standardize.axis(labs.acres,
                                indiv.data.orig$Acres)


labs.bee.abund2 <- (pretty(c(0,
                             indiv.data.orig$MeanBeeAbund[
                                                 !is.na(indiv.data.orig$AnyParasite)]),
                           n=10))
axis.bee.abund2 <-  standardize.axis(labs.bee.abund2,
                                     indiv.data.orig$MeanBeeAbund)


labs.flower.div <- (pretty(indiv.data.orig$FlowerDiversity, n=10))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     indiv.data.orig$FlowerDiversity)


labs.flower.div.repro <- (pretty(repro.block.orig$FlowerDiversity, n=10))
axis.flower.div.repro <-  standardize.axis(labs.flower.div.repro,
                                           repro.block.orig$FlowerDiversity)

labs.bee.abund <- (pretty(c(0, indiv.data.orig$MeanBeeAbund), n=10))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                    indiv.data.orig$MeanBeeAbund)

labs.bee.div <- (pretty(c(0, indiv.data.orig$MeanBeeDiversity), n=10))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  indiv.data.orig$MeanBeeDiversity)

## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************

data.par <- indiv.data[indiv.data$WeightsPar == 1, ]

p1.parasite  <- fit %>%
    spread_draws(
        b_MeanBloomAbund_Intercept,
        b_FlowerDiversity_Intercept,
        b_MeanBeeAbund_Intercept,
        b_MeanBeeDiversity_Intercept,
        b_AnyParasite_Intercept,
        b_MeanBloomAbund_Age_LandTrendr,
        b_MeanBloomAbund_Elev,
        b_MeanBloomAbund_OwnerOwnerB,
        b_MeanBloomAbund_OwnerOwnerC,
        b_MeanBloomAbund_OwnerODF,
        b_FlowerDiversity_Age_LandTrendr,
        b_FlowerDiversity_IAge_LandTrendrE2,
        b_FlowerDiversity_Elev,
        b_FlowerDiversity_OwnerOwnerB,
        b_FlowerDiversity_OwnerOwnerC,
        b_FlowerDiversity_OwnerODF,
        b_MeanBeeAbund_MeanBloomAbund,
        b_MeanBeeAbund_Age_LandTrendr,
        b_MeanBeeAbund_Acres,
        b_MeanBeeDiversity_FlowerDiversity,
        b_MeanBeeDiversity_Age_LandTrendr,
        b_MeanBeeDiversity_Acres,
        b_AnyParasite_MeanBeeAbund,
        b_AnyParasite_FlowerDiversity,
        b_AnyParasite_MeanBeeDiversity,
        b_AnyParasite_MeanBloomAbund,
        `b_MeanBeeAbund_Age_LandTrendr:Acres`,
        `b_MeanBeeDiversity_Age_LandTrendr:Acres`
    ) %>%
    mutate(MeanBeeDiversity =
               list(seq(min(data.par$MeanBeeDiversity),
                        max(data.par$MeanBeeDiversity),
                        0.1)),
           FlowerDiversity= mean(data.par$FlowerDiversity),
           MeanBeeAbund= mean(data.par$MeanBeeAbund),
           MeanBloomAbund= mean(data.par$MeanBloomAbund),
           Elev= mean(data.par$Elev),
           Acres= mean(data.par$Acres),
           Age_LandTrendr= mean(data.par$Age_LandTrendr)
           ) %>%
    unnest(MeanBeeDiversity) %>%
    mutate(pred = exp(
               ## b_MeanBloomAbund_Intercept+
               ## b_FlowerDiversity_Intercept+
               ## b_MeanBeeAbund_Intercept+
               ## b_MeanBeeDiversity_Intercept+
               b_AnyParasite_Intercept+
               ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBloomAbund_Elev*Elev +
               ## b_MeanBloomAbund_OwnerOwnerB +
               ## b_MeanBloomAbund_OwnerOwnerC +
               ## b_MeanBloomAbund_OwnerODF +
               ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
               ## b_FlowerDiversity_Elev*Elev +
               ## b_FlowerDiversity_OwnerOwnerB +
               ## b_FlowerDiversity_OwnerOwnerC +
               ## b_FlowerDiversity_OwnerODF +
               ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
               ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeAbund_Acres*Acres +
               ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
               ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeDiversity_Acres*Acres +
               ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
               b_AnyParasite_FlowerDiversity*FlowerDiversity +
               b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
               b_AnyParasite_MeanBloomAbund*MeanBloomAbund
           )/
               (1+exp(
                      ## b_MeanBloomAbund_Intercept+
                      ## b_FlowerDiversity_Intercept+
                      ## b_MeanBeeAbund_Intercept+
                      ## b_MeanBeeDiversity_Intercept+
                      b_AnyParasite_Intercept+
                      ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBloomAbund_Elev*Elev +
                      ## b_MeanBloomAbund_OwnerOwnerB +
                      ## b_MeanBloomAbund_OwnerOwnerC +
                      ## b_MeanBloomAbund_OwnerODF +
                      ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
                      ## b_FlowerDiversity_Elev*Elev +
                      ## b_FlowerDiversity_OwnerOwnerB +
                      ## b_FlowerDiversity_OwnerOwnerC +
                      ## b_FlowerDiversity_OwnerODF +
                      ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
                      ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeAbund_Acres*Acres +
                      ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
                      ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeDiversity_Acres*Acres +
                      ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
                      b_AnyParasite_FlowerDiversity*FlowerDiversity +
                      b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
                      b_AnyParasite_MeanBloomAbund*MeanBloomAbund
                  ))) %>%
    group_by(MeanBeeDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBeeDiversity, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Parasite Prevalence") +
    xlab("Bee community diversity") +
    scale_x_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #    #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=ParasitismRate, x=MeanBeeDiversity, colour="white"))




p2.parasite  <- fit %>%
    spread_draws(
        b_MeanBloomAbund_Intercept,
        b_FlowerDiversity_Intercept,
        b_MeanBeeAbund_Intercept,
        b_MeanBeeDiversity_Intercept,
        b_AnyParasite_Intercept,
        b_MeanBloomAbund_Age_LandTrendr,
        b_MeanBloomAbund_Elev,
        b_MeanBloomAbund_OwnerOwnerB,
        b_MeanBloomAbund_OwnerOwnerC,
        b_MeanBloomAbund_OwnerODF,
        b_FlowerDiversity_Age_LandTrendr,
        b_FlowerDiversity_IAge_LandTrendrE2,
        b_FlowerDiversity_Elev,
        b_FlowerDiversity_OwnerOwnerB,
        b_FlowerDiversity_OwnerOwnerC,
        b_FlowerDiversity_OwnerODF,
        b_MeanBeeAbund_MeanBloomAbund,
        b_MeanBeeAbund_Age_LandTrendr,
        b_MeanBeeAbund_Acres,
        b_MeanBeeDiversity_FlowerDiversity,
        b_MeanBeeDiversity_Age_LandTrendr,
        b_MeanBeeDiversity_Acres,
        b_AnyParasite_MeanBeeAbund,
        b_AnyParasite_FlowerDiversity,
        b_AnyParasite_MeanBeeDiversity,
        b_AnyParasite_MeanBloomAbund,
        `b_MeanBeeAbund_Age_LandTrendr:Acres`,
        `b_MeanBeeDiversity_Age_LandTrendr:Acres`
    ) %>%
    mutate(MeanBeeAbund =
               list(seq(min(data.par$MeanBeeAbund),
                        max(data.par$MeanBeeAbund),
                        0.1)),
           FlowerDiversity= mean(data.par$FlowerDiversity),
           MeanBeeDiversity= mean(data.par$MeanBeeDiversity),
           MeanBloomAbund= mean(data.par$MeanBloomAbund),
           Elev= mean(data.par$Elev),
           Acres= mean(data.par$Acres),
           Age_LandTrendr= mean(data.par$Age_LandTrendr)
           ) %>%
    unnest(MeanBeeAbund) %>%
    mutate(pred = exp(
               ## b_MeanBloomAbund_Intercept+
               ## b_FlowerDiversity_Intercept+
               ## b_MeanBeeAbund_Intercept+
               ## b_MeanBeeDiversity_Intercept+
               b_AnyParasite_Intercept+
               ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBloomAbund_Elev*Elev +
               ## b_MeanBloomAbund_OwnerOwnerB +
               ## b_MeanBloomAbund_OwnerOwnerC +
               ## b_MeanBloomAbund_OwnerODF +
               ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
               ## b_FlowerDiversity_Elev*Elev +
               ## b_FlowerDiversity_OwnerOwnerB +
               ## b_FlowerDiversity_OwnerOwnerC +
               ## b_FlowerDiversity_OwnerODF +
               ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
               ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeAbund_Acres*Acres +
               ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
               ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeDiversity_Acres*Acres +
               ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
               b_AnyParasite_FlowerDiversity*FlowerDiversity +
               b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
               b_AnyParasite_MeanBloomAbund*MeanBloomAbund
           )/
               (1+exp(
                      ## b_MeanBloomAbund_Intercept+
                      ## b_FlowerDiversity_Intercept+
                      ## b_MeanBeeAbund_Intercept+
                      ## b_MeanBeeDiversity_Intercept+
                      b_AnyParasite_Intercept+
                      ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBloomAbund_Elev*Elev +
                      ## b_MeanBloomAbund_OwnerOwnerB +
                      ## b_MeanBloomAbund_OwnerOwnerC +
                      ## b_MeanBloomAbund_OwnerODF +
                      ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
                      ## b_FlowerDiversity_Elev*Elev +
                      ## b_FlowerDiversity_OwnerOwnerB +
                      ## b_FlowerDiversity_OwnerOwnerC +
                      ## b_FlowerDiversity_OwnerODF +
                      ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
                      ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeAbund_Acres*Acres +
                      ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
                      ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeDiversity_Acres*Acres +
                      ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
                      b_AnyParasite_FlowerDiversity*FlowerDiversity +
                      b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
                      b_AnyParasite_MeanBloomAbund*MeanBloomAbund
                  ))) %>%
    group_by(MeanBeeAbund) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBeeAbund, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Parasite Prevalence") +
    xlab("Bee abundance") +
    scale_x_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #    #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=ParasitismRate, x=MeanBeeAbund, colour="white"))


ggsave(p1.parasite, file="figures/parasite_beeDiv_dark.pdf",
       height=4, width=5)

ggsave(p2.parasite, file="figures/parasite_beeAbund_dark.pdf",
       height=4, width=5)

parasite.all <- grid.arrange(p1.parasite, p2.parasite, ncol=2)

ggsave(parasite.all, file="figures/all_parasite_dark.pdf",
       height=4, width=10)


## ***********************************************************************
## bee offspring
## ***********************************************************************

offspring.quantiles <- quantile(repro.block$SumOffspring, c(0.025, 0.975))

## flower diversity
p1.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept,
                 b_SumOffspring_FlowerDiversity) %>%
    mutate(FlowerDiversity =
               list(seq(min(repro.block$FlowerDiversity),
                        max(repro.block$FlowerDiversity),
                        0.1))) %>%
    unnest(FlowerDiversity) %>%
    mutate(pred = exp(b_SumOffspring_Intercept +
                      b_SumOffspring_FlowerDiversity*FlowerDiversity)) %>%
    group_by(FlowerDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = FlowerDiversity, y = pred_m)) +
    geom_line(color="white") +
    scale_x_continuous(
        breaks = axis.flower.div.repro,
        labels =  labs.flower.div.repro) +
    ## scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    ##               labels = trans_format("log10", math_format(10^.x)),
    ##               limits=c(1,10^4))+
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Osmia offspring") +
    xlab("Floral community diversity") +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=repro.block[repro.block$SumOffspring >=
                                offspring.quantiles[[1]] &
                                repro.block$SumOffspring <= offspring.quantiles[[2]],],
               aes(y=SumOffspring, x=FlowerDiversity, color="white"))

ggsave(p1.offspring, file="figures/offspring_floraldiv.pdf",
       height=4, width=5)

## bee abund
p2.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept,
                 b_SumOffspring_MeanBeeAbund) %>%
    mutate(MeanBeeAbund =
               list(seq(round(range(repro.block$MeanBeeAbund)[1],1),
                        round(range(repro.block$MeanBeeAbund)[2], 1),
                        1))) %>%
    unnest(MeanBeeAbund) %>%
    mutate(pred = exp(b_SumOffspring_Intercept +
                      b_SumOffspring_MeanBeeAbund*MeanBeeAbund)) %>%
    group_by(MeanBeeAbund) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBeeAbund, y = pred_m)) +
    geom_line(color="white") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(1,10^4)) +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("") +
    xlab("Bee abundance") +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=repro.block,
               aes(y=SumOffspring, x=MeanBeeAbund, color="white"))

ggsave(p2.offspring, file="figures/offspring_beeAbund.pdf",
       height=4, width=5)

offspring <- grid.arrange(p1.offspring, p2.offspring, ncol=2)

ggsave(offspring, file="figures/offspring.pdf",
       height=4, width=10)


## ***********************************************************************
## bee community- floral div, acres
## ***********************************************************************

## acres
p1.bee  <- fit %>%
    spread_draws(
        b_MeanBloomAbund_Intercept,
        b_FlowerDiversity_Intercept,
        b_MeanBeeAbund_Intercept,
        b_MeanBeeDiversity_Intercept,
        b_AnyParasite_Intercept,
        b_MeanBloomAbund_Age_LandTrendr,
        b_MeanBloomAbund_Elev,
        b_MeanBloomAbund_OwnerOwnerB,
        b_MeanBloomAbund_OwnerOwnerC,
        b_MeanBloomAbund_OwnerODF,
        b_FlowerDiversity_Age_LandTrendr,
        b_FlowerDiversity_IAge_LandTrendrE2,
        b_FlowerDiversity_Elev,
        b_FlowerDiversity_OwnerOwnerB,
        b_FlowerDiversity_OwnerOwnerC,
        b_FlowerDiversity_OwnerODF,
        b_MeanBeeAbund_MeanBloomAbund,
        b_MeanBeeAbund_Age_LandTrendr,
        b_MeanBeeAbund_Acres,
        b_MeanBeeDiversity_FlowerDiversity,
        b_MeanBeeDiversity_Age_LandTrendr,
        b_MeanBeeDiversity_Acres,
        b_AnyParasite_MeanBeeAbund,
        b_AnyParasite_FlowerDiversity,
        b_AnyParasite_MeanBeeDiversity,
        b_AnyParasite_MeanBloomAbund,
        `b_MeanBeeAbund_Age_LandTrendr:Acres`,
        `b_MeanBeeDiversity_Age_LandTrendr:Acres`
    ) %>%
    mutate(Acres =
               list(seq(min(indiv.data$Acres),
                        max(indiv.data$Acres),
                        0.1)),
           FlowerDiversity= mean(indiv.data$FlowerDiversity),
           MeanBeeAbund= mean(indiv.data$MeanBeeAbund),
           MeanBloomAbund= mean(indiv.data$MeanBloomAbund),
           Elev= mean(indiv.data$Elev),
           MeanBeeDiversity= mean(indiv.data$MeanBeeDiversity),
           Age_LandTrendr= mean(indiv.data$Age_LandTrendr)
           ) %>%
    unnest(Acres) %>%
    mutate(pred = exp(
               ## b_MeanBloomAbund_Intercept+
               ## b_FlowerDiversity_Intercept+
               b_MeanBeeAbund_Intercept+
               ## b_MeanBeeDiversity_Intercept+
               ## b_AnyParasite_Intercept+
               ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBloomAbund_Elev*Elev +
               ## b_MeanBloomAbund_OwnerOwnerB +
               ## b_MeanBloomAbund_OwnerOwnerC +
               ## b_MeanBloomAbund_OwnerODF +
               ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
               ## b_FlowerDiversity_Elev*Elev +
               ## b_FlowerDiversity_OwnerOwnerB +
               ## b_FlowerDiversity_OwnerOwnerC +
               ## b_FlowerDiversity_OwnerODF +
               b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
               b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
               b_MeanBeeAbund_Acres*Acres +
               ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
               ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeDiversity_Acres*Acres +
               `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres
               ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               ## b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
               ## b_AnyParasite_FlowerDiversity*FlowerDiversity +
               ## b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
               ## b_AnyParasite_MeanBloomAbund*MeanBloomAbund
           )/
               (1+exp(
                      ## b_MeanBloomAbund_Intercept+
                      ## b_FlowerDiversity_Intercept+
                      b_MeanBeeAbund_Intercept+
                      ## b_MeanBeeDiversity_Intercept+
                      ## b_AnyParasite_Intercept+
                      ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBloomAbund_Elev*Elev +
                      ## b_MeanBloomAbund_OwnerOwnerB +
                      ## b_MeanBloomAbund_OwnerOwnerC +
                      ## b_MeanBloomAbund_OwnerODF +
                      ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
                      ## b_FlowerDiversity_Elev*Elev +
                      ## b_FlowerDiversity_OwnerOwnerB +
                      ## b_FlowerDiversity_OwnerOwnerC +
                      ## b_FlowerDiversity_OwnerODF +
                      b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
                      b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
                      b_MeanBeeAbund_Acres*Acres +
                      ## b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
                      ## b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeDiversity_Acres*Acres +
                      `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres
                      ## `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      ## b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
                      ## b_AnyParasite_FlowerDiversity*FlowerDiversity +
                      ## b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
                      ## b_AnyParasite_MeanBloomAbund*MeanBloomAbund
                  ))) %>%
    group_by(Acres) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Acres, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee abundance") +
    xlab("Stand area (log acres)") +
    scale_x_continuous(
        breaks = axis.acres,
        labels =  labs.acres) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite) &
                               indiv.data$MeanBeeAbund !=
                               max(indiv.data$MeanBeeAbund),],
               aes(y=MeanBeeAbund, x=Acres, color="white"))

## FloralDiversity
p2.bee  <- fit %>%
    spread_draws(
        b_MeanBloomAbund_Intercept,
        b_FlowerDiversity_Intercept,
        b_MeanBeeAbund_Intercept,
        b_MeanBeeDiversity_Intercept,
        b_AnyParasite_Intercept,
        b_MeanBloomAbund_Age_LandTrendr,
        b_MeanBloomAbund_Elev,
        b_MeanBloomAbund_OwnerOwnerB,
        b_MeanBloomAbund_OwnerOwnerC,
        b_MeanBloomAbund_OwnerODF,
        b_FlowerDiversity_Age_LandTrendr,
        b_FlowerDiversity_IAge_LandTrendrE2,
        b_FlowerDiversity_Elev,
        b_FlowerDiversity_OwnerOwnerB,
        b_FlowerDiversity_OwnerOwnerC,
        b_FlowerDiversity_OwnerODF,
        b_MeanBeeAbund_MeanBloomAbund,
        b_MeanBeeAbund_Age_LandTrendr,
        b_MeanBeeAbund_Acres,
        b_MeanBeeDiversity_FlowerDiversity,
        b_MeanBeeDiversity_Age_LandTrendr,
        b_MeanBeeDiversity_Acres,
        b_AnyParasite_MeanBeeAbund,
        b_AnyParasite_FlowerDiversity,
        b_AnyParasite_MeanBeeDiversity,
        b_AnyParasite_MeanBloomAbund,
        `b_MeanBeeAbund_Age_LandTrendr:Acres`,
        `b_MeanBeeDiversity_Age_LandTrendr:Acres`
    ) %>%
    mutate(FlowerDiversity =
               list(seq(min(indiv.data$FlowerDiversity),
                        max(indiv.data$FlowerDiversity),
                        0.1)),
           Acres= mean(indiv.data$Acres),
           MeanBeeAbund= mean(indiv.data$MeanBeeAbund),
           MeanBloomAbund= mean(indiv.data$MeanBloomAbund),
           Elev= mean(indiv.data$Elev),
           MeanBeeDiversity= mean(indiv.data$MeanBeeDiversity),
           Age_LandTrendr= mean(indiv.data$Age_LandTrendr)
           ) %>%
    unnest(FlowerDiversity) %>%
    mutate(pred = exp(
               ## b_MeanBloomAbund_Intercept+
               ## b_FlowerDiversity_Intercept+
               ## b_MeanBeeAbund_Intercept+
               b_MeanBeeDiversity_Intercept+
               ## b_AnyParasite_Intercept+
               ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBloomAbund_Elev*Elev +
               ## b_MeanBloomAbund_OwnerOwnerB +
               ## b_MeanBloomAbund_OwnerOwnerC +
               ## b_MeanBloomAbund_OwnerODF +
               ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
               ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
               ## b_FlowerDiversity_Elev*Elev +
               ## b_FlowerDiversity_OwnerOwnerB +
               ## b_FlowerDiversity_OwnerOwnerC +
               ## b_FlowerDiversity_OwnerODF +
               ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
               ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
               ## b_MeanBeeAbund_Acres*Acres +
               b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
               b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
               b_MeanBeeDiversity_Acres*Acres +
               ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
               `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres
               ## b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
               ## b_AnyParasite_FlowerDiversity*FlowerDiversity +
               ## b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
               ## b_AnyParasite_MeanBloomAbund*MeanBloomAbund
           )/
               (1+exp(
                      ## b_MeanBloomAbund_Intercept+
                      ## b_FlowerDiversity_Intercept+
                      ## b_MeanBeeAbund_Intercept+
                      b_MeanBeeDiversity_Intercept+
                      ## b_AnyParasite_Intercept+
                      ## b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBloomAbund_Elev*Elev +
                      ## b_MeanBloomAbund_OwnerOwnerB +
                      ## b_MeanBloomAbund_OwnerOwnerC +
                      ## b_MeanBloomAbund_OwnerODF +
                      ## b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
                      ## b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr +
                      ## b_FlowerDiversity_Elev*Elev +
                      ## b_FlowerDiversity_OwnerOwnerB +
                      ## b_FlowerDiversity_OwnerOwnerC +
                      ## b_FlowerDiversity_OwnerODF +
                      ## b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund +
                      ## b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr +
                      ## b_MeanBeeAbund_Acres*Acres +
                      b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity +
                      b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr +
                      b_MeanBeeDiversity_Acres*Acres +
                      ## `b_MeanBeeAbund_Age_LandTrendr:Acres`*Age_LandTrendr*Acres +
                      `b_MeanBeeDiversity_Age_LandTrendr:Acres`*Age_LandTrendr*Acres
                      ## b_AnyParasite_MeanBeeAbund*MeanBeeAbund +
                      ## b_AnyParasite_FlowerDiversity*FlowerDiversity +
                      ## b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity +
                      ## b_AnyParasite_MeanBloomAbund*MeanBloomAbund
                  ))) %>%
    group_by(FlowerDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = FlowerDiversity, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee community diversity") +
    xlab("Floral diversity") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=FlowerDiversity, color="white"))


bee.plots <- grid.arrange(p1.bee, p2.bee, ncol=1)

ggsave(p1.bee, file="figures/bee_acres.pdf",
       height=4, width=5)

ggsave(p2.bee, file="figures/bee_floralDiv.pdf",
       height=4, width=5)

ggsave(bee.plots, file="figures/beeComm.pdf",
       height=8, width=4)


## ***********************************************************************
## bee community- acres-age interaction
## ***********************************************************************
## based on "Statistical Rethinking example"

ages  <- seq(min(indiv.data.orig1$Age_LandTrendr),
                                   max(indiv.data.orig1$Age_LandTrendr),
                                   length.out=10)
ic <-
    list(Age_LandTrendr = standardize(log(ages)))


cond.effects <- conditional_effects(fit,
                                  effects = "Acres:Age_LandTrendr",
                                  int_conditions = ic)

p.interact  <- plot(cond.effects, plot=FALSE, points=TRUE)

cols <- viridis(length(ic[[1]]), end=0.9)

p.interact[[3]] +   ylab("Bee abundance") +
    xlab("Stand acres (log)") +
    scale_x_continuous(
        breaks = axis.acres,
        labels =  labs.acres) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_classic() +
    ## theme_dark_black()+
    scale_colour_discrete(name="Yrs post harvest",
                          labels=ages,
                          type=cols) +
    ## scale_fill_manual(values = alpha(rep("white", length(ic[[1]])),
    ##                                  1)) +
    scale_fill_viridis_d(alpha=0.9, option="viridis", end=0.9) +
    guides(fill = "none")



ic <-
  list(Acres =  axis.acres)

cond.effects <- conditional_effects(fit,
                                  effects = "Age_LandTrendr:Acres",
                                  int_conditions = ic)

p.interact  <- plot(points=T, cond.effects, plot=FALSE)

cols <- viridis(length(axis.acres), end=0.9)

p.interact[[3]] +   ylab("Bee abundance") +
    xlab("Years post harvest (log)") +
    scale_x_continuous(
        breaks = axis.age.x,
        labels =  labs.age.x) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        theme_classic() +
## theme_dark_black()+
    scale_colour_discrete(name="Acres",
                          labels=round(exp(labs.acres), 0),
                          type=cols) +
    scale_fill_viridis_d(alpha=0.01, option="viridis", end=0.9) +
    guides(fill = "none")


## ***********************************************************************
## bee community- age
## ***********************************************************************


p1.flower.age <- fit %>%
    spread_draws(b_MeanBloomAbund_Intercept,
                 b_MeanBloomAbund_Age_LandTrendr) %>%
    mutate(Age_LandTrendr =
               list(seq(min(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        0.1))) %>%
    unnest(Age_LandTrendr) %>%
    mutate(pred = b_MeanBloomAbund_Intercept +
               b_MeanBloomAbund_Age_LandTrendr*Age_LandTrendr) %>%
    group_by(Age_LandTrendr) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age_LandTrendr, y = pred_m)) +
    geom_line(color="black") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Floral abundance") +
    xlab("") +
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
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBloomAbund, x=Age_LandTrendr, color="black"))

p2.flower.age <- fit %>%
    spread_draws(b_FlowerDiversity_Intercept,
                 b_FlowerDiversity_Age_LandTrendr,
                 b_FlowerDiversity_IAge_LandTrendrE2
                 ) %>%
    mutate(Age_LandTrendr =
               list(seq(min(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        0.1))) %>%
    unnest(Age_LandTrendr) %>%
    mutate(pred = b_FlowerDiversity_Intercept +
               b_FlowerDiversity_Age_LandTrendr*Age_LandTrendr +
               b_FlowerDiversity_IAge_LandTrendrE2*Age_LandTrendr) %>%
    group_by(Age_LandTrendr) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age_LandTrendr, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Floral Diversity") +
    xlab("") +
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
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=FlowerDiversity, x=Age_LandTrendr, color="white"))


p3.bee.age <- fit %>%
    spread_draws(b_MeanBeeAbund_Intercept,
                 b_MeanBeeAbund_Age_LandTrendr) %>%
    mutate(Age_LandTrendr =
               list(seq(min(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        0.1))) %>%
    unnest(Age_LandTrendr) %>%
    mutate(pred = b_MeanBeeAbund_Intercept +
               b_MeanBeeAbund_Age_LandTrendr*Age_LandTrendr) %>%
    group_by(Age_LandTrendr) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age_LandTrendr, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee abundance") +
    xlab("Years post harvest (log)") +
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
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeAbund, x=Age_LandTrendr, color="white"))

p4.bee.age <- fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_MeanBeeDiversity_Age_LandTrendr) %>%
    mutate(Age_LandTrendr =
               list(seq(min(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        0.1))) %>%
    unnest(Age_LandTrendr) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept +
               b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr) %>%
    group_by(Age_LandTrendr) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age_LandTrendr, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee Diversity") +
    xlab("Years post harvest (log)") +
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
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=Age_LandTrendr, color="white"))

age.plots <- grid.arrange(p1.flower.age, p2.flower.age, p3.bee.age,
                          p4.bee.age,
                          ncol=2)

ggsave(age.plots, file="figures/standage.pdf",
       height=7, width=8)




p5.parasite.age <- ggplot(data= indiv.data[indiv.data$Weights == 1,],
                          aes(x=Age_LandTrendr, y=ParasitismRate)) +
    geom_point() +
    geom_smooth(method=lm)


fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_MeanBeeDiversity_Age_LandTrendr) %>%
    mutate(Age_LandTrendr =
               list(seq(min(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age_LandTrendr[!is.na(indiv.data$AnyParasite)]),
                        0.1))) %>%
    unnest(Age_LandTrendr) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept +
               b_MeanBeeDiversity_Age_LandTrendr*Age_LandTrendr) %>%
    group_by(Age_LandTrendr) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age_LandTrendr, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee Diversity") +
    xlab("Years post harvest") +
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
                                        #theme_classic() +
    theme_dark_black()+
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=Age_LandTrendr, color="white"))


## ***********************************************************************
## Landowner
## ***********************************************************************
## not that useful of a figure because there "age" is not accounted
## for, but a little useful

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
    ## geom_boxplot(width=0.1) +
    scale_fill_brewer(palette="Blues") +     #theme_classic() +
    theme_dark_black()+
    theme(legend.position="none") +
    ylab("Floral community diversity") +
    xlab("") + theme_dark_black()


ggsave(p1.owner, file="figures/flowerDiv_owner.pdf",
       height=4, width=5)


p2.owner <- ggplot(site.data , aes(x=AgeBin, y=MeanBloomAbund)) +
    geom_bar(aes(fill = Owner),
             stat = "identity", position = position_dodge(0.8),
             width = 0.7) +
    coord_flip() +
    ## geom_boxplot(width=0.1) +
    scale_fill_brewer(palette="Blues") +     #theme_classic() +
    theme_dark_black()+
    ##theme(legend.position="none") +
    ylab("Floral abundance") +
    xlab("") + theme_dark_black()


owner.all <- grid.arrange(p1.owner, p2.owner, ncol=2)

ggsave(owner.all, file="figures/all_owner.pdf",
       height=4, width=10)


library(sp)
library(rgdal)
spatial <- site.data

coordinates(spatial) <- cbind(spatial$CenterLon, spatial$CenterLat)
proj4string(spatial) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

writeOGR(obj = spatial,
         dsn="../../forestosmia_saved/spatial/stands.shp", layer="spatial", driver="ESRI Shapefile", overwrite=T)
