setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')
setwd("analyses")
rm(list=ls())
library(ggplot2)
library(tidyr)
library(viridis)
library(brms)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(gridExtra)
library(grid)
library(scales)
source("src/makeMultiLevelData.R")
source("src/misc.R")

load("../data/indivdata.Rdata")
indiv.data.orig <- indiv.data

load(file="saved/offspringFitMod.Rdata")
load(file="saved/parasiteFitMod.Rdata")

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************

## indiv.data.orig <- indiv.data.orig[!is.na(indiv.data.orig$AnyParasite),]

indiv.data.orig$Age <- log(indiv.data.orig$Age + 1)

labs.age.x <- (pretty(c(indiv.data.orig$Age),
                      n=10))

## labs.age.x <- (pretty(log(indiv.data.orig$Age +1),
##                       n=6))

axis.age.x <-  standardize.axis(labs.age.x, indiv.data.orig$Age)


labs.bloom.abund <- (pretty(c(0, indiv.data.orig$MeanBloomAbund), n=6))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      indiv.data.orig$MeanBloomAbund)


labs.bee.abund2 <- (pretty(c(0,
                               indiv.data.orig$MeanBeeAbund[
                               !is.na(indiv.data.orig$AnyParasite)]),
                             n=6))
axis.bee.abund2 <-  standardize.axis(labs.bee.abund2,
                                      indiv.data.orig$MeanBeeAbund)


labs.flower.div <- (pretty(c(0, indiv.data.orig$FlowerDiversity), n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     indiv.data.orig$FlowerDiversity)

labs.bee.abund <- (pretty(c(0, indiv.data.orig$MeanBeeAbund), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                    indiv.data.orig$MeanBeeAbund)

labs.bee.div <- (pretty(c(0, indiv.data.orig$MeanBeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  indiv.data.orig$MeanBeeDiversity)

## ***********************************************************************
## company
## ***********************************************************************

## not that useful of a figure because there "age" is not accounted
## for

p1.owner <- ggplot(indiv.data[indiv.data$Weights==1,],
                   aes(x=Owner, y=FlowerDiversity)) +
    geom_violin() + coord_flip() +
    geom_boxplot(width=0.1) +
    scale_fill_brewer(palette="Blues") + theme_classic() +
    theme(legend.position="none") +
    scale_y_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    ylab("Floral community diversity") +
    xlab("")


ggsave(p1.owner, file="figures/flowerDiv_owner.pdf",
       height=4, width=5)

p2.owner <- ggplot(indiv.data[indiv.data$Weights==1 &
                              indiv.data$MeanBloomAbund !=
                              max(indiv.data$MeanBloomAbund, na.rm=TRUE) ,],
                   aes(x=Owner, y=MeanBloomAbund)) +
    geom_violin() + coord_flip() +
    geom_boxplot(width=0.1) +
    scale_fill_brewer(palette="Blues") + theme_classic() +
    theme(legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_y_continuous(
        breaks = axis.bloom.abund,
        labels =  labs.bloom.abund ) +
    ylab("Floral abundance")


owner.all <- grid.arrange(p1.owner, p2.owner, ncol=2)

ggsave(owner.all, file="figures/all_owner.pdf",
       height=4, width=10)


## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************

p1.parasite  <- fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_AnyParasite_Intercept,
                 b_AnyParasite_MeanBeeDiversity) %>%
    mutate(MeanBeeDiversity =
               list(seq(min(indiv.data$MeanBeeDiversity[
                                           indiv.data$WeightsPar == 1]),
                        max(indiv.data$MeanBeeDiversity[
                                           indiv.data$WeightsPar == 1]),
                        0.01))) %>%
    unnest(MeanBeeDiversity) %>%
    mutate(pred = exp(b_MeanBeeDiversity_Intercept +
                      b_AnyParasite_Intercept +
                      b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity)/
               (1+exp(b_MeanBeeDiversity_Intercept +
                      b_AnyParasite_Intercept +
                      b_AnyParasite_MeanBeeDiversity*MeanBeeDiversity))) %>%
    group_by(MeanBeeDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBeeDiversity, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("Parasite Prevalence") +
    xlab("Bee community diversity") +
    scale_x_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                               indiv.data$WeightsPar == 1,],
               aes(y=ParasitismRate, x=MeanBeeDiversity))

p2.parasite  <- fit %>%
    spread_draws(b_MeanBeeAbund_Intercept,
                 b_AnyParasite_Intercept,
                 b_AnyParasite_MeanBeeAbund) %>%
    mutate(MeanBeeAbund =
               list(seq(min(indiv.data$MeanBeeAbund[
                                           indiv.data$WeightsPar == 1]),
                        max(indiv.data$MeanBeeAbund[
                                           indiv.data$WeightsPar == 1]),
                        0.01))) %>%
    unnest(MeanBeeAbund) %>%
    mutate(pred = exp(b_MeanBeeAbund_Intercept +
                      b_AnyParasite_Intercept +
                      b_AnyParasite_MeanBeeAbund*MeanBeeAbund)/
               (1+exp(b_MeanBeeAbund_Intercept +
                      b_AnyParasite_Intercept +
                      b_AnyParasite_MeanBeeAbund*MeanBeeAbund))) %>%
    group_by(MeanBeeAbund) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBeeAbund, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("") +
    xlab("Bee abundance") +
    ylim(0,1)  +
    scale_x_continuous(
        breaks = axis.bee.abund2,
        labels =  labs.bee.abund2) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                              indiv.data$WeightsPar == 1,],
               aes(y=ParasitismRate, x=MeanBeeAbund))

parasite.all <- grid.arrange(p1.parasite, p2.parasite, ncol=2)

ggsave(parasite.all, file="figures/all_parasite.pdf",
       height=4, width=10)


## ***********************************************************************
## bee offspring
## ***********************************************************************

## flower diversity
p1.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept,
                 b_SumOffspring_FlowerDiversity) %>%
    mutate(FlowerDiversity =
               list(seq(min(repro.block$FlowerDiversity),
                        max(repro.block$FlowerDiversity),
                        0.01))) %>%
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
    geom_line() +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1,10^4))+
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("Osmia offspring (log)") +
    xlab("Floral community diversity") +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        theme_classic() +
    geom_point(data=repro.block,
               aes(y=SumOffspring, x=FlowerDiversity))

## bee abund
p2.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept,
                 b_SumOffspring_MeanBeeAbund) %>%
    mutate(MeanBeeAbund =
               list(seq(range(repro.block$MeanBeeAbund)[1],
                        range(repro.block$MeanBeeAbund)[2],
                        0.01))) %>%
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
    geom_line() +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(1,10^4)) +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("") +
    xlab("Bee abundance") +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        theme_classic() +
    geom_point(data=repro.block,
               aes(y=SumOffspring, x=MeanBeeAbund))

offspring <- grid.arrange(p1.offspring, p2.offspring, ncol=2)

ggsave(offspring, file="figures/offspring.pdf",
       height=4, width=10)


## ***********************************************************************
## bee community- floral div
## ***********************************************************************

p1.bee <- fit %>%
    spread_draws(b_MeanBeeAbund_Intercept, b_FlowerDiversity_Intercept,
                 b_MeanBeeAbund_FlowerDiversity) %>%
    mutate(FlowerDiversity =
               list(seq(min(indiv.data$FlowerDiversity),
                        max(indiv.data$FlowerDiversity),
                        0.01))) %>%
    unnest(FlowerDiversity) %>%
    mutate(pred = b_MeanBeeAbund_Intercept + b_FlowerDiversity_Intercept +
               b_MeanBeeAbund_FlowerDiversity*FlowerDiversity) %>%
    group_by(FlowerDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = FlowerDiversity, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("Bee abundance") +
    xlab("") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite) &
                               indiv.data$MeanBeeAbund !=
                               max(indiv.data$MeanBeeAbund),],
               aes(y=MeanBeeAbund, x=FlowerDiversity))

p2.bee <- fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept, b_FlowerDiversity_Intercept,
                 b_MeanBeeDiversity_FlowerDiversity) %>%
    mutate(FlowerDiversity =
               list(seq(min(indiv.data$FlowerDiversity),
                        max(indiv.data$FlowerDiversity),
                        0.01))) %>%
    unnest(FlowerDiversity) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept + b_FlowerDiversity_Intercept +
               b_MeanBeeDiversity_FlowerDiversity*FlowerDiversity) %>%
    group_by(FlowerDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = FlowerDiversity, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
    ylab("Bee diversity") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                               !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=FlowerDiversity))

bee.plots <- grid.arrange(p1.bee, p2.bee, ncol=1)

ggsave(bee.plots, file="figures/beeComm.pdf",
       height=8, width=4)


## ***********************************************************************
## bee community- age
## ***********************************************************************


p1.flower.age <- fit %>%
    spread_draws(b_MeanBloomAbund_Intercept,
                 b_MeanBloomAbund_Age) %>%
    mutate(Age =
               list(seq(min(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        0.01))) %>%
    unnest(Age) %>%
    mutate(pred = b_MeanBloomAbund_Intercept +
               b_MeanBloomAbund_Age*Age) %>%
    group_by(Age) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                              !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBloomAbund, x=Age))

p2.flower.age <- fit %>%
    spread_draws(b_FlowerDiversity_Intercept,
                 b_FlowerDiversity_Age) %>%
    mutate(Age =
               list(seq(min(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        0.01))) %>%
    unnest(Age) %>%
    mutate(pred = b_FlowerDiversity_Intercept +
               b_FlowerDiversity_Age*Age) %>%
    group_by(Age) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                             !is.na(indiv.data$AnyParasite),],
               aes(y=FlowerDiversity, x=Age))


p3.bee.age <- fit %>%
    spread_draws(b_MeanBeeAbund_Intercept,
                 b_MeanBeeAbund_Age) %>%
    mutate(Age =
               list(seq(min(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        0.01))) %>%
    unnest(Age) %>%
    mutate(pred = b_MeanBeeAbund_Intercept +
               b_MeanBeeAbund_Age*Age) %>%
    group_by(Age) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1  &
                              !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeAbund, x=Age))

p4.bee.age <- fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_MeanBeeDiversity_Age) %>%
    mutate(Age =
               list(seq(min(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        0.01))) %>%
    unnest(Age) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept +
               b_MeanBeeDiversity_Age*Age) %>%
    group_by(Age) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                             !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=Age))

age.plots <- grid.arrange(p1.flower.age, p2.flower.age, p3.bee.age,
                          p4.bee.age,
                          ncol=2)

ggsave(age.plots, file="figures/standage.pdf",
       height=7, width=8)




p5.parasite.age <- ggplot(data= indiv.data[indiv.data$Weights == 1,],
                          aes(x=Age, y=ParasitismRate)) +
    geom_point() +
    geom_smooth(method=lm)


    fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_MeanBeeDiversity_Age) %>%
    mutate(Age =
               list(seq(min(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        max(indiv.data$Age[!is.na(indiv.data$AnyParasite)]),
                        0.01))) %>%
    unnest(Age) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept +
               b_MeanBeeDiversity_Age*Age) %>%
    group_by(Age) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Age, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="dodgerblue") +
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
    theme_classic() +
    geom_point(data=indiv.data[indiv.data$Weights == 1 &
                             !is.na(indiv.data$AnyParasite),],
               aes(y=MeanBeeDiversity, x=Age))



## ***********************************************************************
## age splines
## ***********************************************************************

## smooth.fit1 <- conditional_smooths(fit)

## age.p1.flower <- smooth.fit1$"mu_MeanBloomAbund: s(Age)" %>%
##     ggplot(aes(x = Age, y = estimate__)) +
##     geom_line() +
##     geom_ribbon(aes(ymin = pmax(axis.bee.abund[1], lower__),
##                     ymax = upper__), alpha=0.2,
##                 fill="darkolivegreen") +
##     ylab("Bloom abundance") +
##     xlab("Years since clearcut") +
##     scale_x_continuous(
##         breaks = axis.age.x,
##         labels =  labs.age.x ) +
##     scale_y_continuous(
##         breaks = axis.bloom.abund,
##         labels =  labs.bloom.abund ) +
##    theme(axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##           text = element_text(size=16)) +
##     theme_classic() +
##     geom_point(data=indiv.data[indiv.data$Weights == 1 &
##                                indiv.data$MeanBloomAbund !=
##                                max(indiv.data$MeanBloomAbund) ,],
##                aes(y=MeanBloomAbund, x=Age))


## age.p2.flower <- smooth.fit1$"mu_FlowerDiversity: s(Age)" %>%
##     ggplot(aes(x = Age, y = estimate__)) +
##     geom_line() +
##     geom_ribbon(aes(ymin = pmax( axis.flower.div[1], lower__),
##                     ymax = upper__), alpha=0.2,
##                 fill="darkolivegreen") +
##     ylab("Flower diversity") +
##     xlab("Years since clearcut") +
##     scale_x_continuous(
##         breaks = axis.age.x,
##         labels =  labs.age.x ) +
##     scale_y_continuous(
##         breaks = axis.flower.div,
##         labels =  labs.flower.div ) +
##   theme(axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##           text = element_text(size=16)) +
##     theme_classic() +
##     geom_point(data=indiv.data[indiv.data$Weights == 1,],
##                aes(y=FlowerDiversity, x=Age))


## flower.age.plots <- grid.arrange(age.p1.flower, age.p2.flower, ncol=1)

## ggsave(flower.age.plots, file="figures/floralCommAge.pdf",
##        height=8, width=4)


## age.p3.bee <- smooth.fit1$"mu_MeanBeeAbund: s(Age)" %>%
##     ggplot(aes(x = Age, y = estimate__)) +
##     geom_line() +
##     geom_ribbon(aes(ymin = pmax(axis.bee.abund[1], lower__),
##                     ymax = upper__), alpha=0.2,
##                 fill="darkolivegreen") +
##     ylab("Bee abundance") +
##     xlab("Years since clearcut") +
##     scale_x_continuous(
##         breaks = axis.age.x,
##         labels =  labs.age.x ) +
##     scale_y_continuous(
##         breaks = axis.bee.abund,
##         labels =  labs.bee.abund ) +
##     theme(
##           axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##         text = element_text(size=16)) +
##         theme_classic() +
##     geom_point(data=indiv.data[indiv.data$Weights == 1,],
##                aes(y=MeanBeeAbund, x=Age))



## age.p4.bee <- smooth.fit1$"mu_MeanBeeDiversity: s(Age)" %>%
##     ggplot(aes(x = Age, y = estimate__)) +
##     geom_line() +
##     geom_ribbon(aes(ymin = pmax(axis.bee.div[1], lower__),
##                     ymax = upper__), alpha=0.2,
##                 fill="darkolivegreen") +
##     ylab("Bee diversity") +
##     xlab("Years since clearcut") +
##     scale_x_continuous(
##         breaks = axis.age.x,
##         labels =  labs.age.x ) +
##     scale_y_continuous(
##         breaks = axis.bee.div,
##         labels =  labs.bee.div ) +
##     theme(
##           axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##         text = element_text(size=16)) +
##         theme_classic() +
##     geom_point(data=indiv.data[indiv.data$Weights == 1,],
##                aes(y=MeanBeeDiversity, x=Age))


## age.plots <- grid.arrange(age.p1.flower, age.p2.flower,
##                           age.p3.bee, age.p4.bee, ncol=1)

## ggsave(age.plots, file="figures/standage.pdf",
##        height=11, width=4)


