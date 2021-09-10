setwd('/Volumes/bombus/Dropbox (University of Oregon)/forestosmia')
setwd("analyses")
rm(list=ls())
library(ggplot2)
library(viridis)
library(brms)
library(bayesplot)
library(tidybayes)
library(tidyverse)
library(gridExtra)
library(grid)

load(file="saved/offspringFitMod.Rdata")
load(file="saved/parasiteFitMod.Rdata")

## ***********************************************************************
## plotting
## ***********************************************************************

stanplot(fit,
         type = "areas",
         prob = 0.95)


## msms <- marginal_smooths(fit)
## plot(msms)

## ***********************************************************************
## bee community diversity and parasitism
## ***********************************************************************

p1.parasite  <- fit %>%
  spread_draws(b_MeanBeeDiversity_Intercept,
               b_AnyParasite_Intercept,
               b_AnyParasite_MeanBeeDiversity) %>%
  mutate(MeanBeeDiversity =
             list(seq(range(indiv.data$MeanBeeDiversity)[1],
                      range(indiv.data$MeanBeeDiversity)[2] + 0.25,
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
  fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                 fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                 fill="darkolivegreen") +
    ylab("Parasite Prevalence") +
      xlab("Bee community diversity") +
    ylim(0,1)  +
    xlim(range(indiv.data$MeanBeeDiversity))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_point(data=indiv.data[indiv.data$Weights == 1,],
               aes(y=ParasitismRate, x=MeanBeeDiversity))


ggsave("figures/parasite_beeDiversity.pdf",
       height=4, width=5)

## ***********************************************************************
## bee community abundance and parasitism
## ***********************************************************************

p2.parasite  <- fit %>%
  spread_draws(b_MeanBeeAbund_Intercept,
               b_AnyParasite_Intercept,
               b_AnyParasite_MeanBeeAbund) %>%
  mutate(MeanBeeAbund =
             list(seq(range(indiv.data$MeanBeeAbund)[1],
                      range(indiv.data$MeanBeeAbund)[2] + 0.25,
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
  fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                 fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                 fill="darkolivegreen") +
    ylab("Parasite Prevalence") +
      xlab("Bee abundance") +
    ylim(0,1)  +
    xlim(range(indiv.data$MeanBeeAbund))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_point(data=indiv.data[indiv.data$Weights == 1,],
               aes(y=ParasitismRate, x=MeanBeeAbund))


ggsave("figures/parasite_beeAbund.pdf",
       height=4, width=5)


parasite.all <- grid.arrange(p1.parasite, p2.parasite, ncol=2)

ggsave(parasite.all, file="figures/all_parasite.pdf",
       height=4, width=10)


## ***********************************************************************
## bee offspring
## ***********************************************************************

## flower diversity
p1.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept, b_FlowerDiversity_Intercept,
                 b_SumOffspring_FlowerDiversity) %>%
  mutate(FlowerDiversity =
             list(seq(min(repro.block$FlowerDiversity),
                      max(repro.block$FlowerDiversity),
                      0.01))) %>%
  unnest(FlowerDiversity) %>%
    mutate(pred = exp(b_SumOffspring_Intercept + b_FlowerDiversity_Intercept +
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
  geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
  fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                 fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                 fill="darkolivegreen") +
    ylab("Osmia offspring") +
    xlab("Floral community diversity") +
    ylim(range(c(3000, repro.block$SumOffspring)))  +
    xlim(range(repro.block$FlowerDiversity)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
      geom_point(data=repro.block,
                 aes(y=SumOffspring, x=FlowerDiversity))

## bee abund
p2.offspring <- fit2 %>%
    spread_draws(b_SumOffspring_Intercept, b_MeanBeeAbund_Intercept,
                 b_SumOffspring_MeanBeeAbund) %>%
    mutate(MeanBeeAbund =
               list(seq(range(repro.block$MeanBeeAbund)[1],
                        range(repro.block$MeanBeeAbund)[2],
                        0.01))) %>%
    unnest(MeanBeeAbund) %>%
    mutate(pred = exp(b_SumOffspring_Intercept + b_MeanBeeAbund_Intercept +
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
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="darkolivegreen") +
    ylab("") +
    xlab("Bee abundance") +
    ylim(range(c(3000, repro.block$SumOffspring)))  +
    xlim(range(repro.block$MeanBeeAbund)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
      geom_point(data=repro.block,
                 aes(y=SumOffspring, x=MeanBeeAbund))

offspring <- grid.arrange(p1.offspring, p2.offspring, ncol=2)

ggsave(offspring, file="figures/offspring.pdf",
       height=4, width=10)






## ***********************************************************************
## bee community
## ***********************************************************************

## bee abundance and bloom abundance

p0.bee <- fit %>%
    spread_draws(b_MeanBeeAbund_Intercept, b_MeanBloomAbund_Intercept,
                 b_MeanBeeAbund_MeanBloomAbund) %>%
    mutate(MeanBloomAbund =
               list(seq(range(indiv.data$MeanBloomAbund)[1],
                        range(indiv.data$MeanBloomAbund)[2],
                        0.01))) %>%
    unnest(MeanBloomAbund) %>%
    mutate(pred = b_MeanBeeAbund_Intercept + b_MeanBloomAbund_Intercept +
               b_MeanBeeAbund_MeanBloomAbund*MeanBloomAbund) %>%
    group_by(MeanBloomAbund) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanBloomAbund, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="darkolivegreen") +
    ylab("Bee abundance") +
    xlab("Floral abundance") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_point(data=indiv.data[indiv.data$Weights == 1,],
               aes(y=MeanBeeAbund, x=MeanBloomAbund))

ggsave("figures/beeAbund_bloomAbund.pdf",
       height=4, width=5)

## ***********************************************************************
## bee diversity and tree richness

p1.bee <- fit %>%
    spread_draws(b_MeanBeeDiversity_Intercept,
                 b_MeanBeeDiversity_TreeRichness) %>%
    mutate(TreeRichness =
               list(seq(range(indiv.data$TreeRichness)[1],
                        range(indiv.data$TreeRichness)[2],
                        0.01))) %>%
    unnest(TreeRichness) %>%
    mutate(pred = b_MeanBeeDiversity_Intercept +
               b_MeanBeeDiversity_TreeRichness*TreeRichness) %>%
    group_by(TreeRichness) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = TreeRichness, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="darkolivegreen") +
    ylab("Bee community diversity") +
    xlab("") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_point(data=indiv.data[indiv.data$Weights == 1,],
               aes(y=MeanBeeDiversity, x=TreeRichness))

## ***********************************************************************
## flower diversity and tree richness

p1.flower <- fit %>%
    spread_draws(b_FlowerDiversity_Intercept,
                 b_FlowerDiversity_TreeRichness) %>%
    mutate(TreeRichness =
               list(seq(range(indiv.data$TreeRichness)[1],
                        range(indiv.data$TreeRichness)[2],
                        0.01))) %>%
    unnest(TreeRichness) %>%
    mutate(pred = b_FlowerDiversity_Intercept +
               b_FlowerDiversity_TreeRichness*TreeRichness) %>%
    group_by(TreeRichness) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = TreeRichness, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="darkolivegreen") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="darkolivegreen") +
    ylab("Floral diversity") +
    xlab("Tree richness") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_point(data=indiv.data[indiv.data$Weights == 1,],
               aes(y=FlowerDiversity, x=TreeRichness))

tree.rich <- grid.arrange(p1.bee, p1.flower, ncol=1)

ggsave(tree.rich, file="figures/bee_floral_treeRichness.pdf",
       height=8, width=5)
