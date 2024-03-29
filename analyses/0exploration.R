## set your working directory to the forestosmia git repo
setwd('forestosmia')

setwd("analyses")
rm(list=ls())

library(ggplot2)
library(lemon)
library(egg)

source("src/misc.R")
source("src/init.R")
source("src/ggplotThemes.R")

load("../data/indivdata.Rdata")
load("../data/reproblock.Rdata")

## parasite totals
par.counts <- table(indiv.data$ParasiteRichness[indiv.data$Apidae == 1])
par.counts <- par.counts/sum(par.counts)
par.counts <- as.data.frame(par.counts)

p <- ggplot(data=par.counts, aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", fill="darkolivegreen") +
    theme_dark_black() + coord_flip() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) +
              xlab("Number of parasites/individual") +
              ylab("Proportion of screened bees")

ggsave(p, file="../../forestosmia_saved/figures/summary_parasites.pdf",
        height=4, width=5)


## Screened individuals per stand
stands <- table(indiv.data$Stand[indiv.data$Apidae == 1])
length(stands)

## number of stands with nest block data
length(unique(repro.block$Stand))

## stands with nest block data but no parasite data
blocks.with.par <- unique(repro.block$Stand[!is.na(repro.block$ParasitismRate)])
names(stands)[!names(stands) %in% blocks.with.par]

