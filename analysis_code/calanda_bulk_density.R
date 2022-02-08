################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Check that bulk density doesn't change
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    25 January 2022
#### ---------------------------------------------------------------------------

#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())
# R session options (no factors, bias against scientific #s)
options(
  stringsAsFactors = F,
  scipen = 6
)

## Libraries ----
# standard library set
library(nlme)
library(emmeans)
library(multcomp)
library(tidyverse)
source("./r_code/mm2in.R")


#### DATA ----------------------------------------------------------------------

## Load ----
bd <- data.table::fread("./data/field_calanda_bd.csv", data.table = F)

## Format ----
bd <- bd %>%
  mutate(treatment = ifelse(treatment == "I", "WL", treatment)) %>%
  mutate(treatment = factor(treatment))


#### ANALYSE -------------------------------------------------------------------

## Model ----
m1 <- lme(
  BD ~ treatment, 
  random = ~ 1 | block,
  data = bd,
  method = "ML"
)

## Validate ----
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ bd$treatment)
hist(r1)

#### Significance ----
anova(m1, update(m1, ~.- treatment))
m1reml <- update(m1, method = "REML")
emmeans(m1reml, pairwise ~ treatment)
summary(glht(m1reml, mcp(treatment = "Tukey")))


#### PLOT ----------------------------------------------------------------------

## Create dataset ----
sumBD <- bd %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(BD, na.rm = T),
    se = sd(BD, na.rm = T)/sqrt(n())
  )

#### Plot ----
bdPlot <- ggplot(sumBD) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0.5, 0.675)) +
  aes(x = treatment, 
      y = mean, 
      ymax = mean + se, 
      ymin = mean - se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste("Bulk density (g ", cm^{-3}, ")", sep = "")))
bdPlot

#### Write file ----

postscript(
  file = "./figure_builds/bulk_density.eps",
  width = mm2in(40),
  height = mm2in(60)
)
bdPlot
dev.off()

