################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Glasshouse pot-level traits 
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    29 November 2021
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
library(tidyverse)
# custom functions
source("./r_code/mm2in.R")

#### DATA ----------------------------------------------------------------------

## Load ----
ghPots <- drake::readd(gh_soil)$pot_plants

## Format plot summary data ----
sumPots <- ghPots %>%
  # remove bare treatment for plotting (ribbon later)
  filter(treatment != "B") %>%
  # group by treatment to summarise
  group_by(treatment) %>%
  # get mean and SE
  summarise(
    AGB_mean = mean(AGB.g, na.rm = T),
    AGB_se = sd(AGB.g, na.rm = T)/sqrt(n()),
    BGB_mean = mean(BGB.g, na.rm = T),
    BGB_se = sd(BGB.g, na.rm = T)/sqrt(n()),
    R2S_mean = mean(R2S, na.rm = T),
    R2S_se = sd(R2S, na.rm = T)/sqrt(n()),
    SLA_mean = mean(cSLA.cm2.g, na.rm = T),
    SLA_se = sd(cSLA.cm2.g, na.rm = T)/sqrt(n()),
    A_mean = mean(cA.umol.m2.s, na.rm = T),
    A_se = sd(cA.umol.m2.s, na.rm = T)/sqrt(n()),
    G_mean = mean(cGS.umol.m2.s, na.rm = T),
    G_se = sd(cGS.umol.m2.s, na.rm = T)/sqrt(n())
  ) %>%
  ungroup


#### ANALYSE PLANT TRAITS ------------------------------------------------------

## SLA ----
# build model
m1 <- lme(
  cSLA.cm2.g ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ ghPots$treatment)
hist(r1)
# test main effects
anova(m1, update(m1, ~.- treatment))

## A max ----
# build model
m2 <- lme(
  cA.umol.m2.s ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r2 <- residuals(m2, type = "normalized")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
boxplot(r2 ~ ghPots$treatment)
hist(r2)
# test main effects
anova(m2, update(m2, ~.- treatment))

## gs Max ----
# build model
m3 <- lme(
  log10(cGS.umol.m2.s) ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r3 <- residuals(m3, type = "normalized")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
boxplot(r3 ~ ghPots$treatment)
hist(r3)
# test main effects
anova(m3, update(m3, ~.- treatment))

## AGB ----
# build model
m5 <- lme(
  AGB.g ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r5 <- residuals(m5, type = "normalized")
par(mfrow = c(1, 3))
plot(r5 ~ fitted(m5))
boxplot(r5 ~ ghPots$treatment)
hist(r5)
# test main effects
anova(m5, update(m5, ~.- treatment))

## BGB ----
# build model
m6 <- lme(
  BGB.g ~ treatment, 
  random = ~ 1 | block, 
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r6 <- residuals(m6, type = "normalized")
par(mfrow = c(1, 3))
plot(r6 ~ fitted(m6))
boxplot(r6 ~ ghPots$treatment)
hist(r6)
# test main effects
anova(m6, update(m6, ~.- treatment))

## R2S ----
# build model
m7 <- lme(
  R2S ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = ghPots,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r7 <- residuals(m7, type = "normalized")
par(mfrow = c(1, 3))
plot(r7 ~ fitted(m7))
boxplot(r7 ~ ghPots$treatment)
hist(r7)
# test main effects
anova(m7, update(m7, ~.- treatment))


#### PLOT ----------------------------------------------------------------------

sla <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = SLA_mean, 
      ymin = SLA_mean - SLA_se, 
      ymax = SLA_mean + SLA_se) +
  coord_cartesian(ylim = c(150, 175)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab(expression(paste("SLA (c", m^{2}, " ", g^{-1}, ")")))
amax <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = A_mean, 
      ymin = A_mean - A_se, 
      ymax = A_mean + A_se) +
  coord_cartesian(ylim = c(2.5, 8)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab(expression(paste(A[max], " (µmol ", m^{2}, " ", s^{-1}, ")")))
gs <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = G_mean, 
      ymin = G_mean - G_se, 
      ymax = G_mean + G_se) +
  coord_cartesian(ylim = c(20, 100)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab(expression(paste(g[S], " (µmol ", m^{2}, " ", s^{-1}, ")")))
agb <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = AGB_mean, 
      ymin = AGB_mean - AGB_se, 
      ymax = AGB_mean + AGB_se) +
  coord_cartesian(ylim = c(0.5, 2.5)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab("AGB (g)")
bgb <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = BGB_mean, 
      ymin = BGB_mean - BGB_se, 
      ymax = BGB_mean + BGB_se) +
  coord_cartesian(ylim = c(0.5, 2.5)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab("BGB (g)")
r2s <- ggplot(sumPots) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = treatment, 
      y = R2S_mean, 
      ymin = R2S_mean - R2S_se, 
      ymax = R2S_mean + R2S_se) +
  coord_cartesian(ylim = c(0.3, 1)) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  xlab("") +
  ylab("Root:Shoot")

#### BUILD PLOTS ---------------------------------------------------------------

postscript(
  file = "./figure_builds/glasshouse_traits.eps",
  width = mm2in(180),
  height = mm2in(60)
)
cowplot::plot_grid(
  agb, bgb, r2s, sla, amax, gs,
  nrow = 1, 
  align = "hv"
)
dev.off()

