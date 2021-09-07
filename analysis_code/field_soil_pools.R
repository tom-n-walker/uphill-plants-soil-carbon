################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Field soil pools analysis
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
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
library(tidyverse)
library(nlme)
library(emmeans)
library(multcomp)


#### DATA ----------------------------------------------------------------------

## Load from Drake plan ----
soil <- drake::readd(field_data)
fluxes <- drake::readd(flux_data)

## Basic formatting ----
# soil data both sites
soil <- soil %>%
  select(site, treatments, soil_pools) %>% 
  unnest(cols = c(treatments, soil_pools)) %>%
  as.data.frame %>%
  # make site-level blocking factor
  mutate(site_block = tolower(paste0(substr(site, 1, 1), block)))
# subset soil data for detailed microbial measures
soilDeep <- filter(soil, site == "lavey")
# add blocking factor fluxes
fluxes <- fluxes %>%
  mutate(site_block = tolower(paste0(substr(site, 1, 1), block)))


#### ANALYSE -------------------------------------------------------------------

## Ecosystem respiration ----
# build model
m1 <- lme(
  ER ~ treatment * site + Temp.C + Soil.WVC, 
  random = ~ 1 | site_block/date, 
  data = fluxes,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ fluxes$treatment)
hist(r1)
# test main effects
m1a <- update(m1, ~.- treatment:site)
anova(m1, m1a)
anova(m1, update(m1a, ~.- treatment))
anova(m1, update(m1a, ~.- site))
# post-hoc
m1reml <- update(m1, method = "REML")
m1areml <- update(m1a, method = "REML")
emmeans(m1reml, pairwise ~ treatment | site)
emmeans(m1areml, pairwise ~ treatment | site)

## Microbial biomass C ----
# build model
m2 <- lme(
  Cmic ~ treatment * site, 
  random = ~ 1 | site_block, 
  data = soil,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r2 <- residuals(m2, type = "normalized")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
boxplot(r2 ~ soil$treatment)
hist(r2)
# test main effects
m2a <- update(m2, ~.- treatment:site)
anova(m2, m2a)
anova(m2, update(m2a, ~.- treatment))
anova(m2, update(m2a, ~.- site))
# post-hoc
m2areml <- update(m2a, method = "REML")
emmeans(m2areml, pairwise ~ treatment | site)

## Microbial per-gram growth (west Alps only) ----
# build model
m3 <- lme(
  Gm ~ treatment, 
  random = ~ 1 | site_block, 
  data = soilDeep,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r3 <- residuals(m3, type = "normalized")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
boxplot(r3 ~ soilDeep$treatment)
hist(r3)
# test main effects
anova(m3, update(m3, ~.- treatment))
# post-hoc
m3reml <- update(m3, method = "REML")
emmeans(m3reml, pairwise ~ treatment)

## Microbial per-gram respiration (west Alps only) ----
# build model
m4 <- lme(
  Rm ~ treatment, 
  random = ~ 1 | site_block, 
  data = soilDeep,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r4 <- residuals(m4, type = "normalized")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
boxplot(r4 ~ soilDeep$treatment)
hist(r4)
# test main effects
anova(m4, update(m4, ~.- treatment))

## Microbial per-capita growth (west Alps only) ----
# build model
m5 <- lme(
  GmM ~ treatment, 
  random = ~ 1 | site_block, 
  data = soilDeep,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r5 <- residuals(m5, type = "normalized")
par(mfrow = c(1, 3))
plot(r5 ~ fitted(m5))
boxplot(r5 ~ soilDeep$treatment)
hist(r5)
# test main effects
anova(m5, update(m5, ~.- treatment))
# post-hoc
m5reml <- update(m5, method = "REML")
emmeans(m5reml, pairwise ~ treatment)
summary(glht(m5reml, mcp(treatment = "Tukey")))

## Microbial per-capita respiration (west Alps only) ----
# build model
m6 <- lme(
  RmM ~ treatment, 
  random = ~ 1 | site_block, 
  data = soilDeep,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r6 <- residuals(m6, type = "normalized")
par(mfrow = c(1, 3))
plot(r6 ~ fitted(m6))
boxplot(r6 ~ soilDeep$treatment)
hist(r6)
# test main effects
anova(m6, update(m6, ~.- treatment))
# post-hoc
m6reml <- update(m6, method = "REML")
emmeans(m6reml, pairwise ~ treatment)
summary(glht(m6reml, mcp(treatment = "Tukey")))

## Microbial CUE (west Alps only) ----
# build model
m7 <- lme(
  CUE ~ treatment, 
  random = ~ 1 | site_block, 
  data = soilDeep,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r7 <- residuals(m7, type = "normalized")
par(mfrow = c(1, 3))
plot(r7 ~ fitted(m7))
boxplot(r7 ~ soilDeep$treatment)
hist(r7)
# test main effects
anova(m7, update(m7, ~.- treatment))
# post-hoc
m7reml <- update(m7, method = "REML")
emmeans(m7reml, pairwise ~ treatment)

