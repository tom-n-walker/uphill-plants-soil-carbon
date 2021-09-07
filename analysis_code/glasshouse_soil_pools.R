################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Glasshouse soil 
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
library(nlme)
library(emmeans)
library(tidyverse)


#### DATA ----------------------------------------------------------------------

## Load from Drake plan ----
soil <- drake::readd(gh_soil)
# select data for analysis
soilPools <- soil$pot_soil
respTime <- soil$mic_resp %>%
  # build categorical hours variable
  mutate(fac_hours = paste0("h", Hours))


#### BARE TREATMENT BOUNDS -----------------------------------------------------

## Calculate ribbon for bare treatment on plots ----
soilPools %>%
  filter(treatment == "B") %>%
  dplyr::select(Cmic.ugC.g:b) %>%
  # calculate upper and lower bounds for all variables
  summarise(
    across(
      everything(), 
      function(x){
        # summary statistics
        mu <- mean(x, na.rm = T)
        n <- sum(!is.na(x))
        sdev <- sd(x, na.rm = T)
        se <- sdev/sqrt(n)
        # calculate upper and lower confidence intervals
        high <- mu + se
        low <- mu - se
        out <- c(high, low)
        return(out)
      }
    )
  )

## Remove bare treatment (only relevant as reference) ----
soilPools <- soilPools %>%
  filter(treatment != "B")


#### ANALYSE RESPIRATION -------------------------------------------------------

# build model
m0 <- lme(
  R.ugC.g.h ~ Treatment * fac_hours, 
  random = ~ 1 | Block,
  data = respTime,
  method = "ML",
  na.action = "na.exclude", 
  weights = NULL
)
# diagnose model
r0 <- residuals(m0, type = "normalized")
par(mfrow = c(1, 4))
plot(r0 ~ fitted(m0))
boxplot(r0 ~ respTime$Treatment)
boxplot(r0 ~ respTime$fac_hours)
hist(r0)
# test main effects
anova(m0, update(m0, ~.- Treatment:fac_hours))
# post-hoc
m0reml <- update(m0, method = "REML")
emmeans(m0, pairwise ~ Treatment | fac_hours)


#### ANALYSE POOLS -------------------------------------------------------------

## DOM a350 ----
# build model
m1 <- lme(
  DOM.a350 ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ soilPools$treatment)
hist(r1)
# test main effects
m1a <- update(m1, ~.- treatment)
anova(m1, m1a)

## DOM total fluorescence ----
# build model
m2 <- lme(
  DOM.Ftot ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r2 <- residuals(m2, type = "normalized")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
boxplot(r2 ~ soilPools$treatment)
hist(r2)
# test main effects
m2a <- update(m2, ~.- treatment)
anova(m2, m2a)

## DOM fluorescence index ----
# build model
m3 <- lme(
  DOM.FI ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r3 <- residuals(m3, type = "normalized")
par(mfrow = c(1, 3))
plot(r3 ~ fitted(m3))
boxplot(r3 ~ soilPools$treatment)
hist(r3)
# test main effects
m3a <- update(m3, ~.- treatment)
anova(m3, m3a)

## Microbial biomass carbon ----
# build model
m4 <- lme(
  log10(Cmic.ugC.g) ~ treatment,  
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r4 <- residuals(m4, type = "normalized")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
boxplot(r4 ~ soilPools$treatment)
hist(r4)
# test main effects
m4a <- update(m4, ~.- treatment)
anova(m4, m4a)

## Fast-decaying pool size ----
# build model
m5 <- lme(
  log10(p) ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r5 <- residuals(m5, type = "normalized")
par(mfrow = c(1, 3))
plot(r5 ~ fitted(m5))
boxplot(r5 ~ soilPools$treatment)
hist(r5)
# test main effects
m5a <- update(m5, ~.- treatment)
anova(m5, m5a)

## Fast-decaying pool rate ----
# build model
m6 <- lme(
  log10(k) ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r6 <- residuals(m6, type = "normalized")
par(mfrow = c(1, 3))
plot(r6 ~ fitted(m6))
boxplot(r6 ~ soilPools$treatment)
hist(r6)
# test main effects
m6a <- update(m6, ~.- treatment)
anova(m6, m6a)

## DOM component 1 (named component #3 in MS) ----
# build model
m7 <- lme(
  DOM.C1 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r7 <- residuals(m7, type = "normalized")
par(mfrow = c(1, 3))
plot(r7 ~ fitted(m7))
boxplot(r7 ~ soilPools$treatment)
hist(r7)
# test main effects
m7a <- update(m7, ~.- treatment)
anova(m7, m7a)

## DOM component 2 (named component #4 in MS) ----
# build model
m8 <- lme(
  DOM.C2 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r8 <- residuals(m8, type = "normalized")
par(mfrow = c(1, 3))
plot(r8 ~ fitted(m8))
boxplot(r8 ~ soilPools$treatment)
hist(r8)
# test main effects
m8a <- update(m8, ~.- treatment)
anova(m8, m8a)

## DOM component 3 (named component #6 in MS) ----
# build model
m9 <- lme(
  DOM.C3 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r9 <- residuals(m9, type = "normalized")
par(mfrow = c(1, 3))
plot(r9 ~ fitted(m9))
boxplot(r9 ~ soilPools$treatment)
hist(r9)
# test main effects
m9a <- update(m9, ~.- treatment)
anova(m9, m9a)

## DOM component 4 (named component #1 in MS) ----
# build model
m10 <- lme(
  DOM.C4 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r10 <- residuals(m10, type = "normalized")
par(mfrow = c(1, 3))
plot(r10 ~ fitted(m10))
boxplot(r10 ~ soilPools$treatment)
hist(r10)
# test main effects
m10a <- update(m10, ~.- treatment)
anova(m10, m10a)

## DOM component 5 (named component #5 in MS) ----
# build model
m11 <- lme(
  DOM.C5 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r11 <- residuals(m11, type = "normalized")
par(mfrow = c(1, 3))
plot(r11 ~ fitted(m11))
boxplot(r11 ~ soilPools$treatment)
hist(r11)
# test main effects
m11a <- update(m11, ~.- treatment)
anova(m11, m11a)

## DOM component 6 (named component #2 in MS) ----
# build model
m12 <- lme(
  DOM.C6 ~ treatment, 
  random = ~ 1 | block, 
  weights = varIdent(form = ~ 1 | treatment),
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r12 <- residuals(m12, type = "normalized")
par(mfrow = c(1, 3))
plot(r12 ~ fitted(m12))
boxplot(r12 ~ soilPools$treatment)
hist(r12)
# test main effects
m12a <- update(m12, ~.- treatment)
anova(m12, m12a)


