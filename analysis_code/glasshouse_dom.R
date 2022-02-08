################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Glasshouse pot data 
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
soil <- drake::readd(gh_soil)

## Format ----
soilPools <- soil$pot_soil


#### FORMAT FOR PLOTS ----------------------------------------------------------

## Soil pools ----
sumSoil <- soilPools %>% 
  # recalculate Cmic in mg g
  mutate(k = k * 24) %>%
  mutate(p = p/1000) %>%
  # group by treatment, calculate means/se, ungroup
  group_by(treatment) %>%
  summarise(a350_mean = mean(DOM.a350, na.rm = T),
            a350_se = sd(DOM.a350, na.rm = T)/sqrt(n()),
            ftot_mean = mean(DOM.Ftot, na.rm = T),
            ftot_se = sd(DOM.Ftot, na.rm = T)/sqrt(n()),
            fi_mean = mean(DOM.FI, na.rm = T),
            fi_se = sd(DOM.FI, na.rm = T)/sqrt(n()),
            dom1_mean = mean(DOM.C1, na.rm = T),
            dom1_se = sd(DOM.C1, na.rm = T)/sqrt(n()),
            dom2_mean = mean(DOM.C2, na.rm = T),
            dom2_se = sd(DOM.C2, na.rm = T)/sqrt(n()),
            dom3_mean = mean(DOM.C3, na.rm = T),
            dom3_se = sd(DOM.C3, na.rm = T)/sqrt(n()),
            dom4_mean = mean(DOM.C4, na.rm = T),
            dom4_se = sd(DOM.C4, na.rm = T)/sqrt(n()),
            dom5_mean = mean(DOM.C5, na.rm = T),
            dom5_se = sd(DOM.C5, na.rm = T)/sqrt(n()),
            dom6_mean = mean(DOM.C6, na.rm = T),
            dom6_se = sd(DOM.C6, na.rm = T)/sqrt(n())) %>%
  ungroup

## get bare lines ----
allBare <- filter(sumSoil, treatment == "B")
sumSoil <- filter(sumSoil, treatment != "B")


#### SOIL BARE TREATMENT BOUNDS ------------------------------------------------

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


#### ANALYSE -------------------------------------------------------------------

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


#### PLOT ----------------------------------------------------------------------

a350 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(20, 140)) +
  aes(x = treatment, 
      y = a350_mean, 
      ymax = a350_mean + a350_se, 
      ymin = a350_mean - a350_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$a350_mean - allBare$a350_se) +
  geom_hline(yintercept = allBare$a350_mean + allBare$a350_se) +
  xlab("") +
  ylab(expression(paste(a[350], " (", m^{-1}, ")", sep = "")))
ftot <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(5, 32)) +
  aes(x = treatment, 
      y = ftot_mean, 
      ymax = ftot_mean + ftot_se, 
      ymin = ftot_mean - ftot_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$ftot_mean - allBare$ftot_se) +
  geom_hline(yintercept = allBare$ftot_mean + allBare$ftot_se) +
  xlab("") +
  ylab(expression(paste(F[tot], " (R. U.)", sep = "")))
fi <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.55, 1.85)) +
  aes(x = treatment, 
      y = fi_mean, 
      ymax = fi_mean + fi_se, 
      ymin = fi_mean - fi_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$fi_mean - allBare$fi_se) +
  geom_hline(yintercept = allBare$fi_mean + allBare$fi_se) +
  xlab("") +
  ylab("Fluorescence Index")
dom1to3 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom1_mean, 
      ymax = dom1_mean + dom1_se, 
      ymin = dom1_mean - dom1_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom1_mean - allBare$dom1_se) +
  geom_hline(yintercept = allBare$dom1_mean + allBare$dom1_se) +
  xlab("") +
  ylab("C3 (humic-like, %)")
dom2to4 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom2_mean, 
      ymax = dom2_mean + dom2_se, 
      ymin = dom2_mean - dom2_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom2_mean - allBare$dom2_se) +
  geom_hline(yintercept = allBare$dom2_mean + allBare$dom2_se) +
  xlab("") +
  ylab("C4 (humic-like, %)")
dom3to6 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom3_mean, 
      ymax = dom3_mean + dom3_se, 
      ymin = dom3_mean - dom3_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom3_mean - allBare$dom3_se) +
  geom_hline(yintercept = allBare$dom3_mean + allBare$dom3_se) +
  xlab("") +
  ylab("C6 (fumic acid-like, %)")
dom4to1 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom4_mean, 
      ymax = dom4_mean + dom4_se, 
      ymin = dom4_mean - dom4_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom4_mean - allBare$dom4_se) +
  geom_hline(yintercept = allBare$dom4_mean + allBare$dom4_se) +
  xlab("") +
  ylab("C1 (protein-like, %)")
dom5to5 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom5_mean, 
      ymax = dom5_mean + dom5_se, 
      ymin = dom5_mean - dom5_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom5_mean - allBare$dom5_se) +
  geom_hline(yintercept = allBare$dom5_mean + allBare$dom5_se) +
  xlab("") +
  ylab("C5 (humic-like, %)")
dom6to2 <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(2, 45)) +
  aes(x = treatment, 
      y = dom6_mean, 
      ymax = dom6_mean + dom6_se, 
      ymin = dom6_mean - dom6_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$dom6_mean - allBare$dom6_se) +
  geom_hline(yintercept = allBare$dom6_mean + allBare$dom6_se) +
  xlab("") +
  ylab("C2 (protein-like, %)")


#### BUILD PANEL ---------------------------------------------------------------

postscript(
    file = "./figure_builds/glasshouse_dom.eps",
    width = mm2in(150),
    height = mm2in(120)
)
cowplot::plot_grid(
  a350, ftot, fi, 
  dom4to1, dom6to2, dom1to3, 
  dom2to4, dom5to5, dom3to6, dom3to6,
  nrow = 2, 
  align = "hv"
)
dev.off()

