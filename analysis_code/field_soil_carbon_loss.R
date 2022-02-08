################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Main soil carbon loss effect
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
library(data.table)
library(nlme)
library(emmeans)


#### DATA ----------------------------------------------------------------------

## Load from Drake plan ----
allData <- drake::readd(field_data)

## Basic combine and unnest ----
soil <- allData %>%
  select(site, treatments, soil_pools) %>% 
  unnest(cols = c(treatments, soil_pools)) %>%
  as.data.frame %>%
  # make site-level blocking factor
  mutate(site_block = tolower(paste0(substr(site, 1, 1), block)))

## Duplicate high site control for 2nd analysis ----
# build dataset with C and W
soilCW <- soil %>%
  filter(treatment == "C" | treatment == "W") %>%
  # account for duplicated control
  mutate(type = "warmed")
# build dataset with C and I
soilCI <- soil %>%
  filter(treatment == "C" | treatment == "I") %>%
  # account for duplicated control
  mutate(type = "warmed+invaded")
# bind them together (duplicates control)
soilCCWI <- bind_rows(soilCW, soilCI) %>%
  select(site, site_block, type, treatment, Soil.temp, Csoil)


#### ANALYSE -------------------------------------------------------------------

## Main soil carbon effect ----
# build model
m1 <- lme(
  Csoil ~ treatment * site, 
  random = ~ 1 | site_block, 
  data = soil,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ soil$treatment)
hist(r1)
# test main effects
m2 <- update(m1, ~.- treatment:site)
anova(m1, m2)
anova(m2, update(m2, ~.- treatment))
anova(m2, update(m2, ~.- site))
# post-hoc
m1reml <- update(m1, method = "REML")
m2reml <- update(m2, method = "REML")
emmeans(m1reml, pairwise ~ treatment | site)
emmeans(m2reml, pairwise ~ treatment | site)

## Soil carbon on temperature ----
# build model
m4 <- lme(
  Csoil ~ Soil.temp * type,
  random = ~ site | site_block,
  data = soilCCWI,
  method = "ML",
  na.action = "na.exclude"
  )
# diagnose model
r4 <- residuals(m4, type = "pearson")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m4))
boxplot(r4 ~ soilCCWI$treatment)
hist(r4)
# test main effects
m5 <- update(m4, ~.- Soil.temp:type)
anova(m4, m5)

## Calculate magnitude of soil carbon loss ----
# calculating errors following standard practice of proliferating error:
# http://www.met.rdg.ac.uk/~swrhgnrj/combining_errors.pdf
# update best model for REML
m4reml <- update(m4, method = "REML")
# coefficients
m4coefs <- intervals(m4reml, which = "fixed")$fixed
# calculate error 
ciWarm <- (m4coefs[2, 3] - m4coefs[2, 1]) / 2
ciWaLo <- (m4coefs[4, 3] - m4coefs[4, 1]) / 2
# calculate estimates
estWarm <- m4coefs[2, 2]
estWaLo <- estWarm + m4coefs[4, 2]
# relative differences
diffAll <- (estWaLo / estWarm - 1) * 100
ciAll <- sqrt(((ciWarm/estWarm)^2) + ((ciWaLo/estWaLo)^2)) * diffAll
# print excess soil C loss ± 95%
diffAll
ciAll


#### BASIC PLOTS ---------------------------------------------------------------

## Soil C loss on treatment ----
# summarise data
soilPlotData <- soil %>%
  group_by(site, treatment) %>%
  summarise(mean = mean(Csoil, na.rm = T),
            se = sd(Csoil, na.rm = T)/sqrt(n())) %>%
  ungroup %>%
  mutate(treatment = fct_relevel(treatment, "C", "W", "I"))
# plot
mainPlot <- ggplot(soilPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(10, 20)) +
  aes(x = treatment, y = mean, ymax = mean + se, ymin = mean - se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  facet_wrap(~site) +
  labs(y = "Soil C", x = "") 

## Soil C loss per ºC ----
siPlot <- ggplot(soilCCWI) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = Soil.temp, y = Csoil, linetype = type, col = site) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  guides(col = "none") +
  labs(y = "Soil C", x = "Soil ºC") 

## Combine ----
cowplot::plot_grid(mainPlot, siPlot)


## Quantify numbers ----
# get means
meanC <- soilPlotData$mean
seC <- soilPlotData$se
# calanda
meanC[2] - meanC[1]
sqrt(seC[2]^2 + seC[1]^2)
# lavey
meanC[5] - meanC[4]
sqrt(seC[5]^2 + seC[4]^2)
# all




