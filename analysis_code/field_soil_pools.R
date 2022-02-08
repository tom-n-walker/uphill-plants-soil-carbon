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
library(nlme)
library(emmeans)
library(multcomp)
library(tidyverse)
source("./r_code/mm2in.R")


#### DATA ----------------------------------------------------------------------

## Load from Drake plan ----
soil <- drake::readd(field_data)
fluxes <- drake::readd(flux_data)

## Basic formatting ----
# soil data both sites
soil <- soil %>%
  dplyr::select(site, treatments, soil_pools) %>% 
  unnest(cols = c(treatments, soil_pools)) %>%
  as.data.frame %>%
  # change treatment labels
  mutate(treatment = ifelse(treatment == "I", "WL", treatment)) %>%
  # make site-level blocking factor
  mutate(site_block = tolower(paste0(substr(site, 1, 1), block)))
# subset soil data for detailed microbial measures
soilDeep <- filter(soil, site == "lavey") %>%
  mutate(treatment = factor(treatment))
# add blocking factor fluxes
fluxes <- fluxes %>%
  mutate(site_block = tolower(paste0(substr(site, 1, 1), block))) %>%
  # change treatment labels
  mutate(treatment = ifelse(treatment == "I", "WL", treatment))


#### FORMAT PLOT DATA ----------------------------------------------------------

## Soil pools ----
sumSoil <- soil %>%
  group_by(treatment) %>%
  summarise(
    DOCmean = mean(DOC, na.rm = T),
    DOCse = sd(DOC, na.rm = T)/sqrt(n()),
    CMICmean = mean(Cmic, na.rm = T),
    CMICse = sd(Cmic, na.rm = T)/sqrt(n()),
    RMmean = mean(Rm, na.rm = T),
    RMse = sd(Rm, na.rm = T)/sqrt(n()),
    GMmean = mean(Gm, na.rm = T),
    GMse = sd(Gm, na.rm = T)/sqrt(n()),
    RMMmean = mean(RmM, na.rm = T),
    RMMse = sd(RmM, na.rm = T)/sqrt(n()),
    GMMmean = mean(GmM, na.rm = T),
    GMMse = sd(GmM, na.rm = T)/sqrt(n()),
    CUEmean = mean(CUE, na.rm = T),
    CUEse = sd(CUE, na.rm = T)/sqrt(n())
  ) %>%
  ungroup

## Fluxes ----
sumNEE <- fluxes %>%
  mutate(NEE = NEE/1000) %>%
  group_by(site, treatment) %>%
  summarise(
    NEEmean = mean(NEE, na.rm = T),
    NEEse = sd(NEE, na.rm = T)/sqrt(n()),
  ) %>%
  ungroup %>%
  mutate(site = factor(site, c("lavey", "calanda")))
sumERGPP <- fluxes %>%
  mutate(ER = ER/1000) %>%
  mutate(GPP = GPP/1000) %>%
  group_by(site, treatment) %>%
  summarise(
    ERmean = mean(ER, na.rm = T),
    ERse = sd(ER, na.rm = T)/sqrt(n()),
    GPPmean = mean(GPP, na.rm = T),
    GPPse = sd(GPP, na.rm = T)/sqrt(n()),
  ) %>%
  ungroup %>%
  mutate(site = factor(site, c("lavey", "calanda")))

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

## NEE ----
# build model
m0 <- lme(
  NEE ~ treatment * site + Temp.C + Soil.WVC + PAR, 
  random = ~ 1 | site_block/date, 
  data = fluxes,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r0 <- residuals(m0, type = "normalized")
par(mfrow = c(1, 3))
plot(r0 ~ fitted(m0))
boxplot(r0 ~ fluxes$treatment)
hist(r0)
# test main effects
m0a <- update(m0, ~.- treatment:site)
anova(m0, m0a)
anova(m0, update(m0a, ~.- treatment))
anova(m0, update(m0a, ~.- site))
# post-hoc
m0reml <- update(m0, method = "REML")
m0areml <- update(m0a, method = "REML")
emmeans(m0reml, pairwise ~ treatment | site)
emmeans(m0areml, pairwise ~ treatment | site)

## GPP ----
# build model
m8 <- lme(
  GPP ~ treatment * site + Temp.C + Soil.WVC + PAR, 
  random = ~ 1 | site_block/date, 
  data = fluxes,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r8 <- residuals(m8, type = "normalized")
par(mfrow = c(1, 3))
plot(r8 ~ fitted(m8))
boxplot(r8 ~ fluxes$treatment)
hist(r8)
# test main effects
m8a <- update(m8, ~.- treatment:site)
anova(m8, m8a)
anova(m8, update(m8a, ~.- treatment))
anova(m8, update(m8a, ~.- site))
# post-hoc
m8reml <- update(m8, method = "REML")
m8areml <- update(m8a, method = "REML")
emmeans(m8reml, pairwise ~ treatment | site)
emmeans(m8areml, pairwise ~ treatment | site)


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


#### PLOT ----------------------------------------------------------------------

## Fluxes ----
nee <- ggplot(sumNEE) +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text.x = element_blank()) +
  aes(x = treatment, 
      y = NEEmean, 
      ymin = NEEmean - NEEse, 
      ymax = NEEmean + NEEse) +
  scale_y_continuous(limits = c(NA, 1)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  facet_wrap(~ site) +
  xlab("") +
  ylab(expression(paste("NEE (g C", O[2], " ", m^{-2}, " ", h^{-1}, ")", sep = "")))
gpp <- ggplot(sumERGPP) +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text.x = element_blank()) +
  aes(x = treatment, 
      y = GPPmean, 
      ymin = GPPmean - GPPse, 
      ymax = GPPmean + GPPse) +
  coord_cartesian(ylim = c(NA, 0.5)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  facet_wrap(~ site) +
  xlab("") +
  ylab(expression(paste("GPP (g C", O[2], " ", m^{-2}, " ", h^{-1}, ")", sep = "")))
er <- ggplot(sumERGPP) +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text.x = element_blank()) +
  aes(x = treatment, 
      y = ERmean, 
      ymin = ERmean - ERse, 
      ymax = ERmean + ERse) +
  coord_cartesian(ylim = c(0.5, 2.6)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", fill = "white", col = "black") +
  facet_wrap(~ site) +
  xlab("") +
  ylab(expression(paste("ER (g C", O[2], " ", m^{-2}, " ", h^{-1}, ")", sep = "")))

## Pools ----
cmic <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(3, 5.5)) +
  aes(x = treatment, 
      y = CMICmean, 
      ymax = CMICmean + CMICse, 
      ymin = CMICmean - CMICse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste(C[mic], " (mg C ", g^{-1}, ")")))
rm <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(6, 8)) +
  aes(x = treatment, 
      y = RMmean, 
      ymax = RMmean + RMse, 
      ymin = RMmean - RMse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste("R (µg C ", g^{-1}, " ", h^{-1}, ")", sep = "")))
gm <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.5, 4)) +
  aes(x = treatment, 
      y = GMmean, 
      ymax = GMmean + GMse, 
      ymin = GMmean - GMse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste("G (µg C ", g^{-1}, " ", h^{-1}, ")", sep = "")))
mgm <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0.5, 0.9)) +
  aes(x = treatment, 
      y = GMMmean, 
      ymax = GMMmean + GMMse, 
      ymin = GMMmean - GMMse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste(G[mass], " (ng C ", µg^{-1}, " ", C[mic], " ", h^{-1}, ")", sep = "")))
mrm <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.5, 2.5)) +
  aes(x = treatment, 
      y = RMMmean, 
      ymax = RMMmean + RMMse, 
      ymin = RMMmean - RMMse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab(expression(paste(R[mass], " (ng C ", µg^{-1}, " ", C[mic], " ", h^{-1}, ")", sep = "")))
cue <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(20, 37)) +
  aes(x = treatment, 
      y = CUEmean, 
      ymax = CUEmean + CUEse, 
      ymin = CUEmean - CUEse) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  xlab("") +
  ylab("CUE (%)")


#### BUILD ---------------------------------------------------------------------

postscript(
  file = "./figure_builds/field_fluxes.eps",
  width = mm2in(180),
  height = mm2in(60)
)
cowplot::plot_grid(
  nee, gpp, er,
  nrow = 1,
  align = "hv",
  axis = "tlbr"
)
dev.off()

postscript(
  file = "./figure_builds/field_soil.eps",
  width = mm2in(120),
  height = mm2in(120)
)
cowplot::plot_grid(
  mrm, mgm, cmic, 
  rm, gm, cue,
  nrow = 2,
  align = "hv"
)
dev.off()
