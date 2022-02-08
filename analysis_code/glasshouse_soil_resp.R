################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Glasshouse soil data 
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
# select data for analysis
soilPools <- soil$pot_soil
respTime <- soil$mic_resp %>%
  # build categorical hours variable
  mutate(fac_hours = paste0("h", Hours))


#### FORMAT FOR PLOTS ----------------------------------------------------------

## Soil pools ----
sumSoil <- soilPools %>% 
  # recalculate Cmic in mg g
  mutate(Cmic = Cmic.ugC.g/1000) %>%
  mutate(k = k * 24) %>%
  # group by treatment, calculate means/se, ungroup
  group_by(treatment) %>%
  summarise(cmic_mean = mean(Cmic, na.rm = T),
            cmic_se = sd(Cmic, na.rm = T)/sqrt(n()),
            a350_mean = mean(DOM.a350, na.rm = T),
            a350_se = sd(DOM.a350, na.rm = T)/sqrt(n()),
            ftot_mean = mean(DOM.Ftot, na.rm = T),
            ftot_se = sd(DOM.Ftot, na.rm = T)/sqrt(n()),
            fi_mean = mean(DOM.FI, na.rm = T),
            fi_se = sd(DOM.FI, na.rm = T)/sqrt(n()),
            p_mean = mean(p, na.rm = T),
            p_se = sd(p, na.rm = T)/sqrt(n()),
            k_mean = mean(k, na.rm = T),
            k_se = sd(k, na.rm = T)/sqrt(n())) %>%
  ungroup

## Respiration ----
sumCO2 <- respTime %>% 
  # group by treatment and hours, mean/se, ungroup
  group_by(Treatment, Hours) %>%
  # get mean and SE
  summarise(mean = mean(R.ugC.g.h, na.rm = T),
            se = sd(R.ugC.g.h, na.rm = T)/sqrt(n())) %>%
  # generate weeks variable
  mutate(weeks = Hours/24/7) %>%
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


#### ANALYSE SOIL POOLS --------------------------------------------------------

## Microbial biomass carbon ----
# build model
m11 <- lme(
  log10(Cmic.ugC.g) ~ treatment,  
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r4 <- residuals(m11, type = "normalized")
par(mfrow = c(1, 3))
plot(r4 ~ fitted(m11))
boxplot(r4 ~ soilPools$treatment)
hist(r4)
# test main effects
m11a <- update(m11, ~.- treatment)
anova(m11, m11a)

## Fast-decaying pool size ----
# build model
m12 <- lme(
  log10(p) ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r5 <- residuals(m12, type = "normalized")
par(mfrow = c(1, 3))
plot(r5 ~ fitted(m12))
boxplot(r5 ~ soilPools$treatment)
hist(r5)
# test main effects
m12a <- update(m12, ~.- treatment)
anova(m12, m12a)

## Fast-decaying pool rate ----
# build model
m13 <- lme(
  log10(k) ~ treatment, 
  random = ~ 1 | block, 
  data = soilPools,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r6 <- residuals(m13, type = "normalized")
par(mfrow = c(1, 3))
plot(r6 ~ fitted(m13))
boxplot(r6 ~ soilPools$treatment)
hist(r6)
# test main effects
m13a <- update(m13, ~.- treatment)
anova(m13, m13a)


#### PLOT ----------------------------------------------------------------------

## CO2 ----
ts <- ggplot(sumCO2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = weeks, y = mean, ymax = mean + se, ymin = mean - se) +
  guides(linetype = "none", shape = "none") +
  scale_shape_manual(values = c(21, 22, 24)) +
  geom_errorbar(width = 0.1) +
  geom_line(aes(linetype = Treatment)) +
  geom_point(aes(shape = Treatment), fill = "white", size = 1) +
  geom_vline(xintercept = 204/24/7) +
  xlab("Weeks") +
  ylab(expression(paste("R (µg C ", g^{-1}, " ", h^{-1}, ")")))

## Soil pools ----
pp <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(100, 1000)) +
  aes(x = treatment, 
      y = p_mean, 
      ymax = p_mean + p_se, 
      ymin = p_mean - p_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$p_mean - allBare$p_se) +
  geom_hline(yintercept = allBare$p_mean + allBare$p_se) +
  xlab("") +
  ylab(expression(paste("Pool size (µg C ", g^{-1}, ")", sep = "")))
kp <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0.09, 0.38)) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3)) +
  aes(x = treatment, 
      y = k_mean, 
      ymax = k_mean + k_se, 
      ymin = k_mean - k_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$k_mean - allBare$k_se) +
  geom_hline(yintercept = allBare$k_mean + allBare$k_se) +
  xlab("") +
  ylab(expression(paste("Decay (µg C ", g^{-1}, " ", d^{-1}, ")", sep = "")))
cmic <- ggplot(sumSoil) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(ylim = c(1.8, 3.5)) +
  aes(x = treatment, 
      y = cmic_mean, 
      ymax = cmic_mean + cmic_se, 
      ymin = cmic_mean - cmic_se) +
  geom_errorbar(width = 0.2) +
  geom_bar(stat = "identity", col = "black", fill = "white") +
  geom_hline(yintercept = allBare$cmic_mean - allBare$cmic_se) +
  geom_hline(yintercept = allBare$cmic_mean + allBare$cmic_se) +
  xlab("") +
  ylab(expression(paste(C[mic], " (mg C ", g^{-1}, ")")))



#### BUILD PLOTS ---------------------------------------------------------------

postscript(
  file = "./figure_builds/glasshouse_soil_resp.eps",
  width = mm2in(150),
  height = mm2in(60)
)
cowplot::plot_grid(
  ts, pp, kp, cmic,
  nrow = 1, 
  rel_widths = c(2, 1, 1, 1),
  align = "hv"
)
dev.off()

