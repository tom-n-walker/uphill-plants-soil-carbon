################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Glasshouse plant effects on soil
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    31 May 2021
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

## Source functions ----
source("./r_code/apply_mice.R")

## Libraries ----
# standard library set
library(nlme)
library(tidyverse)


#### DATA ----------------------------------------------------------------------

## Load from drake cache ----
# all glasshouse data
gh <- readd(gh_soil)
# subset
gh_plants <- gh$pot_plants
gh_soil <- gh$pot_soil
gh_pots <- left_join(gh_plants, gh_soil)


#### PCA of CWM TRAITS ---------------------------------------------------------

# build
pcaPlants <- prcomp(
  apply_mice(select(gh_pots, AGB.g:cGS.umol.m2.s), 5),
  scale = T, center = T
)
# put scores into pot data frame
gh_pots$plantsPC1 <- pcaPlants$x[, 1]
gh_pots$plantsPC2 <- pcaPlants$x[, 2]
# extract loadings
loadings <- pcaPlants$rotation[, 1:2] %>%
  as.data.frame %>%
  rownames_to_column("variable") %>%
  mutate(variable = fct_reorder(variable, PC1, min))


#### ANALYSE -------------------------------------------------------------------

# PCA scores of traits on initial microbial respiration
m1 <- lme(
  intMR.ugC.g.h ~ plantsPC1 + treatment * plantsPC2,
  random = ~ 1 | block,
  method = "ML",
  data = gh_pots,
  na.action = na.exclude
)
# diagnose
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ gh_pots$treatment)
hist(r1)
# main effects
drop1(m1, test = "Chisq")


#### PLOT ----------------------------------------------------------------------

# scores
pc1Plot <- ggplot(gh_pots) +
  aes(x = plantsPC1, y = intMR.ugC.g.h) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none", fill = "none") +
  scale_fill_manual(values = c("#43978D", "#F9AD6A")) +
  scale_color_manual(values = c("#43978D", "#F9AD6A")) +
  geom_point(aes(fill = treatment), shape = 21) +
  geom_smooth(method = "lm", se = F, fullrange = T, col = "black") +
  labs(x = "PC1 score", y = expression(paste("Initial R (µg C ", g^{-1}, " ", h^{-1}, ")")))
pc2Plot <- ggplot(gh_pots) +
  aes(x = plantsPC2, y = intMR.ugC.g.h) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none", fill = "none") +
  scale_fill_manual(values = c("#43978D", "#F9AD6A")) +
  scale_color_manual(values = c("#43978D", "#F9AD6A")) +
  geom_point(aes(fill = treatment), shape = 21) +
  geom_smooth(aes(col = treatment), method = "lm", se = F, fullrange = T) +
  labs(x = "PC2 score", y = expression(paste("Initial R (µg C ", g^{-1}, " ", h^{-1}, ")")))
# loadings
pc1loads <- ggplot(loadings) +
  aes(x = variable, y = PC1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-0.5, 0.7)) +
  geom_segment(aes(yend = 0, xend = variable)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "PC1 loading (49.0% var.)")
pc2loads <- ggplot(loadings) +
  aes(x = variable, y = PC2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-0.5, 0.7)) +
  geom_segment(aes(yend = 0, xend = variable)) +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x = "", y = "PC2 loading (19.3% var.)")
# put together
cowplot::plot_grid(pc1Plot, pc2Plot, pc1loads, pc2loads,
                   nrow = 2, align = "hv")
