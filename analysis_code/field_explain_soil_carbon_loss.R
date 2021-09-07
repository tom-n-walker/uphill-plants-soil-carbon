################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Explaining field soil C loss 
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())
# set random seed
set.seed(3)
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
library(vegan)


#### DATA ----------------------------------------------------------------------

## Load from drake plan ----
# main targets
fieldData <- drake::readd(field_data_subset)
fieldTraits <- drake::readd(trait_data)
ghTraits <- drake::readd(gh_plants)
# separate plant datasets
westFocalsLong <- fieldData$focals$all_long[[2]]
centFocalsLong <- fieldData$focals$all_long[[1]]
westFocalsWide <- fieldData$focals$all_wide[[2]]
centFocalsWide <- fieldData$focals$all_wide[[1]]
responses <- fieldData$plantsRR
# remove block for wide datasets (ease of modelling)
westFocalsWide <- westFocalsWide %>%
  select(-block)
centFocalsWide <- centFocalsWide %>%
  select(-block)
# remove Dactylis, Viola, Plantago and Trifolium from central alps (very rare)
onlyThere <- filter(centFocalsLong, cover > 0)
table(onlyThere$species)
centFocalsLong <- centFocalsLong %>%
  filter(species != "Dactylis.glomerata") %>%
  filter(species != "Plantago.lanceolata") %>%
  filter(species != "Trifolium.montanum") %>%
  filter(species != "Viola.hirta")
toRemove <- colnames(centFocalsWide) %in% c("Dactylis.glomerata", "Plantago.lanceolata", "Trifolium.montanum", "Viola.hirta")
centFocalsWide <- centFocalsWide[, !toRemove]


#### BROAD TRAIT DIFFERENCES ---------------------------------------------------

## Field traits ----
# format traits
ftOnly <- fieldTraits %>%
  select(leaf_area:SLA) %>%
  mutate(leaf_area = log10(leaf_area)) %>%
  mutate(plant_height = log10(plant_height)) %>%
  mutate(seed_mass = log10(seed_mass))
# PCA
ftPCA <- prcomp(ftOnly, center = T, scale = T)
# PERMANOVA
adonis2(ftOnly ~ is_focal, fieldTraits, method = "gower")

## Glasshouse traits ----
# format traits
ghtOnly <- ghTraits %>%
  select(rAGB.mg:gsmax)
# PCA
ghtPCA <- prcomp(ghtOnly, center = T, scale = T)
# PERMANOVA
adonis2(ghtOnly ~ Origin, ghTraits)

## Biplots ----
par(mfrow = c(1,1))
biplot(ftPCA)
biplot(ghtPCA)


#### ANALYSE IDENTITY EFFECTS - RANDOM INTERCEPT -------------------------------

## West Alps ----
# build model
m1 <- lme(
  soilCloss ~ cover, 
  random = ~ 1 | species, 
  data = westFocalsLong,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r1 <- residuals(m1, type = "normalized")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
plot(r1 ~ westFocalsLong$cover)
hist(r1)
# test main effects
m1a <- update(m1, ~.- cover)
anova(m1, m1a)

## Central Alps ----
# build model
m2 <- lme(
  soilCloss ~ cover, 
  random = ~ 1 | species, 
  data = centFocalsLong,
  na.action = "na.exclude",
  method = "ML"
)
# diagnose model
r2 <- residuals(m2, type = "normalized")
par(mfrow = c(1, 3))
plot(r2 ~ fitted(m2))
plot(r2 ~ centFocalsLong$cover)
hist(r2)
# test main effects
m2a <- update(m2, ~.- cover)
anova(m2, m2a)


#### ANALYSE IDENTITY EFFECTS - SPECIES SEPARATE -------------------------------

## West Alps ----
# build full and null models
m3x <- lm(soilCloss ~ 1., westFocalsWide)
m3 <- lm(soilCloss ~ ., westFocalsWide)
# diagnose model fit
par(mfrow = c(1, 4))
plot(m3)
# overall model fit
anova(m3x, m3)
# individual main effects
anova(m3)

## Central Alps ----
# build full and null models
m4x <- lm(soilCloss ~ 1., centFocalsWide)
m4 <- lm(soilCloss ~ ., centFocalsWide)
# diagnose model fit
par(mfrow = c(1, 4))
plot(m4)
# overall model fit
anova(m4x, m4)
# individual main effects
anova(m4)


#### ANALYSE TOTAL LOWLAND COVER -----------------------------------------------

# build model
m5 <- lm(
  Csoil ~ focal_bio + bkgnd_bio + vege_bio + site,
  data = responses
)
# diagnose
plot(m5)
# main effects
anova(m5)


#### ANALYSE BACKGROUND COMMUNITY COMPOSITION ----------------------------------

## Analyse treatment effects on raw NMDS 1 scores ----
# build model
m7 <- lme(
  bgPCC2 ~ treatment * site,
  random = ~ 1 | site_block,
  data = fieldData$plantsFull, 
  method = "ML"
)
# diagnose model
r7 <- residuals(m7, type = "normalized")
par(mfrow = c(1, 3))
plot(r7 ~ fitted(m7))
boxplot(r7 ~ fieldData$plantsFull$treatment)
hist(r7)
# test main effects
m7a <- update(m7, ~.- treatment:site)
anova(m7, m7a)
anova(m7a, update(m7a, ~.- treatment))
anova(m7a, update(m7a, ~.- site))

## Analyse NMDS change effects on soil C loss ----
# build model
m8 <- lm(
  Csoil ~ bgPCC1 + bgPCC2 + site,
  data = responses
)
# diagnose
par(mfrow = c(1, 4))
plot(m8)
# main effects
anova(m8)


#### PLOT ----------------------------------------------------------------------


