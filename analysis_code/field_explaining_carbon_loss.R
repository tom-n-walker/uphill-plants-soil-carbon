################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Contextualising soil C loss
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

## Libraries ----
# standard library set
library(tidyverse)
library(nlme)
library(emmeans)

## Load data ----
allData <- drake::readd(field_data_subset)
tryData <- drake::readd(trait_data)
ghPlants <- drake::readd(gh_plants)


#### LOWLAND PLANT TRAITS ------------------------------------------------------

## Do PCA on traits ----
# check normality
par(mfrow = c(1, 3))
hist(tryData$leaf_area)
hist(tryData$plant_height)
hist(tryData$seed_mass)
# transform
traits <- tryData %>%
  mutate(logLA = log10(leaf_area),
         logPH = log10(plant_height),
         logSM = log10(seed_mass)) %>%
  select(logPH, logSM, logLA, SLA, leaf_C:leaf_N)
# do PCA
traitPCA <- prcomp(traits, center = T, scale = T)
summary(traitPCA)
# add scores to data
traitScores <- bind_cols(
  select(tryData, accepted_name, is_focal),
  as.data.frame(traitPCA$x[, 1:3])
)

## Model alpine/lowland effects on traits ----
# permanova
adonis2(traits ~ is_focal, tryData, method = "gower")
# trait scores
anova(lm(PC1 ~ is_focal, traitScores))
anova(lm(PC2 ~ is_focal, traitScores))
anova(lm(PC3 ~ is_focal, traitScores))

## Rarefy alpine trait data to test robustness ----
# build data sets
focTraits <- traits %>%
  filter(tryData$is_focal == "yes") %>%
  mutate(type = "focal")
bgdTraits <- traits %>%
  filter(tryData$is_focal == "no") %>%
  apply(2, function(x) rnorm(n = 16, mean = mean(x), sd = sd(x))) %>%
  as.data.frame %>%
  mutate(type = "alpine")
rareTraits <- rbind(focTraits, bgdTraits)
# do permanova
adonis2(select(rareTraits, -type) ~ type, rareTraits, method = "gower")

## Plot ----
# create segments and centroids
cent <- aggregate(cbind(PC1, PC2) ~ is_focal, 
                  data = traitScores, FUN = mean)
segm <- merge(traitScores, setNames(cent, c("is_focal", "oPC1", "oPC2")),
              by = "is_focal", sort = F)
# ggbiplot to get arrow information
arrows <- ggplot_build(
  ggbiplot::ggbiplot(
    traitPCA, 
    groups = traitScores$is_focal, 
    circle = T)
)$data[[2]]

# generate hulls
hulls <- traitScores %>%
  group_by(is_focal) %>%
  slice(chull(PC1, PC2))

# plot biplot
pca_plot <- ggplot(traitScores) +
  theme_bw() + theme(panel.grid = element_blank()) +
  guides(col = "none", fill = "none") +
  scale_color_manual(values = c("#9eacb9", "#dab59e")) +
  scale_fill_manual(values = c("#9eacb9", "#dab59e")) +
  aes(x = PC1, y = PC2, fill = is_focal, col = is_focal) +
  # geom_polygon(data = hulls, alpha = 0.2) +
  geom_segment(data = segm, aes(xend = oPC1, yend = oPC2), size = 0.1) +
  geom_point(shape = 21, size = 0.25) +
  geom_point(data = cent, shape = 21, size = 2.5, fill = "white") +
  geom_segment(data = arrows, aes(x = x, xend = xend, 
                                  y = y, yend = yend, 
                                  colour = NULL, fill = NULL), 
               arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "PC2 (33.1% var.)", y = "PC3 (20.1% var.)")




# generate sampled dataset, collate with foc_traits
bgd_sampled <- rbind(foc_traits,
                     apply(bgd_traits, 2, rarefy, n = 16)) %>%
  as.data.frame %>%
  mutate(type = rep(c("foc", "bgd"), each = 16))

# do pca on sampled dataset
samp_pca <- prcomp(bgd_sampled %>% select(-type), scale = T, center = T)
samp_scores <- data.frame(type = bgd_sampled$type, 
                          PC1 = samp_pca$x[, 1],
                          PC2 = samp_pca$x[, 2])

# analysis on these scores
anova(lm(PC1 ~ type, samp_scores))
anova(lm(PC2 ~ type, samp_scores))
adonis2(bgd_sampled %>% select(-type) ~ type, bgd_sampled, method = "gower")



#### SPLIT TRAITS --------------------------------------------------------------

## Do PCA on traits ----
pcaTraits <- prcomp(select(allData$plantsFull, leaf_area:SLA), scale = T, center = T)
summary(pcaTraits)
pcaTraits$rotation

par(mfrow = c(3, 2))
boxplot(pcaTraits$x[, 1] ~ allData$plantsFull$site)
boxplot(pcaTraits$x[, 1] ~ allData$plantsFull$treatment)
boxplot(pcaTraits$x[, 2] ~ allData$plantsFull$site)
boxplot(pcaTraits$x[, 2] ~ allData$plantsFull$treatment)
boxplot(pcaTraits$x[, 3] ~ allData$plantsFull$site)
boxplot(pcaTraits$x[, 3] ~ allData$plantsFull$treatment)

## Background community ----
m1 <- lme(
  fixed = bgPCC2 ~ site * treatment, 
  random = ~ 1 | site_block, 
  method = "ML",
  na.action = "na.exclude",
  data = allData$plantsFull
)
r1 <- residuals(m1, type = "pearson")
par(mfrow = c(1, 3))
plot(r1 ~ fitted(m1))
boxplot(r1 ~ allData$plantsFull$treatment)
hist(r1)
# main effects
m1a <- update(m1, ~.-site:treatment)
anova(m1, m1a)
anova(m1a, update(m1a, ~.- treatment))
# post-hoc
emmeans(m1, pairwise ~ treatment | site)


#### ALPINE PCC ----------------------------------------------------------------

m1 <- lme(Csoil ~ bgPCC1 + bgPCC2, random = ~ site | block, allData$plantsRR, method = "ML")
m2 <- lme(Csoil ~ site*PC3, random = ~ 1 | block, test, method = "ML")
anova(m2)

anova(m1)





## Basic treatment effects ----
# Build model
m1 <- lme(
  Csoil ~ site * bgPCC1 + site * bgPCC2, 
  random = ~ 1 | site_block, 
  data = allData$plantsFull,
  method = "ML"
)



anova(m1)
