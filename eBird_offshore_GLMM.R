# 
# eBird offshore GLMM binary
# Mar 13, 2024
#
# Will Brooks
#
####

setwd("~/Desktop/eBird offshore")

# Clear environment and load packages ----

rm(list = ls()) # clean environment

library(tidyverse)
library(lme4)
library(GGally)
library(AICcmodavg)
library(sjPlot)
library(influence.ME)
library(DHARMa)
library(ncf)
library(stargazer)
library(gridExtra)

set.seed(2000)
theme_set(theme_classic())

# Load and explore data ----

#load data
dat <- read.csv("eBird_covs2023-11-15.csv") %>% 
  subset(species_observed != "NA") %>%
  unique() %>%
  # Removing observations before 1947 because we lack data before that period
  subset(year(ymd(observation_date)) > 1947) %>% 
  # Adding variables for month and year
  mutate(month=month(ymd(observation_date)),
         year=year(ymd(observation_date))) %>% 
  # Removing observations where the species wasn't seen on land and those in winter and summer
  subset(land_obs_freq > 0) %>% 
  subset(as.character(month) %in% c("3","4","5","8","9","10","11")) %>%
  # Mark spring and fall observations
  mutate(season = as.factor(
    if_else(as.character(month) %in% c("3","4","5"),"Spring","Fall"))) %>%
  # Remove species with no observations
  group_by(scientific_name) %>%
  mutate(ever_obs = if_else( sum(species_observed) > 0,TRUE,FALSE)) %>%
  ungroup() %>% 
  subset(ever_obs == "TRUE") %>% 
  # Mark checklists where any bird was seen
  group_by(checklist_id) %>%
  mutate(list_obs = if_else( sum(species_observed) > 0,TRUE,FALSE)) %>%
  ungroup() 

# another I think better approach, just subsampling non-presences for each species?
presences <- subset(dat, species_observed=="TRUE")
absences <- subset(dat, species_observed=="FALSE") 
absences_sub <- absences[sample(nrow(absences),nrow(presences)),]

dat3 <- rbind(presences,absences_sub)

### data exploration
str(dat3)

# Check columns with nas
names(which(colSums(is.na(dat3)) > 0))
# I don't care about any of these variables with NAs

# Check dataset size
table(dat3$species_observed)

# Check distributions and correlations
ggpairs(dat3[, c("uwind", "vwind", "ap_mean","Ra","relh","land_obs_freq","latitude")])

# The distribution of ap looks a little weird
hist(dat3$ap_mean)
# I'm going to leave it for now, clearly skewed but not a massive gap in the data

# Check VIF
dat3 %>% 
  select("uwind", "vwind", "ap_mean", "Ra","relh","land_obs_freq","latitude") %>% 
  as.data.frame() %>%
  usdm::vif()
# Looks great

# Lastly I want to look at the relationship between each predictor and presence and absence. I'll do a series of boxplots

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = uwind)) +
  geom_boxplot()

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = vwind)) +
  geom_boxplot()

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = ap_mean)) +
  geom_boxplot()

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = Ra)) +
  geom_boxplot()

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = relh)) +
  geom_boxplot()

ggplot(data = dat3, aes(x = factor(species_observed), 
                        y = land_obs_freq)) +
  geom_boxplot()
# Not seeing much here

# Lastly lets look at the observations for each possible random effect
table(dat3$scientific_name) 

# Too many dates so I have to do boxplot
boxplot(table(dat3$observation_date))
# Quite the range, seems appropriate for a random effect

# Testing random effect structure ----

# I feel confident I should use an intercept for species, but otherwise I want to test some different structures.

modList.ranef <- list()

modList.ranef[["Species Only"]] <- glmm1 <- 
  glmer(species_observed ~ 
          scale(uwind) + scale(vwind) + season + scale(ap_mean) + scale(Ra) + 
          (1 | scientific_name), 
        data = dat3, 
        family = binomial)

modList.ranef[["Species + Date"]] <- glmm2 <- 
  glmer(species_observed ~ 
          scale(uwind) + scale(vwind) + season + scale(ap_mean) + scale(Ra) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList.ranef[["Species + checklist_id"]] <- glmm3 <- 
  glmer(species_observed ~ 
          scale(uwind) + scale(vwind) + season + scale(ap_mean) + scale(Ra) + 
          (1 | scientific_name) + (1 | checklist_id), 
        data = dat3, 
        family = binomial)

modList.ranef[["Species + Date/checklist_id"]] <- glmm4 <- 
  glmer(species_observed ~ 
          scale(uwind) + scale(vwind) + season + scale(ap_mean) + scale(Ra) + 
          (1 | scientific_name) + (1 | observation_date/checklist_id), 
        data = dat3, 
        family = binomial)
# Failed to converge

modList.ranef[["Species + observer_id"]] <- glmm5 <- 
  glmer(species_observed ~ 
          scale(uwind) + scale(vwind) + season + scale(ap_mean) + scale(Ra) + 
          (1 | scientific_name) + (1 | observer_id), 
        data = dat3, 
        family = binomial)

# compare ranef structures
aictab.ranef <- aictab(modList.ranef)

# Looks like the best model that converged is the observation date ranef

stargazer(aictab.ranef, 
          title = "", 
          summary = F,
          type = "html",
          digits = 2,
          font.size = "normalsize",
          out = "myAICtableranef.doc")

# Model selection ----

modList <- list()

modList[["Intercept only"]] <- glmm6 <- 
  glmer(species_observed ~ 1 + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq"]] <- glmm7 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["relh"]] <- glmm8 <- 
  glmer(species_observed ~ 
          scale(relh) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["uwind"]] <- glmm9 <- 
  glmer(species_observed ~ scale(uwind) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["vwind"]] <- glmm10 <- 
  glmer(species_observed ~ scale(vwind) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["season"]] <- glmm11 <- 
  glmer(species_observed ~ season + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["ap_mean"]] <- glmm12 <- 
  glmer(species_observed ~ scale(ap_mean) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["Ra"]] <- glmm13 <- 
  glmer(species_observed ~ scale(Ra) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

# Now I will build up the parts of the model that accounts for amount of migration overall. This is on land frequency, vwind (north/south), and the interaction between vwind and season. Of course land frequency is important. North and south winds encourage migrants to move, depending on the season. These will be included in all future models, which will aim to test hypotheses for offshore vagrancy.

modList[["land_obs_freq + vwind"]] <- glmm14 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + season"]] <- glmm15 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + season + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind + season"]] <- glmm16 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) + season + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind * season"]] <- glmm17 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

aictab(modList)
# Indeed the land_obs_freq + vwind * season was the best fit for the data. I will now try variables possibly explaining vagrancy.

modList[["land_obs_freq + vwind * season + uwind"]] <- glmm18 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(uwind) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind * season + ap_mean"]] <- glmm19 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(ap_mean) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind * season + Ra"]] <- glmm20 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(Ra) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind * season + relh"]] <- glmm21 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(relh) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

# modList[["land_obs_freq + vwind * season + uwind + ap_mean"]] <- glmm22 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(ap_mean) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + uwind + Ra"]] <- glmm23 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(Ra) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge

modList[["land_obs_freq + vwind * season + uwind + relh"]] <- glmm24 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(relh) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

# modList[["land_obs_freq + vwind * season + ap_mean + Ra"]] <- glmm25 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(ap_mean) + scale(Ra) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + ap_mean + relh"]] <- glmm26 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(ap_mean) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + Ra + relh"]] <- glmm27 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(Ra) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + uwind + ap_mean + Ra"]] <- glmm28 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(ap_mean) + scale(Ra) + 
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + Ra + ap_mean + relh"]] <- glmm29 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(Ra) + scale(ap_mean) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + uwind + Ra + relh"]] <- glmm30 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(Ra) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + uwind + ap_mean + relh"]] <- glmm31 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(ap_mean) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge
# 
# modList[["land_obs_freq + vwind * season + uwind + ap_mean + Ra + relh"]] <- glmm32 <- 
#   glmer(species_observed ~ 
#           scale(land_obs_freq) + scale(vwind) * season + scale(uwind) + scale(ap_mean) + scale(Ra) + scale(relh) +
#           (1 | scientific_name) + (1 | observation_date), 
#         data = dat3, 
#         family = binomial)
# # Failed to converge

# Adding a few models I forgot before
modList[["land_obs_freq + uwind * vwind"]] <- glmm33 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * scale(uwind) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + uwind * vwind + season"]] <- glmm34 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * scale(uwind) + season +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["latitude"]] <- glmm35 <- 
  glmer(species_observed ~ 
          scale(latitude) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)

modList[["land_obs_freq + vwind * season + latitude"]] <- glmm36 <- 
  glmer(species_observed ~ 
          scale(land_obs_freq) + scale(vwind) * season + scale(latitude) +
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)


# Compare models
aictab <- aictab(modList)
aictab

# Model agerage
avg.mod <- MuMIn::model.avg(modList)
plot(avg.mod, intercept = FALSE)
avg.mod.coef <- plot(avg.mod, intercept = FALSE) %>% 
  as.data.frame( )
rownames(avg.mod.coef) <- c("land obs freq", "vwind", "season [Spring]", "relh", "vwind x season [Spring]", "latitude", "Ra", "ap_mean", "uwind", "uwind x vwind")

stargazer(aictab, 
          title = "", 
          summary = F,
          type = "html",
          digits = 2,
          font.size = "normalsize",
          out = "myAICtable.doc")

# Extra models ----

glmm.lat <- 
  glmer(species_observed ~ scale(latitude) + 
          (1 | scientific_name) + (1 | observation_date), 
        data = dat3, 
        family = binomial)
summary(glmm.lat)
# No relationship

# Present results of best model ----

bestMod <- glmm20
summary(bestMod)

car::vif(bestMod) 
# No issues with multicollinearity

# Rsquared
performance::r2(bestMod)
# While the model is an adequate fit the fixed effects explain little of the variance

# Plot fixed effects
effects <- plot_model(bestMod,
           type = "est")
(effects <-  effects + 
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + 
  labs(title = "D") +
  geom_hline(aes(yintercept = 1),colour = "gray",alpha = 0.1))
  
(effects.avg <- avg.mod.coef %>%
  ggplot(aes(x = px, y = rownames(.))) +
  geom_vline(aes(xintercept = 0),colour = "gray") +
  geom_errorbarh(aes(xmin = lx0, xmax = lx1), 
                 height = 0, 
                 size = 0.7,
                 col = if_else(avg.mod.coef$px > 0,"#3f7fb4ff","Red2")) +
  geom_point(size = 2.5, 
             col = if_else(avg.mod.coef$px > 0,"#3f7fb4ff","Red2")) +
  labs(title = "E",x = "Full model averaged estimate",
       y = "") +
  theme_classic() +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")))


# Create nice model results table
tab_model(bestMod, 
          show.se = TRUE, 
          show.stat = TRUE, 
          show.icc = FALSE,
          transform = NULL,
          show.p = TRUE,
          string.est = "Beta",
          string.se = "se",
          use.viewer = TRUE,
          show.re.var = FALSE,
          show.ngroups = FALSE,
          dv.labels = c(""))

# Plotting ----

# Plotting for each important variable

# Ra, make new data
new_data_Ra <- expand.grid(Ra = seq(min(dat3$Ra), 
                                                      max(dat3$Ra), 
                                                      length.out = 1000),
                              vwind = median(dat3$vwind),
                              land_obs_freq = median(dat3$land_obs_freq),
                              season = levels(dat3$season))

# Predict new data
pred.glmm.Ra <- data.frame(new_data_Ra, 
                                    predictSE.merMod(bestMod, 
                                                     newdata = new_data_Ra, 
                                                     type = "link")) 

pred.glmm.Ra <- pred.glmm.Ra %>%
  mutate(fitReal = plogis(fit),
         lci = plogis(fit - 1.96 * se.fit),
         uci = plogis(fit + 1.96 * se.fit) # plogis because binomial
  )

#vwind, make new data
new_data_vwind <- expand.grid(vwind = seq(min(dat3$vwind), 
                                          max(dat3$vwind), 
                                          length.out = 1000),
                              Ra = median(dat3$Ra),
                              land_obs_freq = median(dat3$land_obs_freq),
                              season = levels(dat3$season))

# predict new data
pred.glmm.vwind <- data.frame(new_data_vwind, 
                              predictSE.merMod(bestMod, 
                                               newdata = new_data_vwind, 
                                               type = "link")) 

pred.glmm.vwind <- pred.glmm.vwind %>%
  mutate(fitReal = plogis(fit),
         lci = plogis(fit - 1.96 * se.fit),
         uci = plogis(fit + 1.96 * se.fit) # plogis because binomial
  )


# land_obs_freq, make new data
new_data_land_obs_freq <- expand.grid(land_obs_freq	 = seq(min(dat3$land_obs_freq), 
                                          max(dat3$land_obs_freq), 
                                          length.out = 1000),
                              vwind = median(dat3$vwind),
                              Ra = median(dat3$Ra),
                              season = levels(dat3$season))

# predict new data
pred.glmm.land_obs_freq <- data.frame(new_data_land_obs_freq, 
                              predictSE.merMod(bestMod, 
                                               newdata = new_data_land_obs_freq, 
                                               type = "link")) 

pred.glmm.land_obs_freq <- pred.glmm.land_obs_freq %>%
  mutate(fitReal = plogis(fit),
         lci = plogis(fit - 1.96 * se.fit),
         uci = plogis(fit + 1.96 * se.fit) # plogis because binomial
  )


# Making binned frequencies to see model fit

# Ra binned frequencies
testData.Ra <- dat3 %>% 
  mutate(cutVol = 
           cut(.$Ra, 6)) %>% # Cut into six bins
  group_by(cutVol) %>%
  summarize(species_observed = sum(as.numeric(species_observed) == "1"), # Count presences in each bin
            count = n(), # Total counts in bin
            prop = species_observed/count, # Proportion of presence
            Ra_Means = mean(Ra), # Mean of variable in bin
            se = sqrt(prop*(1-prop)/count)) 

# vwind binned frequencies
testData.vwind.fall <- dat3 %>% 
  subset(season == "Fall") %>%
  mutate(cutVol = 
           cut(.$vwind, 6)) %>% # Cut into six bins
  group_by(cutVol) %>%
  summarize(species_observed = sum(as.numeric(species_observed) == "1"), # Count presences in each bin
            count = n(), # Total counts in bin
            prop = species_observed/count, # Proportion of presence
            vwind_Means = mean(vwind), # Mean of variable in bin
            se = sqrt(prop*(1-prop)/count)) 

testData.vwind.spring <- dat3 %>% 
  subset(season == "Spring") %>%
  mutate(cutVol = 
           cut(.$vwind, 6)) %>% # Cut into six bins
  group_by(cutVol) %>%
  summarize(species_observed = sum(as.numeric(species_observed) == "1"), # Count presences in each bin
            count = n(), # Total counts in bin
            prop = species_observed/count, # Proportion of presence
            vwind_Means = mean(vwind), # Mean of variable in bin
            se = sqrt(prop*(1-prop)/count)) 

# Ra binned frequencies
testData.land_obs_freq <- dat3 %>% 
  mutate(cutVol = 
           cut(.$land_obs_freq, 6)) %>% # Cut into six bins
  group_by(cutVol) %>%
  summarize(species_observed = sum(as.numeric(species_observed) == "1"), # Count presences in each bin
            count = n(), # Total counts in bin
            prop = species_observed/count, # Proportion of presence
            land_obs_freq_Means = mean(land_obs_freq), # Mean of variable in bin
            se = sqrt(prop*(1-prop)/count)) 

# Make plots

# Ra

# Merging the dataframes by the common column 'ID'
pred.glmm.Ra.no.season <- merge(subset(pred.glmm.Ra,season == "Fall"), 
                       subset(pred.glmm.Ra,season == "Spring"), by = "Ra")

# Calculating the average values
pred.glmm.Ra.no.season$fitReal <- rowMeans(pred.glmm.Ra.no.season[, c("fitReal.x", "fitReal.y")])
pred.glmm.Ra.no.season$lci <- rowMeans(pred.glmm.Ra.no.season[, c("lci.x", "lci.y")])
pred.glmm.Ra.no.season$uci <- rowMeans(pred.glmm.Ra.no.season[, c("uci.x", "uci.y")])

plot.Ra <- ggplot(data = pred.glmm.Ra.no.season, 
                           mapping = aes(x = Ra, 
                                         y = fitReal * nrow(dat3)/nrow(dat) )) +
  geom_line() +
  geom_ribbon(aes(ymin = lci * nrow(dat3)/nrow(dat) , 
                  ymax = uci * nrow(dat3)/nrow(dat) ),
              alpha = 0.25) +
  labs(x = "Solar radiation (relative sunspot abundance)", 
       y = "Probability of offshore passerine",
       title = "B") +
  geom_segment(data = testData.Ra, 
               aes(x = Ra_Means, 
                   y = (prop - se)  * nrow(dat3)/nrow(dat) , 
                   xend = Ra_Means, 
                   yend = (prop + se) * nrow(dat3)/nrow(dat) )) + # Binned frequency SE
geom_point(data = testData.Ra, 
           aes(x = Ra_Means, 
               y = prop * nrow(dat3)/nrow(dat) ), 
           pch = 16, 
           cex = 2) # Add binned frequencies
plot.Ra

# vwind
plot.vwind <- ggplot(data = pred.glmm.vwind, 
                     mapping = aes(x = vwind, 
                                   y = fitReal * nrow(dat3)/nrow(dat) )) +
  geom_line(aes(linetype = season)) +
  geom_ribbon(aes(ymin = lci * nrow(dat3)/nrow(dat) , 
                  ymax = uci * nrow(dat3)/nrow(dat) , 
                  fill = season),
              alpha = 0.25) +
  labs(x = "v-wind (m/s)", 
       y = "Probability of offshore passerine",
       title = "C") +
  geom_segment(data = testData.vwind.fall, 
               aes(x = vwind_Means, 
                   y = (prop - se) * nrow(dat3)/nrow(dat) , 
                   xend = vwind_Means, 
                   yend = (prop + se) * nrow(dat3)/nrow(dat) )) + # Binned frequency SE 
  geom_point(data = testData.vwind.fall, 
             aes(x = vwind_Means, 
                 y = prop * nrow(dat3)/nrow(dat) ), 
             pch = 16, 
             cex = 2,
             col = "red") + # Add binned frequencies
  geom_segment(data = testData.vwind.spring, 
               aes(x = vwind_Means, 
                   y = (prop - se) * nrow(dat3)/nrow(dat) , 
                   xend = vwind_Means, 
                   yend = (prop + se) * nrow(dat3)/nrow(dat) )) + # Binned frequency SE 
geom_point(data = testData.vwind.spring, 
           aes(x = vwind_Means, 
               y = prop * nrow(dat3)/nrow(dat) ), 
           pch = 16, 
           cex = 2,
           col = "blue") + # Add binned frequencies
  annotate("text", x=10, y=-.0008, label= "South wind",color = "red") +
  annotate("text", x=-15, y=-.0008, label= "North wind",color = "red") +
  coord_cartesian(ylim = c(0, 0.004), clip = "off")
plot.vwind

# on land obs freq

# Merging the dataframes by the common column 'ID'
pred.glmm.land_obs_freq.no.season <- merge(subset(pred.glmm.land_obs_freq,season == "Fall"), 
                                subset(pred.glmm.land_obs_freq,season == "Spring"), by = "land_obs_freq")

# Calculating the average values
pred.glmm.land_obs_freq.no.season$fitReal <- rowMeans(pred.glmm.land_obs_freq.no.season[, c("fitReal.x", "fitReal.y")])
pred.glmm.land_obs_freq.no.season$lci <- rowMeans(pred.glmm.land_obs_freq.no.season[, c("lci.x", "lci.y")])
pred.glmm.land_obs_freq.no.season$uci <- rowMeans(pred.glmm.land_obs_freq.no.season[, c("uci.x", "uci.y")])


plot.land_obs_freq <- ggplot(data = pred.glmm.land_obs_freq.no.season, 
                     mapping = aes(x = land_obs_freq, 
                                   y = fitReal * nrow(dat3)/nrow(dat) )) +
  geom_line() +
  geom_ribbon(aes(ymin = lci * nrow(dat3)/nrow(dat) , 
                  ymax = uci * nrow(dat3)/nrow(dat)),
              alpha = 0.25) +
  labs(x = "On-land observation frequency", 
       y = "Probability of offshore passerine",
       title = "A") +
  geom_segment(data = testData.land_obs_freq, 
               aes(x = land_obs_freq_Means, 
                   y = (prop - se) * nrow(dat3)/nrow(dat) , 
                   xend = land_obs_freq_Means, 
                   yend = (prop + se) * nrow(dat3)/nrow(dat) )) + # Binned frequency SE
  geom_point(data = testData.land_obs_freq, 
             aes(x = land_obs_freq_Means, 
                 y = prop * nrow(dat3)/nrow(dat) ), 
             pch = 16, 
             cex = 2) # Add binned frequencies
plot.land_obs_freq

# Plot the two variables side by side
binomial_predictions <- arrangeGrob(plot.land_obs_freq,plot.Ra,plot.vwind,effects,effects.avg, ncol = 3)

ggsave("Binomial_predictions.jpeg",
       plot = binomial_predictions,
       device = "jpeg",
       width = 3500,
       height = 2000,
       units = "px")

# Model diagnostics: residuals ----

mean(residuals(bestMod)) # mean is pretty close to zero

# Plotting binned fitted values against residuals
source("~/Desktop/CONS625/functions/diagnosticPlotBinomialGLM.r")
par(mfrow=c(1,1))

plotDiagBinom(bestMod, 
              xvar = fitted.values(bestMod), 
              xlab ="fitted values", 
              main = "binned residual means", 
              bin.size = 100) 
# there is some non-linearity here

# Plot against vwind predictor
plotDiagBinom(bestMod, 
              xvar = dat3$vwind, 
              xlab ="vwind", 
              main = "binned residual means", 
              bin.size = 25) 

# Plot against land_obs_freq predictor
plotDiagBinom(bestMod, 
              xvar = dat3$land_obs_freq, 
              xlab ="land_obs_freq", 
              main = "binned residual means", 
              bin.size = 25) 

# Plot against Ra predictor
plotDiagBinom(bestMod, 
              xvar = dat3$Ra, 
              xlab ="Ra", 
              main = "binned residual means", 
              bin.size = 25) 

# Simulate residuals with DHARMa

res <- DHARMa::simulateResiduals(bestMod, 
                                 plot = T) 
# Again slight deviance from normality

DHARMa::plotResiduals(res, 
                      form = dat3$uwind) 

DHARMa::plotResiduals(res, 
                      form = dat3$land_obs_freq) 

DHARMa::plotResiduals(res, 
                      form = dat3$Ra) 

# Again, some non-linearity of residuals

# Check in sjplot

plot_model(bestMod, type = "diag")

# Lets see what the root of some of these irregularities is. First, I'll compare model residuals to unused predictors to see if I can find any problems

# number of ap_mean
plotDiagBinom(bestMod, 
              xvar = dat3$ap_mean, 
              xlab ="ap_mean", 
              main = "binned residual means", 
              bin.size = 25) 

# number of uwind
plotDiagBinom(bestMod, 
              xvar = dat3$uwind, 
              xlab ="uwind", 
              main = "binned residual means", 
              bin.size = 25) 

# latitude
plotDiagBinom(bestMod, 
              xvar = dat3$latitude, 
              xlab ="latitude", 
              main = "binned residual means", 
              bin.size = 25) 

# longitude
plotDiagBinom(bestMod, 
              xvar = dat3$longitude, 
              xlab ="longitude", 
              main = "binned residual means", 
              bin.size = 25) 

# year
plotDiagBinom(bestMod, 
              xvar = dat3$year, 
              xlab ="year", 
              main = "binned residual means", 
              bin.size = 25) 

# month
plotDiagBinom(bestMod, 
              xvar = dat3$month, 
              xlab ="month", 
              main = "binned residual means", 
              bin.size = 25) 

# distance offshore
plotDiagBinom(bestMod, 
              xvar = dat3$dist_offshore, 
              xlab ="distance offshore", 
              main = "binned residual means", 
              bin.size = 25) 

# checklist duration
plotDiagBinom(bestMod, 
              xvar = dat3$duration_minutes, 
              xlab ="checklist duration", 
              main = "binned residual means", 
              bin.size = 25) 
# Possible trend here

# number of observers
plotDiagBinom(bestMod, 
              xvar = dat3$number_observers, 
              xlab ="fitted values", 
              main = "binned residual means", 
              bin.size = 25) 
# Nothing with any of these

# The irregularities in residuals could also be due to spatial autocorrelation now

dat3 %>% 
  mutate(resid = residuals(bestMod),
         resid_abs = abs(residuals(bestMod)),
         sign = if_else(resid >= 0, "pos",
                        "neg")) %>% 
  ggplot(aes(x = longitude,
             y = latitude,
             size = abs(resid),
             col = sign)) + 
  geom_point()
# This looks fairly random. This is telling me that there is no spatial autocorrelation.

# Model residuals: check influential points ----

# Cooks distance of individual data points
cd <- cooks.distance(bestMod)

# plot it as vertical bars
plot(cd, type = "h", ylab = "Cook's distance")
abline(h = 10,  lty = 4, col = "red")

cd.idx <- cd >  10
cd.vidx <- cd >  1000

sum(cd.idx) # 126 influential points
# There also look to be 5 very influential points that peak over 100000, I'll test those out too

very_inf <- dat3 %>% filter(cd > 1000)

# Rerun without influential points

bestMod_noinf <- update(bestMod,data = dat3 %>% 
                          filter(cd < 10))

bestMod_noinf2 <- update(bestMod,data = dat3 %>% 
                          filter(cd < 1000))

fixef(bestMod)
fixef(bestMod_noinf) 
fixef(bestMod_noinf2) 
# In all cases the model estimates are similar

# Lets check on the plot for uwind
ggplot(data = dat3, 
       aes(x = uwind, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3[cd.idx,], 
             aes(x = uwind, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for vwind
ggplot(data = dat3, 
       aes(x = vwind, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3[cd.idx,], 
             aes(x = vwind, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for Ra
ggplot(data = dat3, 
       aes(x = Ra, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3[cd.idx,], 
             aes(x = Ra, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for dist_offshore
ggplot(data = dat3, 
       aes(x = dist_offshore, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3[cd.idx,], 
             aes(x = dist_offshore, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for land_obs_freq
ggplot(data = dat3, 
       aes(x = land_obs_freq, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3[cd.idx,], 
             aes(x = land_obs_freq, 
                 y = species_observed), 
             col = "red", 
             size = 1) 
# I don't see anything here

# Check influential days with influence.ME

infRE <- influence.ME::influence(bestMod, group = "observation_date") 

# calculate cook's distance by random effect group
cd2 <- cooks.distance(infRE,
                      sort = TRUE)

cd2.idx <- cd2 > 3 * mean(cd2)
sum(cd2.idx) # 62 influential sites

# Plot and sort these values to identify those above the threshold
plot(infRE, 
     which = "cook",
     cutoff = 3 * mean(cd2), 
     sort = TRUE,
     xlab = "Cook's Distance",
     ylab = "date")

# Lets check on the plot for uwind
ggplot(data = dat3, 
       aes(x = uwind, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3 %>% 
               filter(observation_date %in% 
                        names(cd2.idx[cd2.idx == TRUE,])), 
             aes(x = uwind, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for vwind
ggplot(data = dat3, 
       aes(x = vwind, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3 %>% 
               filter(observation_date %in% 
                        names(cd2.idx[cd2.idx == TRUE,])), 
             aes(x = vwind, 
                 y = species_observed), 
             col = "red", 
             size = 1) 

# Lets check on the plot for Ra
ggplot(data = dat3, 
       aes(x = Ra, 
           y = species_observed)) +
  geom_point() +
  geom_point(data = dat3 %>% 
               filter(observation_date %in% 
                        names(cd2.idx[cd2.idx == TRUE,])), 
             aes(x = Ra, 
                 y = species_observed), 
             col = "red", 
             size = 1) 
# I don't see anything here

# Lets rerun the model without those points included.
bestMod_noinf2 <- update(bestMod,data = dat3 %>% 
                           filter(observation_date %in% 
                                    names(cd2.idx[cd2.idx == FALSE,])))

fixef(bestMod)
fixef(bestMod_noinf2) 
# Nothing massive changes

