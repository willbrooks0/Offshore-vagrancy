#
# Offshore eBird species predictors 
#
# Will Brooks
# Feb 22, 2023
#
#####################

library(car)
library(ggplot2)
library(GGally)
library(tidyverse)
library(ggpubr)
library(AICcmodavg)
library(interactions)
library(stargazer)

setwd("~/Desktop/eBird offshore")

#load data
eBird_offshore_species_info <- read.csv(file = "eBird_offshore_species_info.csv",sep=",")

#load extracted eBird data

eBird_covs <- read.csv("eBird_covs2023-11-15.csv") %>% #load extracted data
  subset(land_obs_freq > 0) %>% 
  mutate(weight = 
           if_else(species_observed == "TRUE", 1 - land_obs_freq, land_obs_freq)) %>% 
  subset(land_obs_freq>0) %>%
  mutate(month=month(ymd(observation_date)),
         year=year(ymd(observation_date))) %>% 
  subset(as.character(month) %in% c("3","4","5","8","9","10","11")) %>%
  # Mark spring and fall observations
  mutate(season = as.factor(
    if_else(as.character(month) %in% c("3","4","5"),"Spring","Fall"))) 

# create offshore vagrancy propensity as a weighted average by on land frequency
offshore_freq <- eBird_covs %>% 
  group_by(scientific_name) %>% 
  reframe(freq_adj = weighted.mean(species_observed, weight),
          freq_raw = mean(species_observed)) %>%
  merge(eBird_offshore_species_info,by="scientific_name") #combine with species info

#stacked barplot of offshore frequency for each species
season2 <- eBird_covs %>%  
  group_by(scientific_name, season) %>% 
  summarise(sightings = sum(species_observed)) %>% 
  group_by(scientific_name) %>% 
  mutate(total=sum(sightings)) %>% ungroup() %>%
  mutate(percent=sightings/total)

stacked_data <- merge(season2,offshore_freq,by="scientific_name") %>% 
  mutate(season_freq = percent * freq_adj,
         season_freq_raw = percent * freq_raw)

plot.freqadj <- ggplot(stacked_data, aes(fill=season, y=season_freq, x=Common_name)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + 
  labs(x="Species",y="Adjusted observation frequency") + 
  coord_flip() + 
  scale_x_discrete(limits=rev) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(text = element_text(size=20),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(0.1,5,0.1,0.1, "cm"))

plot.freqraw <- ggplot(stacked_data, aes(fill=season, y=season_freq_raw, x=Common_name)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + 
  labs(x="Species",y="Raw observation frequency") + 
  coord_flip() + 
  scale_x_discrete(limits=rev) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(text = element_text(size=20),
        legend.position = "none",
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))

species_freq <- arrangeGrob(plot.freqraw, plot.freqadj, ncol = 2)

  ggsave("species_freq.jpeg",
       plot = species_freq,
       device = "jpeg",
       width = 5000,
       height = 4000,
       units = "px")


# compare to morphology ----
avonet <- read.csv("AVONET2.csv")
elton <- read.csv("EltonTraits_guilds.csv")
traits <- avonet %>% subset(Species1 %in% offshore_freq$Avonet)%>%
  group_by(Species1) %>%
  summarize(mass = mean(Mass,na.rm=TRUE),
            tarsus = mean(Tarsus.Length,na.rm=TRUE),
            handwing = mean(Hand.Wing.Index,na.rm=TRUE)) %>% 
  ungroup() %>%
  right_join(offshore_freq,by=c("Species1"="Avonet")) %>%
  left_join(elton,c("Common_name" = "English"))

# explore data
ggpairs(traits[, c("freq_adj", "handwing","Migration_dist",
                   "mass","ForStrat.understory","ForStrat.midhigh",
                   "ForStrat.canopy","ForStrat.aerial")])

par(mfrow = c(1,1))
boxplot(freq_adj ~ Noctural, data = traits)

traits2 <- traits %>%
  # Transforming the response variable with log 10, the lowest non-zero value divided by two was added to prevent undefined values.
  mutate(freq_adj = log10(freq_adj+sort(freq_adj)[12]/2))

ggpairs(traits2[, c("freq_adj", "handwing","Migration_dist",
                   "mass","ForStrat.ground","ForStrat.understory","ForStrat.midhigh",
                   "ForStrat.canopy","ForStrat.aerial")])

# Check VIF
traits2[, c("freq_adj", "handwing","Migration_dist",
            "mass","ForStrat.ground","ForStrat.understory","ForStrat.midhigh",
            "ForStrat.canopy","ForStrat.aerial")] %>%
  as.data.frame() %>%
  usdm::vif()
# really high, removing ground

traits2[, c("freq_adj", "handwing","Migration_dist",
            "mass","ForStrat.understory","ForStrat.midhigh",
            "ForStrat.canopy","ForStrat.aerial")] %>%
  as.data.frame() %>%
  usdm::vif()
# Looks great

hist(traits2$Migration_dist)
# I'm going to transform this one

#setting this for plotting
traits2$Hand_wing <- traits2$handwing

# multiple linear regression ----

modList.sp <- list()

modList.sp[["Null"]] <- m1 <- 
  lm(freq_adj ~ 1, data = traits2)
modList.sp[["Migration_dist"]] <- m2 <- 
  lm(freq_adj ~ Migration_dist, data = traits2)
modList.sp[["Hand_wing"]] <- m3 <- 
  lm(freq_adj ~ Hand_wing, data = traits2)
modList.sp[["Hand_wing + Migration_dist"]] <- m4 <- 
  lm(freq_adj ~ Hand_wing + Migration_dist, data = traits2)
modList.sp[["Hand_wing * Migration_dist"]] <- m5 <- 
  lm(freq_adj ~ Hand_wing * Migration_dist, data = traits2)

modList.sp[["mass"]] <- m6 <- 
  lm(freq_adj ~ mass, data = traits2)
modList.sp[["ForStrat.canopy"]] <- m7 <- 
  lm(freq_adj ~ ForStrat.canopy, data = traits2)
modList.sp[["ForStrat.understory"]] <- m8 <- 
  lm(freq_adj ~ ForStrat.understory, data = traits2)
modList.sp[["ForStrat.midhigh"]] <- m9 <- 
  lm(freq_adj ~ ForStrat.midhigh, data = traits2)
modList.sp[["ForStrat.aerial"]] <- m10 <- 
  lm(freq_adj ~ ForStrat.aerial, data = traits2)

modList.sp[["mass + Hand_wing * Migration_dist"]] <- m11 <- 
  lm(freq_adj ~ mass + Hand_wing * Migration_dist, data = traits2)
modList.sp[["ForStrat.understory + Hand_wing * Migration_dist"]] <- m12 <- 
  lm(freq_adj ~ ForStrat.understory + Hand_wing * Migration_dist, data = traits2)
modList.sp[["ForStrat.canopy + Hand_wing * Migration_dist"]] <- m13 <- 
  lm(freq_adj ~ ForStrat.canopy + Hand_wing * Migration_dist, data = traits2)
modList.sp[["ForStrat.midhigh + Hand_wing * Migration_dist"]] <- m14 <- 
  lm(freq_adj ~ ForStrat.midhigh + Hand_wing * Migration_dist, data = traits2)
modList.sp[["ForStrat.aerial + Hand_wing * Migration_dist"]] <- m15 <- 
  lm(freq_adj ~ ForStrat.aerial + Hand_wing * Migration_dist, data = traits2)

aictab <- aictab(modList.sp)
aictab

stargazer(aictab, 
          title = "", 
          summary = F,
          type = "html",
          digits = 2,
          font.size = "normalsize",
          out = "myAICtable_species.doc")

bestMod <- m12

summary(bestMod)

# Ploting ----


# Generate new data for prediction

# ForStrat.understory

newData <- expand.grid(ForStrat.understory = seq(min(traits2$ForStrat.understory), 
                                      max (traits2$ForStrat.understory), 
                                      length.out = 1000), 
                       Hand_wing = mean(traits2$Hand_wing),
                       Migration_dist = mean(traits2$Migration_dist))

predData <- data.frame(newData, 
                       predict(bestMod, 
                               newdata = newData, 
                               interval = "confidence"))
# Plotting

(ForStrat.understory.plot <- ggplot(data = predData, 
       mapping = aes(x = ForStrat.understory, 
                     y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr), 
              alpha = 0.1) +
  geom_point(data = traits2, 
             mapping = aes(x = ForStrat.understory, 
                           y = freq_adj), 
             pch = 16, 
             size = 1.5) +
  labs(x = "% of foraging strategy in the understory", 
       y = "Vagrancy propensity (log)") )

# Handwing
newData <- expand.grid(Hand_wing = seq(min(traits2$Hand_wing), 
                                                 max (traits2$Hand_wing), 
                                                 length.out = 1000), 
                       ForStrat.understory = mean(traits2$ForStrat.understory),
                       Migration_dist = mean(traits2$Migration_dist))

predData <- data.frame(newData, 
                       predict(bestMod, 
                               newdata = newData, 
                               interval = "confidence"))

# Plot best model results

(Hand_wing.plot <- ggplot(data = predData, 
       mapping = aes(x = Hand_wing, 
                     y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr), 
              alpha = 0.1) +
  geom_point(data = traits2, 
             mapping = aes(x = Hand_wing, 
                           y = freq_adj), 
             pch = 16, 
             size = 2) +
  labs(x = "Hand-wing Index", 
       y = "Vagrancy propensity (log)") )

species_predictions <- arrangeGrob(ForStrat.understory.plot,Hand_wing.plot, ncol = 2)

ggsave("species_predictions.jpeg",
       plot = species_predictions,
       device = "jpeg",
       width = 2000,
       height = 1000,
       units = "px")

# Interaction plot
interact_plot(bestMod, 
              pred = Migration_dist, 
              modx = Hand_wing, 
              plot.points = TRUE,
              interval = TRUE,
              line.thickness = 2,
              point.size = 3) +
  theme_bw() +
  theme(legend.box.background = element_rect(colour = "black")) +
  labs(x="Migration Distance (km)",y="Vagrancy propensity (log)") +
  theme(text = element_text(size=20)) 

ggsave("trait_freq.jpeg",
       device = "jpeg",
       width = 2500,
       height = 2000,
       units = "px")

#for figure image
mean(traits2$Hand_wing) + sd(traits2$Hand_wing)
mean(traits2$Hand_wing) - sd(traits2$Hand_wing)

# Figure out where trendlines should end
max(subset(traits2, Hand_wing > mean(traits2$Hand_wing))$Migration_dist)
min(subset(traits2, Hand_wing > mean(traits2$Hand_wing))$Migration_dist)

max(subset(traits2, Hand_wing < mean(traits2$Hand_wing))$Migration_dist)
min(subset(traits2, Hand_wing < mean(traits2$Hand_wing))$Migration_dist)

# Check model assumptions ----

# Residuals

mean(bestMod$residuals)

par(mfrow = c(2,2))
plot(bestMod)

par(mfrow = c(1,1))
plot(bestMod$residuals ~ traits2$Hand_wing)
abline(h = 0)

plot(bestMod$residuals ~ traits2$Migration_dist)
abline(h = 0)

# Multicollinearity

car::vif(bestMod, type = "predictor")

# Influential points

cd <- cooks.distance(bestMod)
cd.idx <- cd >  3* mean(cd)
sum(cd.idx) # 3 points

plot(cd, 
     type = "h", 
     ylab = "Cook's distance")

abline(h = 3* mean(cd),  
       lty = 4, 
       col = "red") 

# rerun model with influential points removed
coef(bestMod)
coef(update(bestMod, 
               data = traits2[-cd.idx, ]))
# similar estimates

# Phylogenetic analysis ----

library(ape)
library(phytools)
library(picante)
library(tidyverse)

# Downloaded data from birdtree.org
#Load treefile
treefile <-read.nexus('output.nex')

#########################################################
# Code done with this tutorial: 
# https://brettkvo.com/how-to-create-a-ultrametric-dichotomous-phylogeny-in-r-using-vertlife-org-or-birdtree-org/


#Create a consensus tree
consensus<-consensus.edges(treefile, consensus.tree=consensus(treefile,p=0.5))  
plotTree(consensus, fsize=0.6)

#check for issues in edge length
print(consensus$edge.length) 

#correct edge length issues
n <- length(consensus$edge.length)  
#reduce the no. elements by 1 and add very small amount to make sure no branch length is zero   
consensus$edge.length <- consensus$edge.length[1:n-1] + .0000001  
print(consensus$edge.length) 

#Make ultrametric
consensus_ultra=chronos(consensus, lambda=0)  
#Check if ultrametric
is.ultrametric(consensus_ultra) 

#Make the tree dichotomous (so each node splits into only two branches)   
tree.unmatched <- multi2di(consensus_ultra, random=TRUE) 

#plot new tree
plotTree(tree.unmatched,fsize=0.6) 

####################################################
# Analysis from this tutorial:
# http://www.phytools.org/Cordoba2017/ex/5/Cont-char-models.html

#import offshore frequency data
offshore_freq2 <- traits %>% 
  subset(Common_name !="Varied Thrush")
offshore_freq2$Phylo <- sub(" ", "_", offshore_freq2$Phylo)

# Add names to tree
ordered_names <- left_join(as.data.frame(tree.unmatched$tip.label), offshore_freq2, 
                           by = c("tree.unmatched$tip.label" = "Phylo"))

tree.unmatched$tip.label <- ordered_names$Common_name

#set offshore frequency variable
var<-setNames(offshore_freq2$freq_adj,offshore_freq2$Common_name)
var

#test phylogenetic signal
phylosig(tree.unmatched,var,test=TRUE)

#plot trait on tree
obj<-contMap(tree.unmatched, var)


