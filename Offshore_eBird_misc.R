#
# Offshore eBird misc analysis and summary stats
# Mar 13, 2024
# Will Brooks
#
###########################

setwd("~/Desktop/eBird offshore")

require(tidyverse)
require(lubridate)
require(sf)
require(raster)
require(rasterVis)
require(RColorBrewer)
library(ggplot2)

#Load species info
eBird_offshore_species_info <- read.csv(file = "eBird_offshore_species_info.csv",sep=",")

#load extracted eBird data
eBird_offshore <- read.csv("eBird_covs2023-11-15.csv") %>% 
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

#load eBird data
offshore <- eBird_offshore %>% 
  dplyr::select("sampling_event_identifier","longitude","latitude") %>% 
  distinct(sampling_event_identifier, .keep_all = TRUE) 

#create raster
r <- raster(ncols=500, nrows=500)
xy <- offshore[,c(2,3)]
offshore_raster <- rasterize(xy, r, fun = "count")

plot(offshore_raster, 
     xlim = c(-130,-117),ylim = c(32,49))

#Set up point data for map
observations <- eBird_offshore %>% 
  subset(species_observed=="TRUE") %>%
  dplyr::select("locality_id","longitude","latitude","scientific_name") %>% 
  distinct(locality_id, .keep_all = TRUE) 

#Load land shapefile
box = c(xmin = -140, ymin = 30, xmax = -112, ymax = 50)
land_sf <- st_read("land/ne_10m_land.shp") %>% st_crop(box)

#plot with base R

#set color palette for raster
cuts=c(seq(1, 200, length.out = 10),500,1000,1450) #set breaks
col <- colorRampPalette(brewer.pal(9,"Oranges"))(12)

#plots
plot(offshore_raster, #base plot
     xlim = c(-130,-112),ylim = c(32.5,49),
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.sub = 1.5)
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     col = "lightskyblue2") # Color background (ocean) blue
plot(offshore_raster, breaks=cuts,col = col,
     xlim = c(-130,-112),ylim = c(32.5,49),add=TRUE) #create raster
plot(land_sf,col = "gray",
     xlim = c(-130,-112),ylim = c(32.5,49),add=TRUE)  #add land
points(observations$longitude,observations$latitude,
       pch = 24,
       cex=0.5,col="black", bg="white") #add observation points

#Summary stats
n_distinct(eBird_offshore$scientific_name) #number of species

n_tot <- eBird_offshore %>% 
  distinct(sampling_event_identifier, .keep_all = TRUE) 

nrow(n_tot) #number of checklists

sum(n_tot$effort_distance_km,na.rm=TRUE) #effort distance

sum(n_tot$duration_minutes,na.rm=TRUE)/60 #effort time in hours

n_obs <- eBird_offshore %>% subset(species_observed=="TRUE") 

max(n_obs$dist_offshore) #furthest distance offshore

mean(n_obs$dist_offshore) #mean distance offshore

sd(n_obs$dist_offshore) #standard deviation distance offshore

lists <- eBird_offshore %>% group_by(locality_id) %>%       #lists per location
  summarise(lists=n_distinct(sampling_event_identifier)) 

nrow(lists) #number of unique locations

lists100 <- subset(lists,lists>100) #locations with more than 100 lists

sum(lists100$lists) #how many were included hotspots

nrow(n_obs) #number of offshore passerine observations

n_obs$observation_count<-gsub("X","1",as.character(n_obs$observation_count))

sum(as.numeric(n_obs$observation_count)) #number of offshore passerine observations

nrow(offshore) #number of offshore checklists

nrow(observations)/nrow(n_tot)*100 #percent of checklists with at least one passerine

# plot histogram of observation distance
(dist_plot <- ggplot(n_obs,aes(x = dist_offshore)) +
  geom_histogram() +
    theme(text = element_text(size=30)) + 
  labs(x = "Distance offshore (m)",y = "Count",title = "A B"))

ggsave("dist_plot.jpeg",
       plot = dist_plot,
       device = "jpeg",
       width = 4000,
       height = 3000,
       units = "px")

##### Text analysis #######

# importing ebird data
files <- list.files(path = "ebd", full.names = T) #create list of all files

offshore_comments <- data.frame()

for (ebd in files){ # Loop
  comments <- read.delim(ebd) %>% 
    subset(SPECIES.COMMENTS != "") %>%
    mutate(checklist_id = SAMPLING.EVENT.IDENTIFIER,
           scientific_name = SCIENTIFIC.NAME) %>%
    dplyr::select("checklist_id","scientific_name","SPECIES.COMMENTS","AGE.SEX")
  
  join <- inner_join(comments,eBird_offshore,by=c("checklist_id","scientific_name"))
  
  offshore_comments <- rbind(offshore_comments,join)
}

# create variables representing presence of keywords
keywords <- offshore_comments %>%  
  mutate(juvenile = grepl("juvenile|juv|imm|immature|young|hatch year", ignore.case = TRUE, paste(SPECIES.COMMENTS,AGE.SEX)),
         adult = grepl("adult", ignore.case = TRUE, paste(SPECIES.COMMENTS,AGE.SEX)),
         female = grepl("female", ignore.case = TRUE, paste(SPECIES.COMMENTS,AGE.SEX)),
         male = grepl("male", ignore.case = TRUE, paste(SPECIES.COMMENTS,AGE.SEX)),
         north = grepl("north", ignore.case = TRUE, SPECIES.COMMENTS),
         south = grepl("south ", ignore.case = TRUE, SPECIES.COMMENTS),
         west = grepl("west", ignore.case = TRUE, SPECIES.COMMENTS),
         east = grepl("east", ignore.case = TRUE, SPECIES.COMMENTS))

keywords2 <- keywords %>% group_by(scientific_name) %>% summarise(
  juvenile = sum(juvenile),
  adult = sum(adult),
  female = sum(female),
  male = sum(male),
  north = sum(north),
  south = sum(south),
  west = sum(west),
  east = sum(east)
) 

write.csv(keywords2,"keywords.csv")

print(keywords2,n=30)
(sum(keywords2$juvenile)-sum(keywords2$adult))/sum(keywords2$adult)*100 #% change from adults to juveniles
