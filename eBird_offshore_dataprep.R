#
# Offshore eBird data prep
# 8 Nov 2023
# Will Brooks
#
###########################
setwd("~/Desktop/eBird offshore")

#load packages
library(sf)
library(tidyverse)
library(auk)
library(ggplot2)
library(car)

######## Designate points as land or ocean and add distance from shore data

# Load land shapefile
box = c(xmin = -140, ymin = 30, xmax = -112, ymax = 50)
land_sf <- st_read("land/ne_10m_land.shp") %>% st_crop(box)
land_trans <- st_transform(land_sf, 32610) #transform to UTM zone10

points_sf <- read.delim("complete.sed.txt") %>% 
  select("LOCALITY.ID","LONGITUDE","LATITUDE") %>% 
  distinct(LOCALITY.ID, .keep_all = TRUE) %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = "EPSG:4326") %>% #Convert ebird data to simple feature
  st_transform(32610) #Transform to UTM zone 10

points_cats <- st_join(points_sf,land_trans)
ocean_points <- subset(points_cats, is.na(featurecla))
land_points <- subset(points_cats, featurecla == "Land")

# distance offshore 
distances <- st_geometry(obj = land_trans) %>% 
  st_distance(y = st_as_sf(ocean_points))
dist_offshore <- apply(distances, MARGIN = 2, FUN = min) %>% 
  cbind(as.data.frame(ocean_points)[,1]) %>% 
  as.data.frame()
colnames(dist_offshore) <- c("dist_offshore","locality_id")
dist_offshore$dist_offshore <- as.numeric(dist_offshore$dist_offshore)

# distance inland
distances <- st_geometry(obj = land_trans) %>% 
  st_cast(to = 'MULTILINESTRING') %>% 
  st_distance(y = st_as_sf(land_points))
dist_inland <- apply(distances, MARGIN = 2, FUN = min) %>% 
  cbind(as.data.frame(land_points)[,1]) %>% 
  as.data.frame() 
colnames(dist_inland) <- c("dist_inland","locality_id")
dist_inland$dist_inland <- as.numeric(dist_inland$dist_inland)
dist_inland <- subset(dist_inland, dist_inland < 10000) # select checklists within 10 km of ocean

#clear out some space
rm(land_sf,land_trans,points_sf,points_cats,ocean_points,land_points,distances)

###### importing ebird data
files <- list.files(path = "ebd", full.names = T) #create list of all files
sed <- c("complete.sed.txt")

# set vectors for cutting, I matched this to the scale of the weather data
breaks <- breaks <- c(32.5,seq(33.75, 48, by = 2.5),48)
labels <- paste0("cut", seq(32.5, 48, by = 2.5))

# create empty df for output
eBird_offshore <- data.frame()
# loop through species 
for (ebd in files){

  # Zerofill  
  zf <- auk_zerofill(ebd, sed, collapse = TRUE) %>% 
    mutate(lat_cut=cut(latitude, breaks=breaks, labels=labels))
  
  # Select ocean and land points
  ocean <- right_join(zf,dist_offshore, by="locality_id") %>% # select offshore checklists
    subset(dist_offshore > 10000) %>% # select 10 km offshore and beyond
    drop_na(scientific_name) # remove some rows where a location with no data was selected
  land <- right_join(zf,dist_inland, by="locality_id") %>% #select on land checklists
    drop_na(scientific_name) # remove some rows where a location with no data was selected
  
  # mark observations
  ocean_obs <- subset(ocean,species_observed=="TRUE")
    
  # add on land frequencies to ocean observations

  ocean <- land %>% 
    group_by(lat_cut,observation_date) %>%  #group land by latitude group and date
    reframe(land_obs_freq = mean(as.logical(species_observed))) %>% #calculate obs freq on land for each day
    unique() %>%
    right_join(ocean, #add on land obs freq to offshore observations
                              by=c("lat_cut","observation_date")) 
  
  eBird_offshore <- rbind(eBird_offshore,ocean)
}

# remove accidental duplicate data
eBird_offshore2 <- eBird_offshore %>% 
  unique() %>% 
  subset(scientific_name != "Loxia curvirostra") #leaving out red crossbill because it looks to be irruptive rather than migratory

write.csv(eBird_offshore,file=paste0("eBird_offshore",Sys.Date(),".csv"))

########## Do covariates now

setwd("~/Desktop/eBird offshore")

library(RNCEP)
library(tidyverse)
library(abind)

eBird_offshore <- read.csv("eBird_offshore2023-11-15.csv")[,-1]
eBird_offshore$observation_date <- ymd(eBird_offshore$observation_date)

#select years to get weather data
years <- unique(year(eBird_offshore$observation_date)) %>% sort()
years <- years[years > 1947] # select data after first year of weather data
#grouping_variable <- cumsum(c(TRUE, diff(years) != 1)) #this code is just to split non-sequential years
#result_list <- split(years, grouping_variable)       #all of our years are sequential so it isn't used here

#define the necessary arguments
month_range_1 <- c(1,12)     #period of months
year_range_1 <- c(min(years),max(years)-1)
month_range_2 <- c(1,9)     #period of months
year_range_2 <- c(max(years),max(years))

lat_range <- c(32, 48)      #latitude range
lon_range <- c(-130, -118)     #longitude range

# get weather data
uwind_1 <- NCEP.gather("uwnd.sig995",    #name of the variable
                     "surface", 
                    month_range_1,year_range_1,
                    lat_range,lon_range,
                    return.units = TRUE,
                    reanalysis2=FALSE)

uwind_2 <- NCEP.gather("uwnd.sig995",    #name of the variable
                     "surface", 
                     month_range_2,year_range_2,
                     lat_range,lon_range,
                     return.units = TRUE,
                     reanalysis2=FALSE)

uwind <- abind(uwind_1,uwind_2, along = 3)

vwind_1 <- NCEP.gather("vwnd.sig995",    #name of the variable
                     "surface", 
                     month_range_1,year_range_1,
                     lat_range,lon_range,
                     return.units = TRUE,
                     reanalysis2=FALSE)

vwind_2 <- NCEP.gather("vwnd.sig995",    #name of the variable
                     "surface", 
                     month_range_2,year_range_2,
                     lat_range,lon_range,
                     return.units = TRUE,
                     reanalysis2=FALSE)

vwind <- abind(vwind_1,vwind_2, along = 3)

rel_hum_1 <- NCEP.gather("rhum.sig995",    #name of the variable
               "surface", 
               month_range_1,year_range_1,
               lat_range,lon_range,
               return.units = TRUE,
               reanalysis2=FALSE)

rel_hum_2 <- NCEP.gather("rhum.sig995",    #name of the variable
                       "surface", 
                       month_range_2,year_range_2,
                       lat_range,lon_range,
                       return.units = TRUE,
                       reanalysis2=FALSE)

rel_hum <- abind(rel_hum_1,rel_hum_2, along = 3)

# get daily averages

# uwind
date_time <- dimnames(uwind)[[3]] # make date a grouping variable
date_time <- ymd_h(date_time)
group <- date(date_time)

uwind_daily_avg <- aperm( # get averages per day
  apply(
    uwind, #our data
    c(1,2), #apply to each time series 1:row, 2:column a the mean( ) function
    by, #group by
    group, #months
    function(x)ifelse(all(is.na(x)),NA,mean(x))),
  c(2,3,1)) #reorder to get an array like the original

# Create a data frame with all combinations of the three variables
combinations <- expand.grid(
  lat=dimnames(uwind_daily_avg)[[1]],
  lon=dimnames(uwind_daily_avg)[[2]],
  date=dimnames(uwind_daily_avg)[[3]]
)

# Convert the three-dimensional array to a data frame
uwind_df <- as.data.frame(matrix(uwind_daily_avg, ncol = 1)) %>% 
  cbind(combinations) %>% 
  rename(uwind=V1) %>% 
  mutate(lon=as.numeric(as.character(lon))-360)

uwind_df <- as.data.frame(lapply(uwind_df, as.character))

# vwind
date_time <- dimnames(vwind)[[3]] # make date a grouping variable
date_time <- ymd_h(date_time)
group <- date(date_time)

vwind_daily_avg <- aperm( # get averages per day
  apply(
    vwind, #our data
    c(1,2), #apply to each time series 1:row, 2:column a the mean( ) function
    by, #group by
    group, #months
    function(x)ifelse(all(is.na(x)),NA,mean(x))),
  c(2,3,1)) #reorder to get an array like the original

# Create a data frame with all combinations of the three variables
combinations <- expand.grid(
  lat=dimnames(vwind_daily_avg)[[1]],
  lon=dimnames(vwind_daily_avg)[[2]],
  date=dimnames(vwind_daily_avg)[[3]]
)

# Convert the three-dimensional array to a data frame
vwind_df <- as.data.frame(matrix(vwind_daily_avg, ncol = 1)) %>% 
  cbind(combinations) %>% 
  rename(vwind=V1) %>% 
  mutate(lon=as.numeric(as.character(lon))-360)

vwind_df <- as.data.frame(lapply(vwind_df, as.character))

# humidity
date_time <- dimnames(rel_hum)[[3]] # make date a grouping variable
date_time <- ymd_h(date_time)
group <- date(date_time)

relh_daily_avg <- aperm( # get averages per day
  apply(
    rel_hum, #our data
    c(1,2), #apply to each time series 1:row, 2:column a the mean( ) function
    by, #group by
    group, #months
    function(x)ifelse(all(is.na(x)),NA,mean(x))),
  c(2,3,1)) #reorder to get an array like the original

# Create a data frame with all combinations of the three variables
combinations <- expand.grid(
  lat=dimnames(relh_daily_avg)[[1]],
  lon=dimnames(relh_daily_avg)[[2]],
  date=dimnames(relh_daily_avg)[[3]]
)

# Convert the three-dimensional array to a data frame
relh_df <- as.data.frame(matrix(relh_daily_avg, ncol = 1)) %>% 
  cbind(combinations) %>% 
  rename(relh=V1) %>% 
  mutate(lon=as.numeric(as.character(lon))-360)

relh_df <- as.data.frame(lapply(relh_df, as.character))

# merging sampling data with weather
eBird_sampling <- eBird_offshore %>% select("observation_date","latitude","longitude") %>% 
  unique() %>% 
  mutate(latitude_cat = case_when(latitude < 33.75 ~ 32.5,
                                  latitude > 33.75 & latitude < 36.25 ~ 35,
                                  latitude > 36.25 & latitude < 38.75 ~ 37.5,
                                  latitude > 38.75 & latitude < 41.25 ~ 40,
                                  latitude > 41.25 & latitude < 43.75 ~ 42.5,
                                  latitude > 43.75 & latitude < 46.25 ~ 45,
                                  latitude > 46.25 ~ 47.5
  ),
       longitude_cat = case_when(longitude < -128.75 ~ -130,
                                 longitude > -128.75 & longitude < -126.25 ~ -127.5,
                                 longitude > -126.25 & longitude < -123.75 ~ -125,
                                 longitude > -123.75 & longitude < -121.25 ~ -122.5,
                                 longitude > -121.25 ~ -120,
                                 longitude > -118.75 ~ -117.5
  ))

eBird_sampling <- as.data.frame(lapply(eBird_sampling, as.character))

eBird_weather <- left_join(eBird_sampling,uwind_df,
                           by=c("observation_date"="date",
                                "latitude_cat"="lat",
                                "longitude_cat"="lon")) %>% 
                  left_join(vwind_df,
                           by=c("observation_date"="date",
                                "latitude_cat"="lat",
                                "longitude_cat"="lon")) %>%
                  left_join(relh_df,
                           by=c("observation_date"="date",
                                "latitude_cat"="lat",
                                "longitude_cat"="lon"))

eBird_offshore <- as.data.frame(lapply(eBird_offshore, as.character))

######## Geomagnetic disturbance and solar radiation #############

#load geomagnetic disturbance data
geodist <- read.table("Kp_1932-01-01_2023-10-31_D.dat",
                      skip=37)[,c(1,5)] 
colnames(geodist) <- c("date","ap")
geodist$date <- as.Date(geodist$date, format= "%Y-%m-%d")

geodist2 <- geodist %>% group_by(date) %>% summarize(ap_mean = mean(ap))

#load solar radiation data
solar <- read.csv("american_relative_sunspot_number_daily.csv")
solar$time..yyyy.MM.dd. <- as.Date(solar$time..yyyy.MM.dd., format= "%Y %m %d")

#load eBird checklists
lists <- eBird_offshore %>% select("checklist_id","observation_date") %>% 
  distinct(checklist_id, .keep_all = TRUE)
lists$observation_date <- as.Date(lists$observation_date, format= "%Y-%m-%d")

geosolar <- lists %>% left_join(geodist2, by=c("observation_date"="date")) %>% 
  left_join(solar, by=c("observation_date"="time..yyyy.MM.dd.")) %>% 
  select("checklist_id","ap_mean","Ra")

############### combine everything

eBird_covs <- full_join(eBird_offshore,eBird_weather,by=c("observation_date",
                                                          "latitude",
                                                          "longitude")) %>% 
  full_join(geosolar,by="checklist_id")

write.csv(eBird_covs,file=paste0("eBird_covs",Sys.Date(),".csv"))

############### lets just look at the data
eBird_covs <- read.csv("eBird_covs2023-11-15.csv")

ggplot(eBird_covs,aes(y=species_observed,x=uwind)) + geom_boxplot()
ggplot(eBird_covs,aes(y=species_observed,x=vwind)) + geom_boxplot()
ggplot(eBird_covs,aes(y=species_observed,x=relh)) + geom_boxplot()
ggplot(eBird_covs,aes(y=species_observed,x=ap_mean)) + geom_boxplot()
ggplot(eBird_covs,aes(y=species_observed,x=Ra)) + geom_boxplot()

