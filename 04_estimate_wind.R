### Using the 1 Hz soaring data, estimate wind vectors folloing Weinzierl et al., 2016 (https://doi.org/10.1002/ece3.2585)
### Hester Brønnvik
### 12.12.2023
### hbronnvik@ab.mpg.de

library(parallel)
library(move)
library(moveWindSpeed)
library(geosphere)
# library(data.table)
library(tidyverse)
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 12), 
                                text = element_text(size = 15)))

# functions for the estimation of wind speeds and other thermal features
# https://github.com/anflack/Updated-moveWindSpeed-
source("/home/hbronnvik/Documents/chapter2/getWindEstimates_update.R")
source("/home/hbronnvik/Documents/chapter2/thermallingFeaturesFunction.R")
source("/home/hbronnvik/Documents/chapter2/getTrackSegments_updated.R")
source("/home/hbronnvik/Documents/chapter2/getWindEstimate_update.R")

wgs <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")# map projection

# functions:
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}
wind_direction <- function(u, v){(90-(atan2(v, u)*(180/pi)))%%360 }
antiwind_direction <- function(u, v){(270-(atan2(v, u)*(180/pi)))%%360 }
# detach("package:plyr", unload = TRUE)

# identify migrations regardless of success
# the flight data from migrations segmented in 03_segment_flight.R
files <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)

# the migration dates generated in 02_segment_tracks.R
records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds")

# look at metrics of thermal use
thermal_results <- lapply(files, function(f){
  tr <- readRDS(f)
  # reclassify any shallow soaring as circular if it is brief and leads into or follows circular soaring
  tween_behav <- tr %>% 
    arrange(timestamp) %>% 
    group_by(track_flight_seg_id) %>% 
    # duration of the behavior
    mutate(len = n()) %>% 
    slice(1) %>% 
    ungroup() %>% 
    # whether the behavior is adjacent to tight circular soars
    mutate(tween = lag(flight_clust_sm3) == "circular_soaring" | lead(flight_clust_sm3) == "circular_soaring") %>% 
    # only shallow, brief, between behaviors
    filter(flight_clust_sm3 == "shallow_circular_soaring" & len <= 60 & tween == T) %>% 
    # the info needed to re-classify
    dplyr::select(track_flight_seg_id, flight_clust_sm3)
  tr <- tr %>% 
    # smooth a 4th time
    mutate(flight_clust_sm3 = ifelse(track_flight_seg_id %in% tween_behav$track_flight_seg_id,
                                     "circular_soaring", flight_clust_sm3)) %>% 
    filter(flight_clust_sm3 == "circular_soaring") %>% 
    rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
    mutate(ind_burst_id = paste0(individual.id, "_", burst_id))
  tr <- tr %>% 
    group_by(ind_burst_id) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    filter(obs > 29)
  tr
})
# 159 individuals left
thermal_results <- thermal_results[!sapply(thermal_results, function(x) is.null(x))]

# estimate wind vectors based on distortion of thermaling circles
wind_results <- lapply(1:length(thermal_results), function(n){
  print(n)
  tr <- thermal_results[[n]]
  if(nrow(tr) > 0){
    # split the data into separate thermaling events
    tr <- tr %>% 
      mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
             time_diff = ifelse(is.na(time_diff), 1, time_diff),
             min_split = time_diff > 60,
             thermal_event = paste(ind_burst_id, cumsum(min_split)+1, sep = "_")) %>% 
      dplyr::select(-time_diff, -min_split)
    # use long enough thermals to estimate wind vectors
    ind <- lapply(split(tr, tr$thermal_event), function(burst){
      burst <- burst %>% as.data.frame()
      # transform to a move object because the windEstimates require one
      mv_burst <- move(x=burst$location.long, y=burst$location.lat, 
                       time=burst$timestamp,
                       proj="+proj=longlat +datum=WGS84 +no_defs",
                       animal=burst$individual.id,
                       data=burst)
      # identify bursts with 30 or more seconds of soaring
      burst_summary <- burst %>%
        mutate(gap = c(1,timeLag(mv_burst)),
               consistent = gap == 1,
               event = cumsum(consistent==0)) %>%
        group_by(event) %>%
        summarize(obs = n()) %>% 
        mutate(sufficient = obs>29)
      # get the actual estimates
      if(T %in% unique(burst_summary$sufficient)){
        class_burst <- getWindEstimates(mv_burst)
        class_burst <- as.data.frame(class_burst)
        return(class_burst)
      }
    })
    ind <- ind[!sapply(ind, function(x) is.null(x))]
    ind <- data.table::rbindlist(ind, use.names = T)
    # save out
    saveRDS(ind, file = paste0("/home/hbronnvik/Documents/WS_data/wind_thermal_250deg_20s/", unique(ind$individual.id), "_seg_wind_",Sys.Date(),".rds"))
    return(ind)
  }
})
