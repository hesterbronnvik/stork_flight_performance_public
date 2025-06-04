### Extract measures of how storks use thermals before they migrate
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 2024-11-12

library(moveACC)
library(sf)
library(gg3D)
library(tidyverse)
soar_gap <- 60
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
colfunc <- colorRampPalette(c("#fedee3", "#fcc2ce", "#feadbb", "#FE9AAB", "#FE869A", "#FD5D78"))

# the 2024 birds
ids <- c(4081089905, 4081050803, 4081110497, 4081168215, 4081121150, 4081043020, 4081147270, 3967948529,
         3968085304, 3968052282, 4000536122, 4000576074, 4000582268, 4000568914, 4060656577, 4060574425,
         4060647827, 4060597278, 4060632373, 3938980110, 3913673266, 4042698143, 4042709187, 4042735070,
         4042767872, 3939011728, 3939017397, 3967979972, 3968035956, 3968109777, 3968066416, 4000558649,
         4000545882, 3913664321, 3913696334, 3913704596, 3913714927, 4042813691)
# the nest sites, fledging dates, etc.
nests <- read.csv("/home/hbronnvik/Documents/explore_shackleton/nests_2024.csv") %>% 
  rename(nest_long = location.long, 
         nest_lat = location.lat) %>% 
  filter(fledging_date_manual_HB != "unclear") %>% 
  mutate(fledging_date_manual_HB = as.Date(fledging_date_manual_HB))

m_days <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2])) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day)

# from buzz_segmentation.R
seg_fls <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)
# the 2024 birds
seg_fls <- seg_fls[grepl(paste(ids, collapse = "|"), seg_fls)]

pre_migrations <- lapply(seg_fls, function(x){
  ind <- readRDS(x)
  rec <- m_days %>% 
    filter(individual.id == unique(ind$individual.id)) %>% 
    filter(date == min(date))
  ind <- ind %>% 
    left_join(nests) %>% 
    mutate(date = date(timestamp),
           onset = rec$date,
           rel_date = as.integer(difftime(date, onset, units = "days")),
           days_since_fledging = as.integer(difftime(date(timestamp), fledging_date_manual_HB, units = "days"))) %>% 
    filter(date >= fledging_date_manual_HB)
  tween_behav <- ind %>% 
    arrange(timestamp) %>% 
    group_by(track_flight_seg_id) %>% 
    # duration of the behavior
    mutate(len = n()) %>% 
    slice(1) %>% 
    ungroup() %>% 
    # whether the behavior is between tight circular soars
    mutate(tween = lag(flight_clust_sm3) == "circular_soaring" | lead(flight_clust_sm3) == "circular_soaring") %>% 
    # only shallow, brief, between behaviors
    filter(flight_clust_sm3 == "shallow_circular_soaring" & len <= 60 & tween == T) %>% 
    # the info needed to re-classify
    dplyr::select(track_flight_seg_id, flight_clust_sm3)
  ind <- ind %>% 
    # smooth a 4th time
    mutate(flight_clust_sm3 = ifelse(track_flight_seg_id %in% tween_behav$track_flight_seg_id,
                                     "circular_soaring", flight_clust_sm3))
  soar_ind <- ind %>% 
    filter(flight_clust_sm3 == "circular_soaring") %>%
    mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           time_diff = ifelse(is.na(time_diff), 1, time_diff),
           min_split = time_diff > 60,
           thermal_event = paste(individual.id, burst_id, cumsum(min_split)+1, sep = "_")) %>%
    group_by(thermal_event) %>% 
    mutate(thermal_duration = n(),
           vspeed_thermal = (height.above.ellipsoid[n()]-height.above.ellipsoid[1])/thermal_duration,
           turn_sd_thermal = sd(turn_angle, na.rm = T)/thermal_duration,
           turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
           directional_change = ifelse(turn_direction == lag(turn_direction), F, T))  %>%
    ungroup() %>% 
    dplyr::select(timestamp, thermal_event, thermal_duration, vspeed_thermal, turn_sd_thermal)
  ind <- ind %>% 
    left_join(soar_ind)
  gc()
  return(ind)
})
pre_migrations <- pre_migrations[!sapply(pre_migrations, function(x) nrow(x)==0)]
pre_migrations <- data.table::rbindlist(pre_migrations, fill = T)

# use the definitions needed for the wind etc. to define proper thermaling
real_therms <- pre_migrations %>% 
  filter(flight_clust_sm3 == "circular_soaring") %>% 
  group_by(thermal_event) %>% 
  # also add a 10 m height requirement because non-migratory storks do a lot of
  # curvy non-soaring
  mutate(climb = height.above.ellipsoid[n()]-height.above.ellipsoid[1]) %>% 
  filter(climb >= 10 & vspeed_thermal > 0 & thermal_duration >= 30) %>% 
  dplyr::select(thermal_event) %>% 
  deframe()

# count thermals on fledging dates
pre_migrations %>%  
  mutate(flight_clust_sm3 = ifelse(flight_clust_sm3 == "circular_soaring" & 
                                     !thermal_event %in% real_therms, NA, flight_clust_sm3)) %>% filter(flight_clust_sm3 == "circular_soaring") %>% 
  group_by(individual.id) %>% 
  slice(1) %>% 
  summarize(days = (unique(days_since_fledging))) %>% 
  group_by(days) %>% 
  summarize(obs = length(unique(individual.id)))

# generate Figure S3
c_sub <- pre_migrations %>%  
  mutate(flight_clust_sm3 = ifelse(flight_clust_sm3 == "circular_soaring" & 
                                     !thermal_event %in% real_therms, NA, flight_clust_sm3)) %>% filter(flight_clust_sm3 == "circular_soaring") %>% 
  group_by(individual.id) %>% 
  filter(days_since_fledging == 0) %>% 
  filter(thermal_event == unique(thermal_event)[1]) %>%
  mutate(obs = 1:n(),
         climb = height.above.ellipsoid[n()]-height.above.ellipsoid[1]) %>% 
  group_split()

# separate sub-plots
plts <- lapply(c_sub, function(p){
  color <- ifelse(p$climb[1] < 10, "gray50", "black")
  plt <- ggplot(p, aes(x = location.long, y = location.lat, 
                       z = height.above.ellipsoid, color = obs)) +
    axes_3D(phi=10) +
    stat_3D(phi=10, geom="path")  +
    stat_3D(phi=10) +
    scale_color_gradientn(colors = rev(colfunc(60))) +
    facet_wrap(~individual.id, scales = "free") +
    labs(color = "Second") +
    theme_void()
  return(plt)
})
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/day0.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(plts[[1]], plts[[2]], plts[[3]], plts[[4]], plts[[5]],
                  plts[[6]], plts[[7]], plts[[8]], plts[[9]])
# dev.off()


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/nestDist_2024.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
pre_migrations %>%  
  filter(date < onset) %>% 
  group_by(burst_id) %>% slice(1) %>% 
  group_by(individual.id, date) %>% 
  mutate(daily_dist = distVincentyEllipsoid(cbind(location.long[1], location.lat[1]),
                                            cbind(location.long[n()], location.lat[n()]))/1000,
         # daily_dist_sum = sum(distVincentyEllipsoid(cbind(location.long, location.lat),
         #                                            cbind(lag(location.long), lag(location.lat)))/1000,
         #                      na.rm = T),
         nest_dist = distVincentyEllipsoid(cbind(location.long, location.lat),
                                           cbind(nest_long, nest_lat))/1000) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(days_since_fledging) %>% 
  mutate(ids = length(unique(individual.id))) %>% 
  ggplot(aes(days_since_fledging, log(nest_dist))) +
  geom_smooth(color = "gray50") +
  geom_point(aes(color = ids)) +
  scale_x_continuous(n.breaks = 7) +
  scale_y_continuous(n.breaks = 7) +
  labs(x = "Days since fledging", y = "Distance to the nest (km)")
# dev.off()

### process ACC data for flapping flight
# the files containing raw acc data
acc <- list.files("/home/hbronnvik/Documents/WS_data/ACC_data", full.names = T)
# the files of the 2024 individuals with pre-migration GPS data
acc <- acc[grepl(paste(ids, collapse = "|"), acc)]

# the GPS data with ACC
gps_birds <- pre_migrations %>%  
  # only before migration
  filter(date <= onset) %>% 
  # one bird at a time
  group_by(individual.id) %>% 
  group_split()

# associate the ACC data to the GPS data using rounded off time stamps
start_time <- Sys.time()
flight_acc <- lapply(1:length(gps_birds), function(x){
  print(x)
  # load the GPS data
  gex <- gps_birds[[x]] %>% 
    rename(burstID = burst_id)%>% 
    mutate(date = date(timestamp), # add date
           minute = round_date(timestamp, "minute")) # round off the timestamp
           
  # if there are GPS data from migration, and there are ACC data:
  if(length(acc[grepl(unique(gex$individual.id), acc)]) > 0){
    # load the matching ACC data
    aex <- readRDS(acc[grepl(unique(gex$individual.id), acc)])
    # filter ACC data to the migrations
    aex <- aex %>% 
      filter(timestamp <= max(gex$timestamp)) %>% 
      # reduce by taking only ACC bursts following GPS bursts classified as flight
      filter(round_date(timestamp, "minute") %in% unique(round_date(gex$timestamp, "minute"))) %>% 
      # only consider 10 Hz samples
      filter(eobs_acceleration_sampling_frequency_per_axis == 10)
    
    if(nrow(aex) > 0){
      # attach whether this is a classified flight burst to the ACC data
      aex <- lapply(1:nrow(aex), function(n){
        # print(n)
        acc <- aex[n,]
        # use the time from the ACC
        ts <- round_date(acc$timestamp, "minute")
        # extract the GPS burst(s) containing that minute
        burst_oi <- unique(gex$burstID[gex$minute == ts])
        
        # take the GPS burst with the closest timestamp to the ACC
        burst_oi <- gex %>% 
          filter(burstID %in% burst_oi) %>% 
          dplyr::select(timestamp, burstID) %>% 
          mutate(prox = timestamp - acc$timestamp) %>% 
          filter(abs(prox) == min(abs(prox))) %>% 
          dplyr::select(burstID) %>% 
          deframe() %>% 
          unique()
        # the last location before the ACC burst started
        acc$location.long <- gex$location.long[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
        acc$location.lat <- gex$location.lat[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
        acc$gps_timestamp <- gex$timestamp[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
        acc$burstID <- burst_oi
        acc$date <- date(acc$timestamp)
        
        return(acc)
      }) %>% reduce(rbind)
      
      return(aex)
    }
  }
})
Sys.time()-start_time # Time difference of 25.39785 secs
# leaving 153 individuals
flight_acc <- flight_acc[!sapply(flight_acc, function(x) is.null(x))]

# estimate flapping following https://gitlab.com/anneks/moveACC
flapping <- lapply(1:length(flight_acc), function(a){
  print(paste0("Processing flight ACC data, individual ", a, " of ", length(flight_acc), "."), quote = F)
  acc <- flight_acc[[a]]
  if(nrow(acc) > 5){
    # ID <- unique(acc$individual_id)
    locs_prev_GPS <- acc %>% 
      dplyr::select(location.long, location.lat, gps_timestamp, individual_id, timestamp)
    waveDF <- ACCwave(acc, transformedData=F)
    # plot to determine best metrics:
    # clusterPlot(waveDF, cluster=F)
    # clusterPlot(waveDF, cluster=T, forclustering= c("varWaveWingBeat","eigenValue1"))
    # wingBeatsPlot(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"))
    # wingBeatsPlot(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), interactivePlot=F)
    # wingBeatsHist(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), interactivePlot=F)
    # once metrics are chosen:
    wingbeatsDF <- WingBeatsSelection(waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), 
                                      minbeat=2, maxbeat=5)
    wingbeatsDF <- wingbeatsDF %>% 
      # add on last location of the GPS burst preceeding the ACC
      left_join(locs_prev_GPS, by = join_by(timestamp))
    # saveRDS(wingbeatsDF, file = paste0("/home/hbronnvik/Documents/chapter2/flapping_flight/", ID, "_270323.rds"))
    return(wingbeatsDF)
  }
}) %>% reduce(rbind)

# check flapping rates over age
flapping <- flapping %>% 
  rename(individual.id = individual_id) %>% 
  left_join(nests %>% dplyr::select(individual.id, fledging_date_manual_HB)) %>% 
  # add on track IDs
  mutate(date = date(timestamp),
         days_since_fledging = round(as.numeric(difftime(timestamp, fledging_date_manual_HB, units = "days"))))

# confirm a connection between flapping estimate and DBA
flapping %>% 
  mutate(behavior = ifelse(behavior != "Flapping", "Passive", behavior)) %>% 
  ggplot(aes(as.factor(days_since_fledging), odbaMedian, fill = behavior)) +
  geom_boxplot() +
  labs(x = "Age", y = "Median ODBA", fill = "Flight") 
# count individuals in each day
flapping %>% 
  group_by(days_since_fledging) %>% 
  summarize(obs = length(unique(individualID)))
# count observations of each behavior
flapping %>% 
  group_by(behavior) %>% 
  summarize(obs = n())
# count observations on each day
flapping %>% 
  group_by(days_since_fledging) %>% 
  summarize(obs = n()) %>% 
  filter(days_since_fledging %in% c(10, 20, 30, 40, 50))

# visualize flapping over days since tagging
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/flapping_2024.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
# dev.off()


# generate Figure S4
pan1 <- pre_migrations %>%  
  filter(date < onset) %>% 
  mutate(flight_clust_sm3 = ifelse(flight_clust_sm3 == "circular_soaring" & 
                                     !thermal_event %in% real_therms, NA, flight_clust_sm3)) %>% 
  group_by(days_since_fledging) %>% 
  mutate(time_in_flight = sum(time_lag_sec)) %>% 
  ungroup() %>% 
  group_by(days_since_fledging, flight_clust_sm3) %>% 
  summarize(time_in_mode = sum(time_lag_sec)/unique(time_in_flight)) %>% 
  mutate(flight = ifelse(flight_clust_sm3 == "circular_soaring", "Soaring",
                         ifelse(flight_clust_sm3 == "gliding", "Gliding", "Other"))) %>% 
  drop_na(flight) %>% 
  ggplot(aes(days_since_fledging, time_in_mode, color = flight)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(n.breaks = 13, limits = c(0, 55),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#0081A7", "#C5CFB7", "#EE5E53")) +
  labs(x = "Days since fledging", y = "\nAmount of time\nin flight mode", color = "Flight type") 
pan2 <- pre_migrations %>%  
  group_by(rel_date) %>% 
  mutate(time_in_flight = sum(time_lag_sec)) %>% 
  ungroup() %>% 
  group_by(rel_date, flight_clust_sm3) %>% 
  summarize(time_in_mode = sum(time_lag_sec)/unique(time_in_flight)) %>% 
  mutate(flight = ifelse(flight_clust_sm3 == "circular_soaring", "Soaring",
                         ifelse(flight_clust_sm3 == "gliding", "Gliding", "Other"))) %>% 
  drop_na(flight) %>% 
  ggplot(aes(rel_date, time_in_mode, color = flight)) +
  geom_hline(yintercept = c(.2, .4, .6, .8), color = "gray50") +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(n.breaks = 13,
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#0081A7", "#C5CFB7", "#EE5E53")) +
  labs(x = "Days relative to migration", y = "\nAmount of time\nin flight mode", color = "Flight type") 
pan3 <- flapping %>% 
  group_by(days_since_fledging) %>% 
  summarize(total_obs = n(),
            flap_obs = sum(behavior=="Flapping"),
            flap_per = 100*(flap_obs/total_obs)) %>% 
  ggplot(aes(days_since_fledging, flap_per)) +
  geom_hline(yintercept = seq(0, 20, 2.5), color = "gray50") +
  geom_smooth(se = F, color = "grey50", lwd = 1.5) +
  geom_bar(position="stack", stat="identity", fill = "#FEADBB", width = 0.75) +
  scale_x_continuous(expand = c(0, 0), breaks = c(10, 20, 30, 40, 50), 
                     labels = c("10\n(n = 36)", "20\n(n = 112)", 
                                "30\n(n = 113)", "40\n(n = 26)", "50\n(n = 6)")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 31)) +
  labs(x = "Days since fledging", y = "\nPercentage of observations\nclassified as flapping")


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/flight_2024.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(pan1, pan2, pan3+theme(plot.margin = margin(7, 100, 7, 7)), 
                  nrow = 3, labels = "AUTO")
# dev.off()



# visualize where flapping happens
# bounding box
# ext <- st_sfc(st_point(c(5, 45)), st_point(c(10, 50)), crs = 4326)
# # DEM rasters
# library(terra)
# srtm <- list.files("/home/hbronnvik/Documents/Teaching/srtms", full.names = T)
# srtm <- sprc(srtm)
# srtm <- merge(srtm)
# # locations as a spatial object
# acc_locs <- flapping %>% 
#   dplyr::select(odbaMedian, behavior, location.long, location.lat) %>% 
#   st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)
# # map
# ggplot() +
#   tidyterra::geom_spatraster(data = srtm) +
#   tidyterra::scale_fill_grass_c(palette = "grey") +
#   geom_sf(data = acc_locs, aes(color = behavior)) +
#   coord_sf(xlim = st_coordinates(ext)[, "X"],
#            ylim = st_coordinates(ext)[, "Y"], expand = T) +
#   scale_color_manual("Flight", values = c("#EE5E53", "#0081A7")) +
#   labs(x = "Longitude", y = "Latitude", fill = "DEM") +
#   theme(legend.position = "right",
#         legend.text = element_text(size = 12))

# flapping <- flapping %>% 
#   mutate(dem = extract(srtm, vect(., geom=c("location.long", "location.lat"), crs = "+proj=longlat +datum=WGS84"), method = "bilinear")$ALB_elv_msk)
# 
# flapping %>% 
#   ggplot(aes(days_since_fledging, dem)) +
#   geom_point() +
#   geom_smooth(method = "lm", aes(group = as.factor(individual.id)))
