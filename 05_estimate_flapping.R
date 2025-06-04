### Estimate flapping rates for white storks migrating
### Hester Brønnvik
### 2024-03-21
### hbronnvik@ab.mpg.de

library(move)
library(moveACC)
library(sf)
library(tidyverse)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  panel.background = element_rect(fill = "white"),
                  strip.background = element_rect(fill = "white", color = "white"),
                  strip.text = element_text(face = "bold")))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# colfunc <- colorRampPalette(c("#8A2846", "#B9375E", "#E05780", "#FF7AA2", "#FF9EBB", "#FFC2D4", "#FFE0E9"))
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633, 908232414)

# Here, our first step is to use the GPS flight data to identify which ACC bursts were
# recorded in flight. Then, we can use those ACC data to identify flapping flight.

# get metadata for birds with ACC data
info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 2365683 & grepl("acceleration", sensor_type_ids)) %>% 
    mutate(study = x)
  return(info)
}) %>% reduce(rbind)

# the migration dates generated in 02_segment_tracks.R
records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2]),
         season = ifelse(grepl("fall", trackID), "Fall", "Spring")) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day, season)

# the migration dates generated in 02_segment_tracks.R
meta <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2]),
         season = ifelse(grepl("fall", trackID), "Fall", "Spring")) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day, season) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(migration = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()

# the flight data from migrations segmented in 03_segment_flight.R
gps_files <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)#[grepl("173673348|173661590|23460438", list.files("/home/hbronnvik/Documents/chapter2/classified_data", full.names = T))]
# all the ACC data from 01_access_data.R
acc_files <- list.files("/home/hbronnvik/Documents/WS_data/ACC_data", full.names = T)

# the birds with ACC data
acc_birds <- lapply(acc_files, function(f){
  sub("/home/hbronnvik/Documents/WS_data/ACC_data/", "", 
      sub("_20241026_full_ACC_data.rds", "", f))
}) %>% unlist()
# the GPS data with ACC
gps_birds <- data.frame(file = gps_files) %>% 
  rowwise() %>% 
  mutate(bird = sub("_segmented_bursts.rds", "", 
                    sub("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s/", "", 
                        file))) %>% 
  filter(bird %in% acc_birds)

# associate the ACC data to the GPS data using rounded off time stamps
start_time <- Sys.time()
flight_acc <- lapply(1:length(gps_files), function(x){
  print(x)
  # load the GPS data
  gex <- readRDS(gps_files[x]) %>% 
    rename(burstID = burst_id)
  # get the dates of the migration
  rec <- records %>% 
    filter(individual.id == unique(gex$individual.id))
  
  if(nrow(rec) > 0){
    # filter the GPS data to migration and add a rounded-off minute
    gex <- lapply(split(rec, rec$trackID), function(track){
      ex <- gex %>% 
        mutate(date = date(timestamp), # add date
               minute = round_date(timestamp, "minute"), # round off the timestamp
               trackID = unique(track$trackID)) %>% # add the track ID
        # take the GPS data on migratory days for the given track
        filter(date %in% track$date)
      return(ex)
    }) %>% reduce(rbind)
    # if there are GPS data from migration, and there are ACC data:
    if(nrow(gex) > 0){
      if(length(acc_files[grepl(unique(gex$individual.id), acc_files)]) > 0){
        # load the matching ACC data
        aex <- readRDS(acc_files[grepl(unique(gex$individual.id), acc_files)])
        # s_id <- info$study[info$individual.id == unique(gex$individual.id)]
        # aex <- getMovebankNonLocationData(study = s_id, animalName = unique(gex$individual.id),
        #                                   sensorID = 2365683, login = loginStored)
        # filter ACC data to the migrations
        aex <- lapply(split(rec, rec$trackID), function(track){
          ex <- aex %>% 
            filter(date(timestamp) %in% track$date) %>% 
            mutate(trackID = unique(track$trackID))
          return(ex)
        }) %>% reduce(rbind)
        
        aex <- aex %>% 
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
            # if there is a burst, then the ACC followed GPS classified flight, else it did not
            if(length(burst_oi) == 0){
              acc$flight <- F
              acc$location.long <- NA
              acc$location.lat <- NA
              acc$burstID <- NA
            }else{
              # take the GPS burst with the closest timestamp to the ACC
              burst_oi <- gex %>% 
                filter(burstID %in% burst_oi) %>% 
                dplyr::select(timestamp, burstID) %>% 
                mutate(prox = timestamp - acc$timestamp) %>% 
                filter(abs(prox) == min(abs(prox))) %>% 
                dplyr::select(burstID) %>% 
                deframe() %>% 
                unique()
              acc$flight <- T
              # the last location before the ACC burst started
              acc$location.long <- gex$location.long[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
              acc$location.lat <- gex$location.lat[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
              acc$burstID <- burst_oi
            }
            return(acc)
          }) %>% reduce(rbind)
          
          # keep the ACC data that were on migration and in flight
          acc_remaining <- aex %>% 
            filter(flight == T) %>% 
            mutate(date = date(timestamp)) %>% 
            left_join(rec[, c("date", "ld_day")], by = join_by(date)) %>% 
            left_join(meta[, c("trackID", "migration")], by = join_by(trackID)) 
          return(acc_remaining)
        }
      }
    }
  }
})
Sys.time()-start_time # Time difference of 1.181882 hours
# leaving 153 individuals
flight_acc <- flight_acc[!sapply(flight_acc, function(x) is.null(x))]

# saveRDS(flight_acc, file = "/home/hbronnvik/Documents/chapter2/flight_acc_10Hz_20241210.rds")

flight_acc <- readRDS(file = "/home/hbronnvik/Documents/chapter2/flight_acc_10Hz_20241210.rds")

# check <- lapply(flight_acc, function(ch){
#   ch <- ch %>%
#     mutate(yr = year(timestamp)) %>%
#     dplyr::select(eobs_acceleration_sampling_frequency_per_axis, yr, individual_id) %>%
#     group_by(eobs_acceleration_sampling_frequency_per_axis) %>%
#     summarize(yr = unique(yr),
#               individual.id = unique(individual_id))
#   return(ch)
# }) %>% reduce(rbind)

# estimate flapping following https://gitlab.com/anneks/moveACC
flapping <- lapply(1:length(flight_acc), function(a){
  print(paste0("Processing flight ACC data, individual ", a, " of ", length(flight_acc), "."), quote = F)
  acc <- flight_acc[[a]]
  if(nrow(acc) > 5){
    # ID <- unique(acc$individual_id)
    locs_prev_GPS <- acc %>% 
      dplyr::select(location.long, location.lat, timestamp, individual_id, trackID)
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
  # add on track IDs
  mutate(date = date(timestamp)) %>% 
  rename(individual.id = individual_id) %>% 
  left_join(records) %>% 
  # add on number of migrations
  left_join(meta %>% dplyr::select(trackID, migration))

# saveRDS(flapping, file = "/home/hbronnvik/Documents/chapter2/flapping_flight_10Hz_241210.rds")

flapping <- readRDS("/home/hbronnvik/Documents/chapter2/flapping_flight_10Hz_241210.rds") %>% 
  # drop eastern IDs
  filter(!individual.id %in% c(1173989698, 1174005664)) %>% 
  # this individual migrated 732 km and survived, but exclusively in Sub-Sahara
  filter(trackID != "1176038499_spring_2021") 

# df <- getMovebankLocationData(24442409, sensorID = 653,
#                               animalName = 1173982416, login = loginStored)
# df %>% 
#   drop_na(location.long) %>% 
#   mutate(date = date(timestamp)) %>% 
#   filter(between(date, as.Date("2020-08-17"), as.Date("2020-08-30"))) %>% 
#   mutate(date = as.character(date)) %>% 
#   st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
#   mapview::mapview(zcol = "date")

flapping %>% 
  group_by(migration, season) %>% 
  summarize(obs = length(unique(individualID)))

flapping %>% 
  group_by(behavior) %>% 
  summarize(obs = n())

flapping %>% 
  filter(migration < 5) %>% 
  mutate(behavior = ifelse(behavior != "Flapping", "Passive", behavior)) %>% 
  ggplot(aes(as.factor(migration), odbaMedian, fill = behavior)) +
  geom_boxplot() +
  labs(x = "Age", y = "Median ODBA", fill = "Flight") +
  facet_wrap(~season)

library(terra)
file <- "/home/hbronnvik/Documents/storkSSFs/ecmwf/single/boundary_layer_height_land_2022.nc"
blh <- rast(file)
template <- blh[[1]]
terra::plot(template)

# Hobo-Dyer
target_crs <- "+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

library("rnaturalearth")
world <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")
# bounding box
ext <- st_sfc(st_point(c(-20.125, 10.125)), st_point(c(20.125, 58)), crs = 4326) %>% 
  st_transform(target_crs)

base <- world %>% 
  st_transform(target_crs) 

mps <- lapply(c("fall", "spring"), function(s){
  
  if(s == "spring"){
    colfunc <- colorRampPalette(c("#9DD2C9",  "#7fc4b8", "#00AFB9", "#0098b0", "#0081A7", "#006D8F","#004E66"))#, "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
  }else{
    colfunc <- colorRampPalette(c("#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53", "#EB3F33"))
  }
  
  if(s == "spring"){
    lab <- "Spring"
  }else{
    lab <- "Fall"
  }
  
  locs <- flapping %>% 
    filter(season == s) %>% 
    dplyr::select(season, location.long, location.lat, journey_number, behavior) %>% 
    mutate(flapping = behavior=="Flapping") %>% 
    vect(geom=c("location.long", "location.lat"), crs = "+proj=longlat +datum=WGS84") %>% 
    # rasterize using ECMWF .25 degree template by summing the logical of flapping/passive
    terra::rasterize(template, field = "flapping", fun = sum)  %>% 
    project(target_crs)
  # divide by the number of cells
  locs <- locs/as.numeric(global(locs, fun="notNA"))
  
  mp <- ggplot() +
    geom_sf(data = base, fill = "grey60", color = "grey30") +
    tidyterra::geom_spatraster(data = locs) +
    coord_sf(xlim = st_coordinates(ext)[, "X"],
             ylim = st_coordinates(ext)[, "Y"], expand = T) +
    scale_fill_gradientn("Flapping", colors = colfunc(30), na.value = "transparent") +
    # tidyterra::scale_fill_whitebox_c(palette = "muted") +
    # guides(fill = "none") +
    labs(x = "Longitude", y = "Latitude", fill = "Flapping\nproportion") +
    # theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 12))
  return(mp)
})
ggpubr::ggarrange(mps[[1]], mps[[2]])

pl <- lapply(c("spring", "fall"), function(s){
  col <- ifelse(s == "spring", "#007BA7", "#F07268")
  # lab <- ifelse(s == "spring", "Spring", "Fall")
  # plot bars
  prop_losses <- flapping %>%
    filter(migration < 5 & season == s) %>%
    mutate(flapping = behavior == "Flapping") %>% 
    group_by(season, migration) %>% 
    summarize(flaps = sum(flapping),
              obs = n(),
              prop = 100*(flaps/obs),
              label = paste0(unique(migration), "\n(n=", unique(obs), ")")) %>% 
    ungroup() %>% 
    ggplot(aes(x = label, y = prop)) +
    geom_hline(yintercept = seq(0, 20, 2.5), color = "gray50") +
    geom_bar(position="stack", stat="identity", width = 0.5, fill = col) +
    # scale_x_discrete(expand = c(0, 0 )) +
    scale_y_continuous(expand = c(0, 0 ), limits = c(0, 23), n.breaks = 10, labels = scales::label_number(accuracy = 1)) +
    geom_text(aes(label=paste0(sprintf("%1.1f", prop),"%")), y = 22,
              # position=position_stack(vjust=0.5), size = 3.5,
              color = "black", fontface = "bold") +
    # scale_fill_manual(values = c("#EE5E53", "#0081A7", "#F07268", "#009BB1", "#F38979", "#24B5B8", "#FABBA0", "#B5CDB7")) +
    labs(y = "Percentage of observations\nclassified as flapping",
         x = "Age (years)", fill = "Flapping") +
    theme(legend.position = "none")
  return(prop_losses)
})


maps <- ggpubr::ggarrange(mps[[1]], mps[[2]], nrow = 1, labels = "AUTO")
bars <- ggpubr::ggarrange(pl[[2]], pl[[1]], nrow = 1, labels = c("C", "D"))

# create Figure 4
tgrob1 <- ggpubr::text_grob("Fall", size = 20)
tgrob2 <- ggpubr::text_grob("Spring", size = 20)
# Draw the text
plot_0 <- ggpubr::as_ggplot(tgrob1) + theme(plot.margin = margin(0,12,0,0, "cm"))
plot_1 <- ggpubr::as_ggplot(tgrob2) + theme(plot.margin = margin(0,12,0,0, "cm"))

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/flapping_bars_10Hz.png",
#     height = 8.5/2, width = 11, units = "in", res = 500)
ggpubr::ggarrange(ggpubr::ggarrange(plot_0, plot_1, nrow = 1),
                  ggpubr::ggarrange(pl[[2]], pl[[1]], nrow = 1, labels = "AUTO"),
                  nrow = 2, heights = c(1, 10))
# dev.off()

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/flapping_locs2.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(ggpubr::ggarrange(plot_0, plot_1, nrow = 1), 
                  maps, 
                  bars, 
                  nrow = 3, heights = c(1, 5, 5))

# dev.off()

flapping %>%
  filter(migration < 5) %>% 
  group_by(season, migration, behavior) %>% 
  summarize(obs = n()) %>% 
  pivot_wider(names_from = "behavior", values_from = obs) %>% 
  mutate(total = sum(Flapping, NotFlapping),
         prop_flap = Flapping/total)

selected_data <- flapping %>%
  filter(migration < 5 & season == "fall") %>% 
  dplyr::select(season, migration, behavior) 

# Create a contingency table
contingency_table <- table(selected_data$migration, selected_data$behavior)
# View the contingency table
print(contingency_table)
# Perform chi-square test
chi_square_test <- chisq.test(contingency_table)
# View the results
print(chi_square_test)
# Observed counts
observed_counts <- chi_square_test$observed
print(observed_counts)
# Expected counts
expected_counts <- chi_square_test$expected
print(round(expected_counts, 2))
# Pearson residuals
pearson_residuals <- chi_square_test$residuals
print(round(pearson_residuals, 2))
# Calculate contribution to chi-square statistic
contributions <- (observed_counts - expected_counts)^2 / expected_counts
# Calculate percentage contributions
total_chi_square <- chi_square_test$statistic
percentage_contributions <- 100 * contributions / total_chi_square
# Print percentage contributions
print(round(percentage_contributions, 2))

# next, look at flapping in relation to ground speeds
# the cleaned GPS data
clean_locations <- readRDS("/home/hbronnvik/Documents/chapter2/clean_data/clean_locations_2024-11-04.rds")

ground_days <- lapply(clean_locations, function(ind){
  # ind <- clean_locations[[2]]
  if(max(ind$location.long) <= 16){
    # its migration days
    rec <- records %>% 
      filter(individual.id == unique(ind$individual.id))
    ind <- ind %>% 
      as.data.frame() %>% 
      mutate(date = date(timestamp)) %>% 
      # only the migration days
      filter(date %in% rec$date) %>% 
      group_by(date) %>%
      arrange(timestamp) %>% 
      slice(1, n()) %>% 
      mutate(daily_time = as.numeric(difftime(timestamp[2], timestamp[1], units = "days"))) %>% 
      # the distance covered
      dplyr::select(individual.id, date, daily_dist, daily_time) %>% 
      slice(1) %>%
      mutate(daily_speed = (daily_dist/1000)/daily_time) %>%  # ground speed km/day
      ungroup()
    return(ind)
  }
}) %>% reduce(rbind)

rm(clean_locations);gc()

# add daily ground speeds to flapping
# there are some missing speeds?!
daily_flapping <- flapping %>% 
  mutate(flapping = behavior == "Flapping") %>% 
  group_by(date, individual.id) %>% 
  mutate(flaps = sum(flapping),
         obs = n(),
         daily_flaps = flaps/obs) %>% 
  slice(1) %>% 
  ungroup()%>% 
  left_join(ground_days)

summary(daily_flapping$daily_speed)

# flapping <- flapping %>% 
#   mutate(flapping = behavior == "Flapping") %>% 
#   group_by(date) %>% 
#   mutate(flaps = sum(flapping),
#          obs = n(),
#          daily_flaps = flaps/obs) %>% 
#   ungroup()

ggplot(daily_flapping, aes(log(daily_flaps), log(daily_speed))) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_point() +
  labs(x = "log ratio of flapping per day", y = "log ground speed per day") +
  facet_wrap(~season+migration, nrow = 2)

# route straightness
clean_locations <- readRDS("/home/hbronnvik/Documents/chapter2/clean_data/clean_locations_2024-11-06.rds")

cl <- data.table::rbindlist(clean_locations, fill = T)

pride <- lapply(1:length(clean_locations), function(x){
  print(x)
  lgb <- clean_locations[[x]]
  # lgb_sf <- cl %>% filter(individual.id == 24563363) %>%
  #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)
  # mapview::mapview(lgb_sf, zcol = "distance")
  # the tracks
  rec <- records %>% 
    filter(individual.id == unique(lgb$individual.id))
  lgb <- lgb %>% 
    filter(date %in% rec$date)
  # ggplot(lgb, aes(location.long, location.lat, color = distance)) +
  #   geom_point()
  lgb %>% 
    left_join(rec, by = join_by(individual.id, date)) %>% 
    group_by(trackID) %>% 
    # for each track, the total distance between locations (down sampled to ~15 mins) in meters
    mutate(sum_dist_15 = sum(distance)) %>% 
    # just the first and last locations
    slice(1, n()) %>% 
    mutate(displacement = distVincentyEllipsoid(cbind(location.long[2], location.lat[2]), cbind(location.long[1], location.lat[1]))) %>% 
    slice(1) %>% 
    # straightness index D/L
    mutate(straightness = displacement/sum_dist_15) %>% 
    ungroup() %>% 
    dplyr::select(trackID, sum_dist_15, displacement, straightness)
}) %>% reduce(rbind)

pride <- pride %>% 
  left_join(meta %>% dplyr::select(trackID, individual.id, season, migration)) %>% 
  filter(individual.id %in% gps_birds$bird) 

flap_pride <- flapping %>% 
  mutate(flapping = behavior == "Flapping") %>% 
  group_by(trackID) %>% 
  summarize(flaps = sum(flapping),
            obs = n(),
            flap_ratio = flaps/obs) %>% 
  left_join(pride)

flap_pride %>% 
  ggplot(aes(as.factor(migration), straightness)) +
  geom_violin(aes(fill = migration)) +
  geom_boxplot(width = 0.33) +
  scale_fill_gradientn(colors = colfunc(5)) +
  facet_wrap(~season)

# create Figure S2
# png(filename = "//home/hbronnvik/Documents/chapter2/figures/look24/november/flapping_straight.png",
#     height = 8.5, width = 11, units = "in", res = 500)
flap_pride %>% 
  filter(migration < 5) %>% 
  mutate(label = paste0(season, " ", migration),
         season_cor = ifelse(season == "Spring", "Spring\nPearson's r = 0.373",
                             "Fall\nPearson's r = 0.168")) %>% 
  ggplot(aes(straightness, flap_ratio)) +
  geom_smooth(method = "lm", aes(color = label)) +
  geom_point() +
  scale_color_manual(values = c("#EE5E53", "#EF6E64", "#F17B6E", "#F59885",
                                "#0081A7", "#0095AF", "#00A9B7", "#54BDB8")) +
  scale_x_continuous(labels = c("0.4\nTortuous", "0.6", "0.8", "0.9\nStraight"), 
                     breaks = c(0.4, 0.6, 0.8, 0.9), limits = c(0.35, 0.95)) +
  labs(x = "Route straightness index", y = "Ratio of flapping to passive bursts", color = "Age") +
  facet_wrap(~season_cor) +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(face = "bold"))
# dev.off()

flap_pride %>% 
  filter(migration < 5) %>%
  group_by(season, migration) %>%
  summarize(cor = cor(straightness, flap_ratio))




