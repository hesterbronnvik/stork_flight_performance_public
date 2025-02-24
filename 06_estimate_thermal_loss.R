### Extract occasions when thermals were lost and compare across ages
### Hester Brønnvik
### 2023-11-13
### hbronnvik@ab.mpg.de

library(MASS)
library(tidyverse)
soar_gap <- 60
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 10), 
                  text = element_text(size = 15),
                  strip.background = element_rect(fill = "white", color = "white")))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

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
# antiwind_direction <- function(u, v){(270-(atan2(v, u)*(180/pi)))%%360 }

# from wip_seg.R
# fall_migrations <- readRDS("/home/hbronnvik/Documents/chapter2/migration_days_records_2023-11-02.rds")
# records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-05-07.rds")
# 361 individuals
# fall_migrations <- data.table::rbindlist(fall_migrations) %>% 
#   filter(season == "fall")
# the wind speed estimates for the thermals
# classified <- readRDS("/home/hbronnvik/Documents/chapter2/wind_thermal_data/full_track_seg_windspeed_2023-11-13.rds")

m_days <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2])) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day, track_status)

status <- m_days %>% 
  group_by(trackID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(trackID, track_status)

records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  mutate(season = ifelse(grepl("fall", trackID), "fall", "spring")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>%
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup() %>% 
  left_join(status) %>%
  dplyr::select(trackID, season, journey_number, track_status)


# thermal data with wind vectors attached
classified <- lapply(list.files("/home/hbronnvik/Documents/WS_data/wind_thermal_250deg_20s", 
                                full.names = T, pattern = "2024-11-08|2024-11-07"), readRDS)
classified <- classified[!sapply(classified, function(x) nrow(x)==0)]

# classified <- lapply(classified, function(ind){
#   ind <- ind %>% 
#     mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
#            time_diff = ifelse(is.na(time_diff), 1, time_diff),
#            min_split = time_diff > 60,
#            thermal_event = paste(ind_burst_id, cumsum(min_split)+1, sep = "_")) %>% 
#     dplyr::select(-time_diff, -min_split) %>% 
#     group_by(thermal_event) %>% 
#     mutate(thermal_duration = n(),
#            vspeed_thermal = (height_above_ellipsoid[n()]-height_above_ellipsoid[1])/thermal_duration,
#            turn_var_thermal = var(turn_angle, na.rm = T)/thermal_duration) %>% 
#     ungroup()
#   return(ind)
# })

classified <- lapply(classified, function(ind){
  ind <- ind %>%
    mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           time_diff = ifelse(is.na(time_diff), 1, time_diff),
           min_split = time_diff > 60,
           thermal_event = paste(ind_burst_id, cumsum(min_split)+1, sep = "_"),
           bearing = geosphere::bearing(p1 = cbind(location.long, location.lat)),
           bearing = ifelse(between(bearing, 0, 180), bearing, bearing+360),
           wind_speed = wind_speed(windX, windY),
           # north is 0 clockwise to 360 (the same as the heading from the tags)
           wind_direction = wind_direction(windX, windY),
           # rotate the wind compass 180 degrees simply to allow ease of comparison to Harel et al. 2016
           # antiwind = antiwind_direction(windX, windY),
           cross_wind = cross_wind(windX, windY, bearing),
           wind_support = wind_support(windX, windY, bearing)) %>%
    group_by(thermal_event) %>% 
    mutate(thermal_duration = n(),
           vspeed_thermal = (height_above_ellipsoid[n()]-height_above_ellipsoid[1])/thermal_duration,
           turn_sd_thermal = sd(turn_angle, na.rm = T)/thermal_duration,
           turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
           directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
           directional_set = paste(thermal_event, cumsum(c(F, directional_change[2:n()])), sep = "_"),
           n.changes = length(unique(directional_set)),
           n.change_s = n.changes/thermal_duration)  %>%
    ungroup() %>%
    left_join(records, by = join_by(trackID))
  ind <- ind %>%
    dplyr::select(timestamp, location.long, location.lat, height_above_ellipsoid, individual.id, burst_id, 
                  thermal_event, thermal_duration, vspeed_thermal, turn_sd_thermal, vert_speed_smooth,
                  vert_speed, wind_speed, wind_direction, bearing, windX, windY, CircRadius,
                  cross_wind, wind_support, turn_direction, directional_change, directional_set,
                  n.changes, n.change_s, season, journey_number, trackID, track_status)
  return(ind)
})

# get winds etc. from classified thermals to add on to these
# thermal data with wind vectors attached
classified <- data.table::rbindlist(classified, use.names = T, fill = T) 
classified <- classified[!is.na(classified$trackID),]
length(unique(classified$individual.id))

### 2. Pull out the 1 Hz data during migration
# from buzz_segmentation.R
seg_fls <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)

migrations <- lapply(seg_fls, function(x){
  ind <- readRDS(x)
  rec <- m_days %>% 
    filter(individual.id == unique(ind$individual.id))
  ind <- ind %>% 
    mutate(date = date(timestamp)) %>% 
    filter(date %in% rec$date) %>% 
    left_join(rec)
  gc()
  return(ind)
})
migrations <- migrations[!sapply(migrations, function(x) nrow(x)==0)]

### 1 Hz bursts that contain at least 30 seconds of thermal soaring & of another behavior
soar_bursts <- lapply(migrations, function(file){
  burst <- file %>% 
    mutate(burst_ind_ID = paste(individual.id, burst_id, sep = "_")) 
  tr <- burst %>% 
    filter(flight_clust_sm3 == "circular_soaring")  %>% 
    group_by(burst_ind_ID) %>%
    mutate(obs = n()) %>% 
    rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
    ungroup() %>% 
    filter(obs > 30)
  ot <- burst %>% 
    filter(flight_clust_sm3 != "circular_soaring")  %>% 
    group_by(burst_ind_ID) %>%
    mutate(obs = n()) %>% 
    rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
    ungroup() %>% 
    filter(obs > 30)
  if(nrow(tr) > 0 & nrow(ot) > 0){
    burst <- burst %>% 
      filter(burst_ind_ID %in% unique(tr$burst_ind_ID) & burst_ind_ID %in% unique(ot$burst_ind_ID))
  }
  if(nrow(burst)>0){
    saveRDS(burst, file = paste0("/home/hbronnvik/Documents/chapter2/soar_bursts/",
                                 unique(burst$individual.id), "_soar_bursts_",
                                 gsub("-", "", Sys.Date()), ".rds"))
  }
})
rm(migrations)
gc()
# 140 individuals left
soar_bursts <- soar_bursts[!sapply(soar_bursts, function(x) is.null(x))]
# soar_bursts <- data.table::rbindlist(soar_bursts, use.names = T, fill = T)
# saveRDS(soar_bursts, file = "/home/hbronnvik/Documents/chapter2/soar_bursts_20240223.rds")

fls <- list.files("/home/hbronnvik/Documents/chapter2/soar_bursts", pattern = "20241112", full.names = T)

gc()

start_time <- Sys.time()
thermal_other <- lapply(1:length(fls), function(n){
  # n <- 13
  b <- readRDS(fls[n])
  colnames(b)[which(colnames(b) == "azimuth_positive")] <- "heading"
  # b <- soar_bursts %>% filter(burst_ind_ID == "1173978476_6285")
  print(paste0("Extracting information on falling out for individual ", n, " of ", length(fls), "."))
  # for the purposes of fall out, do not count "shallow_circular_soaring" as different from "circular_soaring"
  
  # choose only the migratory dates
  record <- m_days %>%
    filter(trackID %in% unique(b$trackID))
  
  b <- b %>% 
    # mutate(flight_clust_sm3 = ifelse(flight_clust_sm3 == "shallow_circular_soaring", "circular_soaring", 
    #                                     flight_clust_sm3)) %>% 
    # fix the burst ID columns, MuffineJr. has the same burst number in multiple years
    mutate(burst_id = paste0(burst_id, "_", gsub("-", "", date)),
           track_flight_seg_id = paste0(track_flight_seg_id, "_", gsub("-", "", date)))
  
  eastern <- max(b$location.long) > 16
  
  if(nrow(b) > 0 & eastern == F){
    
    # tween_behav <- b %>% 
    #   arrange(timestamp) %>% 
    #   group_by(track_flight_seg_id) %>% 
    #   # duration of the behavior
    #   mutate(len = n()) %>% 
    #   slice(1) %>% 
    #   ungroup() %>% 
    #   # whether the behavior is between tight circular soars
    #   mutate(tween = lag(flight_clust_sm3) == "circular_soaring" & lead(flight_clust_sm3) == "circular_soaring") %>% 
    #   # only shallow, brief, between behaviors
    #   filter(flight_clust_sm3 == "shallow_circular_soaring" & len <= 60 & tween == T) %>% 
    #   # the info needed to re-classify
    #   dplyr::select(track_flight_seg_id, flight_clust_sm3)
    tween_behav <- b %>% 
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
    b <- b %>% 
      # smooth a 4th time
      mutate(flight_clust_sm3 = ifelse(track_flight_seg_id %in% tween_behav$track_flight_seg_id,
                                       "circular_soaring", flight_clust_sm3))
    
      b <- b %>% 
        group_by(burst_id) %>% 
        # find times when the bird changed classified flight behaviors
        mutate(behavioral_switch = flight_clust_sm3 != lag(flight_clust_sm3),
               behavioral_switch = ifelse(is.na(behavioral_switch), F, behavioral_switch),
               # signal whether that change was from or to circular soaring
               circle_switch = behavioral_switch == T & flight_clust_sm3 == "circular_soaring" | behavioral_switch == T & lag(flight_clust_sm3) == "circular_soaring",
               # number the different behaviors
               behavioral_bout = cumsum(circle_switch),
               behavioral_bout = behavioral_bout+1,
               # use the number as a label with the gps burst
               burst_bout = paste0(burst_id, "_", behavioral_bout),
               # count the number of times a bird switched to or from circling
               n.bouts = sum(circle_switch),
               # signal whether an observation is the start or end of a bout of behavior
               bout_pos = ifelse(circle_switch == T & flight_clust_sm3 == "circular_soaring", "after",
                                 ifelse(circle_switch == T & flight_clust_sm3 != "circular_soaring", "before", NA))) %>% 
        ungroup() %>% 
        # only proceed with bursts that contain two different behaviors
        filter(n.bouts > 2)
      
      if(nrow(b) > 0){

      fob_bouts <- b %>% 
        sf::st_drop_geometry() %>% 
        group_by(burst_bout) %>% 
        # take one observation per behavioral bout
        slice(1) %>% 
        ungroup() %>% 
        dplyr::select(timestamp, flight_clust_sm3, burst_id, burst_bout, behavioral_bout, bout_pos, circle_switch) %>% 
        group_by(burst_id) %>% 
        # label each behavioral bout as either being between two circular soaring bouts or not
        mutate(middle = flight_clust_sm3 != "circular_soaring" & 
                 lag(flight_clust_sm3) == "circular_soaring" &
                 lead(flight_clust_sm3) == "circular_soaring") %>% 
        ungroup() %>% 
        # proceed with only the bouts between two circular soaring bouts
        filter(middle == T) %>% 
        dplyr::select(burst_bout, behavioral_bout) 
      
      # b %>%
      #   filter(burst_id %in% 5132) %>% #unique(burst_id)[1:11]) %>%
      #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs =  sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>%
      #   mapview::mapview(zcol = "behavioral_bout", burst = F)
      
      classified_df <- classified %>% 
        # only this bird's wind-classified data
        filter(individual.id == unique(b$individual.id)) %>%
        mutate(burst_ind_ID = paste(individual.id, burst_id, sep = "_")) %>% 
        # only wind-classified data relevant for the bouts between two circular soaring bouts 
        filter(burst_ind_ID %in% unique(b$burst_ind_ID)) %>% 
        dplyr::select(timestamp, thermal_duration, vspeed_thermal, turn_sd_thermal, wind_speed, 
                      wind_direction, cross_wind, wind_support, turn_direction, directional_change, 
                      directional_set, n.changes, n.change_s, season, journey_number, thermal_event)
      
      b <- b %>% 
        left_join(classified_df, by = join_by(timestamp))
      # class_sf <- b %>% 
      #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = "EPSG:4326")
      # mapview::mapview(class_sf, zcol = "flight_clust_sm3")
      
      exits <- b %>% 
        # the behaviors between circular soaring
        filter(burst_bout %in% fob_bouts$burst_bout) %>% 
        group_by(burst_bout) %>% 
        # the duration of each behavioral bout between circular soaring
        mutate(time_spent = as.numeric(difftime(timestamp[n()], timestamp[1], units = "secs")),
               # and the type of behavior
               out_behav = paste(unique(flight_clust_sm3), collapse = " "))  %>% 
        slice(1) %>% 
        ungroup()
      
      if(nrow(exits) > 0){
        # go through each burst separately (slow and steady)
        exits <- exits %>% 
          group_by(burst_id) %>% 
          group_split()
        
        soaring_metrics <- lapply(exits, function(exit){
            # go through each behavior separately
          burst_soaring_metrics <- lapply(unique(exit$behavioral_bout), function(event){
            # Here, we are only interested in behaviors when we are confident in their classification,
            # if there is only one second at the start of a burst, we cannot know what was happening
            pre_exit <- b %>% filter(behavioral_bout == event-1 & burst_id == unique(exit$burst_id)) 
            out <- b %>% filter(behavioral_bout == event & burst_id == unique(exit$burst_id))
            post_exit <- b %>% filter(behavioral_bout == event+1 & burst_id == unique(exit$burst_id)) 
            
            if(nrow(pre_exit) > 29 & nrow(post_exit) > 29 & nrow(out) > 9){
              pre_exit <- pre_exit %>% 
                # count the number of directional changes in the first circular bout
                mutate(n.changes = sum(directional_change, na.rm = T),
                       # label each burst's bout's direction as a subsublabel
                       directional_set = cumsum(c(F, directional_change[2:n()])))
              # True/False the bird changed direction in the 20 seconds before leaving the thermal
              does_turn <- pre_exit %>% 
                slice((n()-19):n()) %>% 
                filter(directional_change == T)
              
              # if there are wind annotations in the last 5 seconds of the thermal,
              if(max(which(!is.na(pre_exit$wind_direction))) > (nrow(pre_exit)-6)){
                # take the last direction of wind
                exit_wind_direction <- pre_exit %>% 
                  sf::st_drop_geometry() %>% 
                  drop_na(wind_direction) %>% 
                  slice(n()) %>% 
                  dplyr::select(wind_direction) %>% 
                  deframe()
                # take the last direction of the bird
                exit_heading <- pre_exit%>% 
                  sf::st_drop_geometry() %>% 
                  drop_na(wind_direction) %>% 
                  slice(n()) %>% 
                  dplyr::select(heading) %>% 
                  deframe()
              }else{ 
        # if there is no wind value in the last five seconds of the thermal it fell out of
                exit_wind_direction <- NA
                exit_heading <- NA
              }
              # for the circling after the behavioral change
              post_exit <- post_exit %>% 
                # count the number of directional changes in the first circular bout
                mutate(n.changes = sum(directional_change, na.rm = T),
                       # label each burst's bout's direction as a subsublabel
                       directional_set = cumsum(c(F, directional_change[2:n()])))
              # True/False the bird changed direction in the 10 seconds after entering the thermal
              turn_again <- post_exit %>% 
                slice(1:20) %>% 
                filter(directional_change == T)
              # collate information about each behavioral bout in each burst
              metrics <- data.frame(individual.id = unique(b$individual.id),
                                    trackID = unique(exit$trackID),
                                    burst_ind_ID = unique(exit$burst_ind_ID),
                                    # the latitude where the bird lost a thermal
                                    location.lat = pre_exit$location.lat[nrow(pre_exit)],
                                    # the number identifying the behavioral bout between two circling bouts
                                    exit = event,
                                    date = unique(date(exit$timestamp)),
                                    # the track ID
                                    track = unique(m_days$trackID[m_days$trackID == unique(exit$trackID)]),
                                    # the migration day
                                    ld_day = unique(m_days$ld_day[m_days$date == unique(date(exit$timestamp)) & 
                                                                    m_days$trackID == unique(exit$trackID)]), 
                                    # the final vertical speed before fall out
                                    exit_vspeed = pre_exit$vert.speed[nrow(pre_exit)], 
                                    # the first vertical speed after fall out
                                    reentry_vspeed = post_exit$vert.speed[1], 
                                    # the vertical speed of the preceding thermal segment
                                    pre_vspeed = (pre_exit$height.above.ellipsoid[nrow(pre_exit)]-pre_exit$height.above.ellipsoid[1])/nrow(pre_exit),
                                    # the vertical speed of the following thermal segment
                                    post_vspeed = (post_exit$height.above.ellipsoid[nrow(post_exit)]-post_exit$height.above.ellipsoid[1])/nrow(post_exit),
                                    # last recorded height before fall out
                                    exit_height = pre_exit$height.above.ellipsoid[nrow(pre_exit)], 
                                    # first recorded height of re-entry
                                    reentry_height = post_exit$height.above.ellipsoid[1], 
                                    # which behaviors were performed between thermals
                                    out_behav = gsub("soaring", "", paste(unique(out$flight_clust_sm3), collapse = " ")), 
                                    # how long the bird was in the behavioral bout between soaring bouts
                                    out_secs = exit$time_spent[which(exit$behavioral_bout == event)],
                                    # how long the bird spent soaring before fall out
                                    pre_duration = as.numeric(difftime(pre_exit$timestamp[nrow(pre_exit)],
                                                                       pre_exit$timestamp[1],
                                                                       units = "secs")), 
                                    # how long the bird spent soaring after fall out
                                    post_duration = as.numeric(difftime(post_exit$timestamp[nrow(post_exit)],
                                                                        post_exit$timestamp[1],
                                                                        units = "secs")), 
                                    # the number of times the bird switched direction before fall out
                                    turns_pre = unique(pre_exit$n.changes), 
                                    # the number of times the bird switched direction after fall out
                                    turns_post = unique(post_exit$n.changes), 
                                    # whether or not the bird turned in the 10 seconds before fall out
                                    exit_turn = nrow(does_turn)>0, 
                                    # whether or not the bird turned in the 10 seconds after fall out
                                    reentry_turn = nrow(turn_again)>0, 
                                    # the mean value of however many wind estimates were available before fall out
                                    avg_wind_speed_pre = mean(pre_exit$wind_speed, na.rm = T), 
                                    # the mean of the wind after fall out
                                    avg_wind_speed_post = mean(post_exit$wind_speed, na.rm = T), 
                                    # the mean value of however many wind estimates were available before fall out
                                    avg_support_pre = mean(pre_exit$wind_support, na.rm = T), 
                                    # the mean of the wind after fall out
                                    avg_cross_pre = mean(pre_exit$cross_wind, na.rm = T), 
                                    # the variance in all of the wind estimates before fall out
                                    var_wind_speed_pre = var(pre_exit$wind_speed, na.rm = T), 
                                    # the variance in all of the wind estimates after fall out
                                    var_wind_speed_post = var(post_exit$wind_speed, na.rm = T), 
                                    var_wind_direction_pre = var(pre_exit$wind_direction, na.rm = T),
                                    var_wind_direction_post = var(post_exit$wind_direction, na.rm = T),
                                    # the direction of the last wind estimate within five seconds before fall out
                                    exit_wind_direction = exit_wind_direction, 
                                    # the heading of the bird at the same point as the last wind estimate within five seconds before fall out
                                    exit_heading = exit_heading, 
                                    # the total height reached after the second soaring bout
                                    final_height = post_exit$height.above.ellipsoid[nrow(post_exit)],
                                    # the thermal ID after the thermal loss
                                    post_thermal_id = unique(post_exit$thermal_event),
                                    # the thermal ID before the thermal loss
                                    pre_thermal_id = unique(pre_exit$thermal_event)) 
              return(metrics)
            }
          }) %>% reduce(rbind)
          return(burst_soaring_metrics)
        }) %>% reduce(rbind)
        return(soaring_metrics)
      }
    }
  }
})
Sys.time()-start_time # Time difference of 26.86288 mins

thermal_other <- thermal_other[!sapply(thermal_other, function(x) is.null(x))]
thermal_other <- thermal_other[!sapply(thermal_other, function(x) nrow(x)==0)]
# 130 IDs meeting the criteria
thermal_other <- data.table::rbindlist(thermal_other, use.names = T)

# saveRDS(thermal_other, file = paste0("/home/hbronnvik/Documents/chapter2/between_thermals_", gsub("-", "", Sys.Date()), ".rds"))
thermal_other <- readRDS("/home/hbronnvik/Documents/chapter2/between_thermals_20241112.rds")

thermal_other  %>% 
  filter(between(out_secs, 10, 60)) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  mutate(exit_sink = exit_height-reentry_height) %>% 
  ggplot(aes(exit_sink)) +
  geom_histogram(color = "black") + 
  labs(x = "Exit height - re-entry height", y = "Count") +
  facet_wrap(~season+journey_number, scales = "free_y", nrow = 2)

thermal_count <- classified %>%
  mutate(date = date(timestamp)) %>% 
  group_by(individual.id, date) %>% 
  mutate(n_thermals = length(unique(thermal_event))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(date, individual.id, n_thermals)

eastern_IDs <- classified %>% 
  filter(location.long >= 16) %>% 
  dplyr::select(individual.id) %>% 
  deframe() %>% 
  unique()

thermal_other <- thermal_other  %>% 
  filter(between(out_secs, 10, 60) & !individual.id %in% eastern_IDs) %>%
  group_by(individual.id, date) %>% 
  mutate(n_falls = length(unique(paste0(burst_ind_ID, "_", exit)))) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(thermal_count, by = join_by(individual.id, date)) %>% 
  mutate(fall_in_ratio = n_falls/n_thermals) %>% 
  filter(journey_number < 5) %>%
  ungroup() %>% 
  mutate(label = paste0(journey_number,
                        ifelse(journey_number == 1, "st ", 
                               ifelse(journey_number == 2, "nd ",
                                      ifelse(journey_number == 3, "rd ", "th "))),
                        sub("fall", "Fall", 
                            sub("spring", "Spring", season))))

# look at the final heights of the whole thermals that the birds fell out of over time
ms <- classified %>% 
  filter(thermal_event %in% unique(thermal_other$post_thermal_id) & journey_number < 5) %>% 
  group_by(thermal_event) %>% 
  slice(n()) %>% 
  mutate(label = paste0(journey_number,
                        ifelse(journey_number == 1, "st ", 
                               ifelse(journey_number == 2, "nd ",
                                      ifelse(journey_number == 3, "rd ", "th "))),
                        sub("fall", "Fall", 
                            sub("spring", "Spring", season)))) %>% 
  group_by(label) %>% 
  summarize(m = median(height_above_ellipsoid)) %>% 
  dplyr::select(m) %>% 
  deframe()
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/final_heights.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
classified %>% 
  filter(thermal_event %in% unique(thermal_other$post_thermal_id)) %>% 
  group_by(thermal_event) %>% 
  slice(n()) %>% 
  mutate(label = paste0(journey_number,
                        ifelse(journey_number == 1, "st ", 
                               ifelse(journey_number == 2, "nd ",
                                      ifelse(journey_number == 3, "rd ", "th "))),
                        sub("fall", "Fall", 
                            sub("spring", "Spring", season)))) %>% 
  dplyr::select(label, height_above_ellipsoid) %>% 
  ggplot(aes(label, height_above_ellipsoid)) +
  geom_hline(yintercept = ms, lty = 2, color = colfunc(8), show.legend = T, lwd = 1) +
  geom_violin(aes(fill = label)) +
  geom_boxplot(width = 0.33) +
  # facet_wrap(~season) +
  scale_fill_manual(values = colfunc(8)) +
  labs(x = "Age", y = "Maximum height reached in a thermal (m)") +
  theme(legend.position = "none")
# dev.off()


### find the thermals that did not have losses
ridden <- classified %>% 
  filter(thermal_duration >= 30 & vspeed_thermal > 0 & !individual.id %in% eastern_IDs) %>% 
  group_by(thermal_event) %>% 
  mutate(avg_wind_sp = mean(wind_speed, na.rm = T),
         avg_vspeed = mean(vert_speed_smooth, na.rm = T),
         sd_wind_sp = sd(wind_speed, na.rm = T),
         # if the thermal had a loss, 1, else 0
         ridden = ifelse(!thermal_event %in% unique(c(thermal_other$pre_thermal_id, thermal_other$post_thermal_id)), 0, 1)) %>% 
  slice(1) %>%
  ungroup() %>% 
  # focus on sufficient sample sizes
  filter(journey_number < 5) %>% 
  # prefer realistic values
  filter(between(avg_vspeed, quantile(avg_vspeed, 0.025), quantile(avg_vspeed, 0.975))) %>% 
  # drop useless thermals
  drop_na(sd_wind_sp) 

x <- ridden$sd_wind_sp
b <- boxcox(lm(x ~ 1))
# Exact lambda
lambda <- b$x[which.max(b$y)]
ridden$sd_wind_sp <- (x ^ lambda - 1) / lambda

x <- ridden$avg_vspeed+abs(min(ridden$avg_vspeed))+0.1
b <- boxcox(lm(x ~ 1))
# Exact lambda
lambda <- b$x[which.max(b$y)]
ridden$avg_vspeed <- (x ^ lambda - 1) / lambda

# table S1 The number of individuals that contributed data in each migration
ridden %>% 
  group_by(journey_number, season) %>% 
  summarize(len = length(unique(individual.id)))
# there are 151 total storks, but only 144 first falls
length(unique(ridden$individual.id))
# which birds did not contribute a first fall?
juvies <- unique(ridden$individual.id[which(ridden$journey_number == 1 & ridden$season == "fall")])
ridden %>% 
  filter(!individual.id %in% juvies) %>% 
  group_by(journey_number, season) %>% 
  reframe(individual.id = unique(individual.id)) %>% 
  arrange(individual.id)

# the total number of thermals detected
length(unique(ridden$thermal_event))
# the number of thermal losses detected
length(unique(ridden$thermal_event[ridden$ridden == T]))

# save out the losses and non-losses for efficiency estimation in script 7
# saveRDS(ridden, file = "/home/hbronnvik/Documents/chapter2/thermal_losses_20241130.rds")

fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("fall", "spring")
c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53")
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/wind_vspeed_loss.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
pls <- lapply(c("spring", "fall"), function(s){
  lab <- ifelse(s == "spring", "Spring", "Fall")
  if(s == "spring"){
    lab2 <- c("C", "D")
  }else{
    lab2 <- c("A", "B")
  }
  if(s == "spring"){
    cols <- c("#0081A7", "#7fc4b8")
  }else{
    cols <- c("#ee5e53",  "#f7a58f")
  }
  ps <- lapply(c("sd_wind_sp", "avg_vspeed"), function(v){
    p <- ridden %>% 
      filter(season == s) %>% 
      dplyr::select(ridden, journey_number, sd_wind_sp, avg_vspeed) %>% 
      # mutate(avg_vspeed = sqrt(avg_vspeed),
      #        # avg_wind_sp = sqrt(avg_wind_sp)
      #        ) %>%
      pivot_longer(cols = contains("sp")) %>% 
      filter(name == v) %>% 
      mutate(name = ifelse(name == "sd_wind_sp", "sd wind speed", "mean thermal strength"),
             ridden = ifelse(ridden == 0, "No loss", "Loss")) %>% 
      ggplot(aes(value, as.character(journey_number), fill = ridden)) +
      # ggridges::geom_density_ridges(panel_scaling = T, alpha = 0.75) +
      ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                    panel_scaling = T, alpha = 0.85) +
      scale_fill_manual(values = cols, name = "") +
      labs(x = "m/s", y = "Age", title = "") +
      facet_wrap(~name, scales = "free", nrow = 1)
    return(p)
  })
  pl <- ggpubr::ggarrange(ps[[2]], ps[[1]], nrow = 1, labels = lab2, 
                          common.legend = T, legend = "top")
  return(pl)
})
loss_dists <- ggpubr::ggarrange(pls[[2]], pls[[1]], nrow = 2)

# scale the predictors
ridden <- ridden %>% 
  mutate(individual.id = as.factor(individual.id)) %>% 
  dplyr::select(individual.id, season, thermal_event, ridden, journey_number, 
                n.change_s, sd_wind_sp, avg_vspeed, avg_wind_sp) %>%
  mutate_at(c("journey_number", "n.change_s", "sd_wind_sp", "avg_vspeed", "avg_wind_sp"), 
            list(z = ~(scale(.)[,1])))

# leaves 49856 thermals without losses and 7150 thermals with
ridden %>% 
  group_by(ridden) %>% 
  summarize(obs = length(unique(thermal_event)))

# thermals they lost had slightly lower vertical speeds and more directional changes
ridden %>% 
  pivot_longer(cols = c("n.change_s_z", "sd_wind_sp_z", "avg_vspeed_z")) %>% 
  ggplot(aes(value, fill = as.character(ridden), group = as.character(ridden))) +
  geom_density(alpha = 0.5) +
  labs(x = "Scaled value", y = "Density", fill = "Thermal loss") +
  facet_wrap(~name+season, ncol = 2)

pl <- lapply(c("spring", "fall"), function(s){
  col <- ifelse(s == "spring", "#007BA7", "#F07268")
  lab <- ifelse(s == "spring", "Spring", "Fall")
  # plot bars
  prop_losses <- ridden %>% 
    filter(season == s) %>% 
    group_by(ridden, journey_number) %>% 
    summarise(obs = n()) %>% 
    group_by(journey_number) %>% 
    mutate(total = sum(obs),
           prop = obs/total) %>% 
    filter(ridden == 1) %>% 
    mutate(#label = paste0(journey_number,
                          # ifelse(journey_number == 1, "st ", 
                          #        ifelse(journey_number == 2, "nd ",
                          #               ifelse(journey_number == 3, "rd ", "th "))),
                          # s,
                          # "\n(n = ", total, ")"),
           label = paste0(journey_number, "\n(n=", total, ")"),
           prop = prop*100) %>% 
    ggplot(aes(x = label, y = prop)) + 
    geom_hline(yintercept = seq(0, 14, 2.5), color = "gray50") +
    geom_bar(position="stack", stat="identity", width = 0.5, fill = col) +
    # scale_x_discrete(expand = c(0, 0 )) +
    scale_y_continuous(expand = c(0, 0 ), limits = c(0, 16))+
    geom_text(aes(label=paste0(sprintf("%1.1f", prop),"%")), y = 14.75,
              # position=position_stack(vjust=0.5), size = 3.5,
              color = "black", fontface = "bold") +
    # scale_fill_manual(values = c("#EE5E53", "#0081A7", "#F07268", "#009BB1", "#F38979", "#24B5B8", "#FABBA0", "#B5CDB7")) +
    labs(x = "Age", y = "\nPercentage of\nthermals lost", title = "") +
    theme(legend.position = "none",
          # plot.margin = margin(5,50,5,50)
          )
  return(prop_losses)
})
ggpubr::ggarrange(pl[[2]], pl[[1]], nrow = 2)

effects <- lapply(c("spring", "fall"), function(s){
  ride_mod <- lme4::lmer(ridden~sd_wind_sp_z+avg_vspeed_z+journey_number_z+(1|individual.id), 
                         data = ridden[ridden$season == s,])
  # lab <- round(performance::performance_rmse(ride_mod), 2)
  col <- ifelse(s == "spring", "#007BA7", "#F07268")
  
  conf_limit_ride <- confint(ride_mod) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "pred") %>%
    left_join(lme4::fixef(ride_mod) %>% 
                as.data.frame() %>% 
                rename(fixef = ".") %>% 
                rownames_to_column(var = "pred"), 
              by = join_by(pred)) %>% 
    filter(!grepl("sig|cept", pred)) %>% 
    rename(lower = "2.5 %",
           upper = "97.5 %") %>% 
    mutate(pred = gsub("journey_number", "Age",
                      sub("avg_vspeed", "mean thermal strength",
                          sub("sd_wind_sp", "sd wind speed",
                              gsub("_z", "", pred)))),
           sig = sign(lower)==sign(upper))
  
  ride_coefs <- ggplot(conf_limit_ride, aes(fixef, fct_reorder(pred, nchar(pred)))) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_linerange(aes(xmin = lower, xmax = upper), lwd = 1.75, color = col) +
    geom_point(size = 4, color = col) +
    scale_x_continuous(limits = c(-0.05, 0.08), breaks = seq(-0.05, 0.1, 0.025)) +
    labs(x = "Fixed effect estimate and 95% CI", y = "", color = "Significant", title = "")
  return(ride_coefs)
})

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/loss_int3_coefs.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(effects[[2]]+labs(title = "Fall"), 
                  effects[[1]]+labs(title = "Spring"), 
                  nrow = 2, labels = "AUTO")
# dev.off()

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/loss_wind_vspeed_coefs.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(pl[[2]], effects[[2]], 
                  pl[[1]], effects[[1]], 
                  nrow = 2, ncol = 2, labels = c("A", "C", "B", "D"))
# dev.off()

# create Figure 2
tgrob1 <- ggpubr::text_grob("Fall", size = 20)
tgrob2 <- ggpubr::text_grob("Spring", size = 20)
# Draw the text
plot_0 <- ggpubr::as_ggplot(tgrob1) + theme(plot.margin = margin(0,12,0,0, "cm"))
plot_1 <- ggpubr::as_ggplot(tgrob2) + theme(plot.margin = margin(0,12,0,0, "cm"))


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/loss_dists_coefs.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(plot_0, plot_1,
                  pls[[2]]+theme(plot.margin = margin(0,0,0,0)), 
                  pls[[1]]+theme(plot.margin = margin(0,0,0,0)),
                  ggpubr::ggarrange(pl[[2]], effects[[2]], labels = c("E", "G"), nrow = 2), 
                  ggpubr::ggarrange(pl[[1]], effects[[1]], labels = c("F", "H"), nrow = 2), 
                  nrow = 3, ncol = 2, heights = c(1, 5, 8), vjust = 6)
# dev.off()

# create Figure S1
p_therms <- classified %>% 
  filter(thermal_duration >= 30 & vspeed_thermal > 0 & !individual.id %in% eastern_IDs) %>% 
  group_by(thermal_event) %>% 
  mutate(avg_wind_sp = mean(wind_speed, na.rm = T),
         avg_vspeed = mean(vert_speed_smooth, na.rm = T),
         avg_rad = mean(CircRadius, na.rm = T),
         sd_wind_sp = sd(wind_speed, na.rm = T),
         # if the thermal had a loss, 1, else 0
         ridden = ifelse(!thermal_event %in% unique(c(thermal_other$pre_thermal_id, thermal_other$post_thermal_id)), 0, 1)) %>% 
  slice(1) %>%
  ungroup() %>%
  # focus on sufficient sample sizes
  filter(journey_number < 5) %>% 
  # prefer realistic values
  filter(between(avg_vspeed, quantile(avg_vspeed, 0.025), quantile(avg_vspeed, 0.975))) %>% 
  # drop useless thermals
  drop_na(sd_wind_sp)  %>% 
  mutate(log_turn = log(turn_sd_thermal),
         sqrt_vspeed = sqrt(vspeed_thermal),
         sqrt_ws = sqrt(avg_wind_sp),
         sqrt_sd_ws = sqrt(sd_wind_sp),
         log_rad = log(avg_rad),
         season = ifelse(season == "fall", "Fall", "Spring"),
         age = paste0(season, " ", journey_number)) %>% 
  group_by(age) %>% 
  mutate(samp_size = n()) %>% 
  ungroup() %>% 
  dplyr::select(journey_number, season, age, samp_size,
                sqrt_vspeed, log_turn, sqrt_ws, sqrt_sd_ws, log_rad) %>% 
  pivot_longer(cols = sqrt_vspeed:log_rad) 

fac_labs <- c("sqrt vertical speed (m/s)", "log turning angle variance (deg/s)", "log mean circling radius (m)", 
              "sqrt mean wind speed (m/s)", "sqrt sd wind speed (m/s)")
names(fac_labs) <- c("sqrt_vspeed", "log_turn", "log_rad", "sqrt_ws", "sqrt_sd_ws")

info <- p_therms %>% 
  group_by(season, journey_number, name) %>% 
  summarize(Mean = round(mean(value, na.rm = T), 2),
            Median = round(median(value, na.rm = T), 2),
            SD = round(sd(value, na.rm = T), 2)) %>% 
  ungroup() %>% 
  rename(Age = journey_number)

leg1 <- ggpubr::get_legend(p_therms %>%
  filter(season == "Fall") %>% 
  mutate(age = paste0(age, ", n = (", samp_size, ")")) %>% 
  ggplot(aes(value, fill = age)) +
  geom_density(alpha = 0.75) + 
  labs(x = "", y = "", fill = "Age") +
  scale_fill_manual(values = c("#EE5E53", "#F07268", "#F38979", "#FABBA0")) +
  facet_wrap(~name+season, ncol = 2, scales = "free", labeller = labeller(name = fac_labs)) 
)
leg2 <- ggpubr::get_legend(p_therms %>%
                             filter(season == "Spring") %>% 
                             mutate(age = paste0(age, ", n = (", samp_size, ")")) %>% 
                             ggplot(aes(value, fill = age)) +
                             geom_density(alpha = 0.75) + 
                             labs(x = "", y = "", fill = "Age") +
                             scale_fill_manual(values = c("#0081A7", "#009BB1", "#24B5B8", "#B5CDB7")) +
                             facet_wrap(~name+season, ncol = 2, scales = "free", labeller = labeller(name = fac_labs)) 
)

panels <- lapply(split(p_therms, list(p_therms$name, p_therms$season)), function(p){
  if(unique(p$season) == "Fall"){
    cols <- c("#EE5E53", "#F07268", "#F38979", "#FABBA0")
  }else{
    cols <- c("#0081A7", "#009BB1", "#24B5B8", "#B5CDB7")
  }
  st <- ifelse(unique(p$name) == "sqrt_vspeed",
               "sqrt vertical speed (m/s)",
               ifelse(unique(p$name) == "log_turn",
                      "log turning angle SD (deg/s)",
                      ifelse(unique(p$name) == "log_rad",
                             "log mean circling radius (m)",
                             ifelse(unique(p$name) == "sqrt_ws",
                                    "sqrt mean wind speed (m/s)",
                                    ifelse(unique(p$name) == "sqrt_sd_ws",
                                           "sqrt SD wind speed (m/s)", NA)))))
  tab <- info %>% 
    filter(season == unique(p$season) & name == unique(p$name)) %>% 
    dplyr::select(Age, Mean, Median, SD)
  pl <- ggplot(p, aes(value, fill = age)) +
    geom_density(alpha = 0.75) + 
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0)) +
    # annotation_custom(gridExtra::tableGrob(tab, rows=NULL, 
    #                                        theme = gridExtra::ttheme_minimal()), 
    #                   xmin=4, xmax=5, ymin=1, ymax=2) +
    labs(x = "", y = "", fill = "Age", subtitle = st) +
    theme(legend.position = "none")
  pl2 <- ggpubr::ggarrange(pl, gridExtra::tableGrob(tab, rows=NULL, 
                                             theme = gridExtra::ttheme_minimal()), 
                           widths = c(2, 1))
  return(pl2)
})

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/densities_full.png",
#     height = 8.5*1.5, width = 11*1.5, units = "in", res = 500)
ggpubr::ggarrange(leg1, leg2,
                  panels[[1]], panels[[6]], 
                  panels[[2]], panels[[7]], 
                  panels[[3]], panels[[8]],
                  panels[[4]], panels[[9]], 
                  panels[[5]], panels[[10]], 
                  ncol = 2, nrow = 6)
# dev.off()

# create Table S2
# model the affect of age on each performance metric
m_therms <- classified %>% 
  filter(thermal_duration >= 30 & vspeed_thermal > 0 & !individual.id %in% eastern_IDs) %>% 
  group_by(thermal_event) %>% 
  mutate(avg_wind_sp = mean(wind_speed, na.rm = T),
         avg_vspeed = mean(vert_speed_smooth, na.rm = T),
         avg_rad = mean(CircRadius, na.rm = T),
         sd_wind_sp = sd(wind_speed, na.rm = T),
         # if the thermal had a loss, 1, else 0
         ridden = ifelse(!thermal_event %in% unique(c(thermal_other$pre_thermal_id, thermal_other$post_thermal_id)), 0, 1)) %>% 
  slice(1) %>%
  ungroup() %>% 
  # focus on sufficient sample sizes
  filter(journey_number < 5) %>% 
  mutate(log_turn = log(turn_sd_thermal),
         sqrt_vspeed = sqrt(vspeed_thermal),
         sqrt_ws = sqrt(avg_wind_sp),
         sqrt_sd_ws = sqrt(sd_wind_sp),
         log_rad = log(avg_rad),
         season = ifelse(season == "fall", "Fall", "Spring")) %>% 
  # prefer realistic values
  filter(between(avg_vspeed, quantile(avg_vspeed, 0.025), quantile(avg_vspeed, 0.975))) %>%
  # drop useless thermals
  drop_na(sd_wind_sp) %>% 
  dplyr::select(individual.id, season, journey_number, sqrt_vspeed, log_turn, sqrt_ws, sqrt_sd_ws, log_rad) %>% 
  mutate(ID = as.factor(individual.id),
         age_z = scale(journey_number)[,1],
         sqrt_vspeed_z = scale(sqrt_vspeed)[,1],
         log_turn_z = scale(log_turn)[,1],
         sqrt_ws_z = scale(sqrt_ws)[,1],
         sqrt_sd_ws_z = scale(sqrt_sd_ws)[,1],
         log_rad_z = scale(log_rad)[,1])

metrics_table <- lapply(c("sqrt_vspeed", "log_turn", "sqrt_ws", "sqrt_sd_ws", "log_rad"), function(v){
  ests <- lapply(c("Fall", "Spring"), function(s){
    # the relevant data
    temp <- m_therms %>% 
      dplyr::select(v, age_z, season, ID) %>% 
      filter(season == s) %>% 
      rename(var_oi = v)
    # the model
    temp_mod <- lme4::lmer(var_oi~age_z+(1|ID), data = temp)
    # the coefficients
    temp_ests <- as.data.frame(summary(temp_mod)$coefficients)
    # the confidence interval (95%)
    temp_ests$lower <- as.data.frame(confint(temp_mod, method = "Wald"))[4, 1]
    temp_ests$upper <- as.data.frame(confint(temp_mod, method = "Wald"))[4, 2]
    # the variable
    temp_ests$variable <- v
    # the season
    temp_ests$season <- s
    # the performance
    temp_ests$rmse <- performance::rmse(temp_mod)
    # clean up
    temp_ests <- temp_ests[2,]
    return(temp_ests)
  }) %>% reduce(rbind)
  return(ests)
}) %>% reduce(rbind)


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/vspeed_losses_rain.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
thermal_other %>% 
  ungroup()%>%
  # mutate(label = factor(label, levels = c("Pre-migration", unique(thermal_other$label)[1:6]))) %>% 
  dplyr::select(label, pre_vspeed, post_vspeed) %>% 
  pivot_longer(cols = contains("vspeed")) %>% 
  ggplot(aes(y = fct_rev(label), x = value, fill = name)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    # width = .6, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0
  ) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA ## `outlier.shape = NA` or `outlier.alpha = 0` works as well
  ) +
  ## add dot plots from {ggdist} package
  # ggdist::stat_dots(
  #   ## orientation to the left
  #   side = "left", 
  #   ## move geom to the left
  #   justification = 1.12
  # ) +
  scale_fill_manual(values = colfunc(4), labels = c("After", "Before")) +
  scale_color_manual(values = colfunc(4), labels = c("After", "Before")) +
  labs(y = " ", x = "Vertical speed (m/s)", fill = "Relative to\nthermal loss")
# dev.off()
