### Identify birds to include, segment their migration tracks, save dates out
### Hester Brønnvik
### 29.09.2022
### hbronnvik@ab.mpg.de

library(lubridate)
library(geosphere)
library(move)
library(stringr)
library(tidyverse)
library(data.table)
library(mapview)
d_thresh <- 40000 # meters
w_thresh <- 6 # weeks
s_thresh <- 30 # days
l_thresh <- 250 # km
g_thresh <- 7 # days (allowable gap in transmission)

# required information
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData") # Movebank credentials
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633) # Movebank study IDs

# Here, birds to include are the ones that migrated at all. We identify successful migrations later, 
# but whether the individual generated 1 Hz data is not considered here.

# get the current names and deployment dates
metad <- lapply(studies, function(x){
  md <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653) %>%
    dplyr::select(animal_id, deploy_on_timestamp, animal_local_identifier, tag_timestamp_end) %>%
    mutate(deploy_on_timestamp = date(deploy_on_timestamp),study = x) %>%
    rename(individual.id = animal_id)
  return(md)
}) %>% 
  reduce(rbind) %>% 
  # remove birds that were captured twice
  filter(!duplicated(individual.id)) %>% 
  mutate(animal_local_identifier = gsub("\\+ |\\- ", "", animal_local_identifier)) 

# dates retrieved by AF + HB visual estimates 10-2024
deaths2 <- read.csv("/home/hbronnvik/Documents/chapter2/death_dates_1024.csv") %>% 
  # combine estimates
  mutate(dod_diff = as.numeric(difftime(visDOD, dod, units = "days")),
         # if there is no official death date, but HB recorded a death, use the visual estimate
         dod = ifelse(is.na(dod) & !is.na(visDOD), visDOD, 
                      # if there are two dates, check the difference, if they are different, use the HB estimates
                      # if there is little difference or HB did not record death, trust the metadata
                      ifelse(!is.na(visDOD) & !is.na(dod) & abs(dod_diff) > 3, date(visDOD), dod)),
         dod = as.Date(dod, tz = "UTC")) %>% 
  dplyr::select(individual.id, death_animal_identifier, dod)

# HB visual estimates 09-2023
visDODs <- readRDS("/home/hbronnvik/Documents/storkSSFs/visual_estimated_deaths_2023-09-05.rds") %>% 
  rename(death_animal_identifier = local_identifier,
         individual.id = individual_id) %>% 
  dplyr::select(individual.id, death_animal_identifier, dod)

deaths <- visDODs %>% 
  filter(!individual.id %in% unique(deaths2$individual.id)) %>% 
  rbind(deaths2) %>% 
  drop_na(dod)

# full metadata for 526 individuals that sent GPS data
ref_info <- lapply(studies, function(i){
  birds <- getMovebankReferenceTable(study = i, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653) %>%
    mutate(study = i)
  birds
})
ref_info <- data.table::rbindlist(ref_info, fill = T) %>%
  rename(individual.id = animal_id)

# determine the identities of the nestling birds (remove any care center adults)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653)
  if("animal_life_stage" %in% colnames(birds)){
    chicks <- birds %>% 
      filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
  }
  return(chicks)
}) %>% reduce(rbind)

# here we access the GPS data downloaded in 01_access_data.R
full_files <- list.files("/home/hbronnvik/Documents/WS_data/GPS_data", full.names = T)
# the date the data were retrieved from Movebank
dwld_date <- "2024-10-26"

# remove errors
start_time <- Sys.time()
clean_locations <- lapply(1:length(full_files), function(x){
  # load the data
  ind <- readRDS(full_files[x]) # "C:/Users/hbronnvik/Documents/storkSSFs/full_data/78031713_2023-01-28.rds"
  # clean the data
  locs_df <- ind %>% 
    drop_na(location.long) %>% 
    mutate(index = row_number())
  # remove duplicated locations because they prevent accurate calculations of distance and speed
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),] %>% 
    filter(is.na(height.above.ellipsoid))
  
  locs_df <- locs_df %>% 
    filter(!index %in% doubles$index) 
  
  # warn if a duplicated timestamp contains information other than location (not usually the case)
  if(nrow(locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp),
                                                       "timestamp"]),]) > 0){print("Duplicates containing HAE and DOP values exist.", quote = F)}
  
  # remake doubles in the event of duplicates that hold values
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),]%>% 
    mutate(event = round_date(timestamp, "minute"))
  
  if(nrow(doubles) > 0){
    # if there is more than one instance of duplicates with information
    check <- lapply(unique(round_date(doubles$timestamp, "minute")), function(q){
      # take the duplicates within one minute
      dd <- doubles %>% 
        filter(event == q)
      # determine the last location before the duplicates (point of interest)
      poi <- locs_df %>% 
        filter(timestamp < dd$timestamp[1]) %>% 
        slice(n())
      # find the distance from each duplicate to the poi, even for true points this may increase
      cc <- lapply(1:nrow(dd), function(p){
        d <- distVincentyEllipsoid(c(poi$location.long, poi$location.lat), c(dd$location.long[p], dd$location.lat[p]))
        return(d)
      }) %>% unlist()
      # calculate the ground speeds from each duplicate to the poi
      dd <- dd %>% 
        mutate(dist_from_unique = cc, 
               time_since_unique = as.numeric(difftime(dd$timestamp[1], poi$timestamp, units = "sec")),
               speed_after_unique = dist_from_unique/time_since_unique) %>% 
        filter(speed_after_unique > s_thresh)
      return(dd)
    }) %>% reduce(rbind)
    # filter out the ground speeds higher than reasonable
    locs_df <- locs_df %>% 
      filter(!index %in% check$index)
  }
  
  if(nrow(locs_df) > 1){
    locs_df <- locs_df  %>% 
      mutate(date = date(timestamp),
             td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
      # remove bursts
      filter(td >= 300) %>%
      # down-sample to 15 minutes
      mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>%
      group_by(seq15) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(-seq15, -td) %>% 
      # calculate ground speeds
      mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
             timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
             ground_speed_15 = distance/timediff) %>% 
      # remove absurd speeds to clean outliers
      filter(ground_speed_15 < 50)
    
    # add daily metrics between days 
    # (if there is only one point in a day, this fails, but we do not want those birds anyway)
    ind <- locs_df %>% 
      filter(!is.na(location.lat)) %>%
      group_by(date) %>% 
      slice(1, n()) %>% 
      mutate(daily_dist = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat)))) %>% 
      st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
      mutate(dist = st_distance(geometry),
             # the bearing between start and end of a day
             daily_direction = c(NA, lwgeom::st_geod_azimuth(geometry)),
             daily_direction = units::set_units(daily_direction, "radians"),
             # put that in degrees -180 to 180
             daily_direction = as.numeric(units::set_units(daily_direction, "degrees")),
             # define anything below the midline as broadly "southward"
             compass_direction = ifelse(abs(daily_direction) >= 90, "southward", "northward")) %>% 
      slice(2) %>% 
      ungroup() %>% 
      st_drop_geometry() %>% 
      dplyr::select(date, daily_dist, daily_direction, compass_direction)
    
    # ind$daily_direction <- NA
    # for (i in 2:nrow(ind)) {
    #   # print(i)
    #   ind$daily_direction[i] <- bearingRhumb(c(ind$location.long[i-1], ind$location.lat[i-1]),
    #                                          c(ind$location.long[i], ind$location.lat[i]))
    # }
    # ind <- ind %>% 
    #   # the Haversine distance between first locations of consecutive dates when the bird transmitted
    #   mutate(daily_dist = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
    #          # the rhumbline bearing between first locations of consecutive dates when the bird transmitted
    #          # daily_direction = bearingRhumb(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
    #          compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward")
    #   ) %>% 
    #   ungroup() %>% 
    #   arrange(timestamp) %>% 
    #   dplyr::select(date, daily_dist, daily_direction, compass_direction)
    # sp <- ind
    # coordinates(sp) <- ~location.long+location.lat
    # proj4string(sp) <- "EPSG:4326"
    # mapView(sp, zcol = "compass_direction")
    
    # add the daily metrics to the full days
    locs_df <- locs_df %>% 
      left_join(ind, by = join_by(date))
    
    # locs_sf <- locs_df %>%
    #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)
    # mapview::mapview(locs_sf, zcol = "compass_direction")
    
    print(paste0("Cleaned data for ", locs_df$individual.id[1], ", ", x, " of ", length(full_files), "."), quote = F)
    
    return(locs_df)
  }
}) # Time difference of of 44.80227 mins
Sys.time() - start_time

# saveRDS(clean_locations, "/home/hbronnvik/Documents/chapter2/clean_data/clean_locations_2024-11-06.rds")
clean_locations <- readRDS("/home/hbronnvik/Documents/chapter2/clean_data/clean_locations_2024-11-06.rds")

# to check on 2024 downloads:
# ids <- c(4081089905, 4081050803, 4081110497, 4081168215, 4081121150, 4081043020, 4081147270, 3967948529,
#          3968085304, 3968052282, 4000536122, 4000576074, 4000582268, 4000568914, 4060656577, 4060574425,
#          4060647827, 4060597278, 4060632373, 3938980110, 3913673266, 4042698143, 4042709187, 4042735070,
#          4042767872, 3939011728, 3939017397, 3967979972, 3968035956, 3968109777, 3968066416, 4000558649,
#          4000545882, 3913664321, 3913696334, 3913704596, 3913714927, 4042813691)

# Here, we use a double ground speed filter to define migration and extract on which days
# each individual migrated. First, we estimate daily displacements, then we define those as 
# spring or fall based on date and travel angle, then we include slower movements if they are
# within a week of the faster movements.

# n <- 91 # to test a single ID
start_time <- Sys.time()
ld_locs <- lapply(1:length(clean_locations), function(n){tryCatch({
  # one individual cleaned of outliers and errors
  x <- clean_locations[[n]]
  # to check on 2024 downloads:
  # if(unique(x$individual.id) %in% ids){
  #   print("Bird from 2024.", quotes = F)
  # }
  print(paste0("Segmenting tracks for ", x$individual.local.identifier[1], ", individual ", 
               n, " of ", length(clean_locations), "."))
  if(max(x$daily_dist, na.rm = T) > d_thresh){
    # reduce the file
    df <- x %>% 
      # choose columns
      dplyr::select(individual.id, timestamp, location.lat, location.long, daily_dist, daily_direction, ground_speed_15) %>% 
      # choose long distance travel days
      filter(daily_dist > d_thresh) %>%
      arrange(timestamp) %>%
      mutate(timeLag = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
             # insert an enormous time lag for the first location to remove the NA
             timeLag = ifelse(is.na(timeLag), 1e5, timeLag),
             newBurst = ifelse(round(timeLag) <= weeks(1), F, T),
             newCluster = ifelse(round(timeLag) <= weeks(6), F, T),
             # take the cumulative sum to act as a unique ID for each burst
             cumu_check_for_burst = cumsum(newBurst),
             cumu_check_for_clust = cumsum(newCluster)) %>% 
      group_by(cumu_check_for_burst) %>% 
      # for the bursts, calculate the time difference between the last and first locations
      mutate(burstLength = length(unique(date(timestamp))),
             # add an ID to each burst
             burstID = as.character(cur_group_id())) %>%
      ungroup() %>% 
      group_by(cumu_check_for_clust) %>% 
      # for the bursts, calculate the time difference between the last and first locations
      mutate(clusterLength = length(unique(date(timestamp))),
             # add an ID to each burst
             clusterID = as.character(cur_group_id())) %>% 
      ungroup() %>% 
      # clean up the sorting columns
      dplyr::select(-"newBurst", -"cumu_check_for_burst", -"newCluster", -"cumu_check_for_clust")
    
    # sp <- x
    # sp$check <- ifelse(sp$timestamp %in% df$timestamp, "yes", "no")
    # coordinates(sp) <- ~location.long+location.lat
    # proj4string(sp) <- "EPSG:4326"
    # mapView(sp, zcol = "check")
    
    class_df <- df %>% 
      group_by(burstID) %>% 
      mutate(burst_angle = ifelse(location.lat[1] > location.lat[n()], "south", "north"),
             burst_season = ifelse(month(timestamp[n()]) %in% c(8:11), "late", 
                                   ifelse(month(timestamp[n()]) %in% c(1:6), "early", "change")),
             burst_class = ifelse(burst_season == "early", "spring", 
                                  ifelse(burst_season == "late", "fall",
                                         ifelse(burst_angle == "north" & burst_season == "change", "spring", 
                                                ifelse(burst_angle == "south" & burst_season == "change", "fall", NA)))),
             burst_year = ifelse(month(timestamp[n()]) == 12 & burst_class == "spring", year(years(1) + timestamp[n()]), year(timestamp[n()])),
             trackID = paste(unique(individual.id), burst_class, burst_year, sep = "_")) %>% 
      ungroup()
    
    # sp <- class_df
    # coordinates(sp) <- ~location.long+location.lat
    # proj4string(sp) <- "EPSG:4326"
    # mapView(sp, zcol = "burstID")
    
    # remove ragged edges of tracks the bird survived
    # if the bird died on a ragged edge, keep those data, the bird might have continued if it had survived
    
    # identify ragged edges as bursts that do not have high speed days in them
    ragged <- class_df %>% 
      group_by(burstID) %>% 
      mutate(high_speed = max(daily_dist > 70000)) %>%
      ungroup() %>% 
      filter(high_speed == F)
    # did the bird die?
    dod <- deaths %>% 
      filter(individual.id == unique(x$individual.id)) %>% 
      dplyr::select(dod) %>% 
      deframe()
    # if there is no date of death
    if(length(dod) == 0){
      # if the date of death is missing, use the last transmission (someone might have ended the deployment without a death),
      last_date <- ref_info$tag_timestamp_end[which(ref_info$individual.id == unique(x$individual.id))]
      if(last_date > as.Date(dwld_date)-days(7)){
        # if the bird survived until the download of the data (or its tag did),
        dod <- NA
      }else{
        dod <- last_date
      }
    }
    if(nrow(ragged) > 0){
      if(!is.na(dod)){
        # did the bird die in a ragged edge?
        cutoff <- date(dod) %in% unique(date(ragged$timestamp))
        if(cutoff == T){
          # how ragged?
          messy <- lapply(2:length(class_df$burstID), function(z){
            temp <- class_df %>%
              filter(burstID == z)
            temp_prev <- class_df %>%
              filter(burstID == z-1)
            temp$prev_long <- !unique(temp_prev$burstID) %in% unique(ragged$burstID)
            temp
          }) %>% reduce(rbind)
          # these are preceded by high speed bursts, but are low speed
          messy <- messy %>% 
            filter(prev_long == T & burstID %in% unique(ragged$burstID))
          # did the bird die in a ragged edge adjacent to a high speed burst?
          keep <- date(dod) %in% unique(date(messy$timestamp))
          if(keep == T){
            edge <- unique(class_df$burstID[date(class_df$timestamp) == date(dod)])
            # discard ragged edges that do not contain a death
            class_df <- class_df %>% 
              filter(!burstID %in% ragged$burstID | burstID == edge)
          }
        }else{
          class_df <- class_df %>% 
            filter(!burstID %in% ragged$burstID)
        }
      }else{
        class_df <- class_df %>% 
          filter(!burstID %in% ragged$burstID)
      }
    }
    
    # sp <- class_df
    # coordinates(sp) <- ~location.long+location.lat
    # proj4string(sp) <- "EPSG:4326"
    # mapView(sp, zcol = "trackID")
    
    if(nrow(class_df) > 0){
      tracks <- class_df %>%
        group_by(trackID) %>% 
        slice(1, n())
      tracks <- lapply(split(tracks, tracks$trackID), function(c){
        ID <- unique(c$trackID)
        track <- x %>% 
          filter(between(timestamp, c$timestamp[1], c$timestamp[2])) %>% 
          mutate(trackID = ID,
                 track_displacement = distVincentyEllipsoid(c(head(location.long, 1), head(location.lat, 1)),
                                                            c(tail(location.long, 1), tail(location.lat, 1)))/1000)
        track
      }) %>% rbindlist()
      
      # ggplot(tracks, aes(location.long, location.lat, color = trackID)) +
      #   borders("world", xlim = c(-10, 10), ylim = c(0, 60), colour = "black") +
      #   geom_point() +
      #   theme_classic() +
      #   theme(legend.position = "none",
      #         text = element_text(size = 20, color = "black"),
      #         axis.line = element_line(color = "black"),
      #         axis.text = element_text(size = 20, color = "black")) +
      #   facet_wrap(~trackID)
      # 
      # x <- x %>% 
      #   rowwise() %>% 
      #   mutate(trackID = ifelse(timestamp %in% tracks$timestamp, tracks$trackID[which(tracks$timestamp == timestamp)], NA)) %>% 
      #   ungroup()
      
      suppressMessages(x <- x %>% 
                         left_join(tracks))
      
      # sp <- x# %>% filter(trackID == "1178289602_spring_2021")
      # coordinates(sp) <- ~location.long+location.lat
      # proj4string(sp) <- "EPSG:4326"
      # mapView(sp, zcol = "trackID")
      
      return(x)
    }else{
      x$trackID <- NA
      return(x)
    }
  }else{
    x$trackID <- NA
    return(x)
  }
}, error = function(e){
  print(geterrmessage())
  if(geterrmessage() == "Must request at least one colour from a hue palette."){dev.off()}
})
})
Sys.time()-start_time # Time difference of 49.60893 secs

# reduce memory burden
rm(clean_locations)
gc()

# compress the list
migration_locations <- ld_locs[lapply(ld_locs, length) > 1]
migration_locations <- data.table::rbindlist(migration_locations, fill = T)
migration_locations <- migration_locations %>% 
  filter(!is.na(trackID))

# get the track IDs from the migratory birds
rs_ids <- migration_locations %>% 
  group_by(individual.id, trackID) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, trackID) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("fall", trackID), "fall", "spring")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# saveRDS(meta, file = "/home/hbronnvik/Documents/chapter2/metadata_tracks_2024-11-06.rds")

# add the number of journeys
migration_locations <- migration_locations %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$trackID == trackID)]) %>% 
  ungroup()
# find any animals that slipped through the cracks for visual estimates of DOD
lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id) %>% 
  filter(!individual.id %in% deaths$individual.id)

# # define birds that took eastern routes as ones that are ever east of 16.5 longitude (East Germany)
# eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 16.5])
# # remove the eastern birds
# migration_locations <- migration_locations %>% 
#   filter(!individual.id %in% eastern_birds) %>% 
#   group_by(trackID) %>%  
#   ungroup() 

migration_locations %>% 
  drop_na(trackID) %>% 
  group_by(trackID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(track_displacement) %>% 
  deframe() %>% 
  hist(breaks = 100)

# filter the short routes out:
migration_locations <- migration_locations %>%
  filter(track_displacement > 40) # only tracks longer than 40 km

# find all the final transmission dates
info <- lapply(studies, function(x){
  df <- getMovebankAnimals(x, loginStored) %>% 
    filter(number_of_events > 0) %>% 
    mutate(timestamp_end = as.POSIXct(sub("\\.000", "", timestamp_end), tz = "UTC"),
           study_id = x) %>% 
    filter(sensor_type_id == 653)
}) %>% 
  reduce(rbind) %>% 
  filter(individual_id %in% migration_locations$individual.id)

# find the success/failure of each migration track
migration_locations_ls <- migration_locations %>% 
  group_by(trackID) %>% 
  group_split()
migration_conditions <- lapply(migration_locations_ls, function(x){
  # take the migratory route
  x <- x %>% 
    arrange(timestamp) 
  
  # if the last GPS time point is greater than the time the data were downloaded, use the download date
  t_end <- ifelse(info$timestamp_end[info$individual_id == unique(x$individual.id)] > dwld_date, 
                  dwld_date, info$timestamp_end[info$individual_id == unique(x$individual.id)]) %>% 
    as.POSIXct(tz = "UTC") %>% 
    date()
  # take either the confirmed DOD or the last transmitted time stamp
  if(unique(x$individual.id) %in% deaths$individual.id){
    loss_time <- deaths %>% 
      filter(individual.id == unique(x$individual.id)) %>% 
      mutate(dod = as.Date(ifelse(is.na(dod), as.character(t_end), as.character(dod)))) %>% 
      dplyr::select(dod) %>% 
      deframe()
  }else{
    loss_time <- t_end
  }
  
  # compare the DOD to the end of the migration
  loss <- max(x$timestamp) > loss_time - days(s_thresh)
  # add a column containing the outcome of the migratory track
  # x <- x %>% 
  #   mutate(track_status = ifelse(loss == T & trackID == unique(x$trackID), "incomplete", "complete"))
  status <- data.frame(trackID = unique(x$trackID),
                       track_status = ifelse(loss == T, "incomplete", "complete"))
  # gc()
  return(status)
}) %>% reduce(rbind)

rm(migration_locations_ls);gc()

migration_locations <- migration_locations %>% 
  left_join(migration_conditions)

# complete_ml <- migration_locations %>% 
#   filter(track_status == "complete")

# save out the data of the tracks longer than 40 km labeled with their outcomes
# saveRDS(migration_locations, file = paste0("/home/hbronnvik/Documents/chapter2/migration_locations_40km70km30daySpeed_", Sys.Date(),".rds"))

migration_locations <- readRDS("/home/hbronnvik/Documents/chapter2/migration_locations_40km70km30daySpeed_2024-11-06.rds")

# facts about the migrations attempted and completed by each animal in each season 
metad <- migration_locations %>%
  group_by(trackID) %>% 
  filter(daily_dist >= 40000) %>% 
  mutate(start_migration = timestamp[1],
         end_migration = timestamp[n()],
         id_date = paste0(individual.id, "_", date)) %>% 
  group_by(id_date) %>% 
  slice(1) %>%
  ungroup() %>% 
  dplyr::select(individual.id, trackID, id_date, track_displacement, track_status, daily_dist, start_migration, end_migration)

# save out the data of the tracks longer than 40 km labeled with their outcomes
# saveRDS(metad, file = paste0("/home/hbronnvik/Documents/chapter2/migration_dates_", Sys.Date(),".rds"))
