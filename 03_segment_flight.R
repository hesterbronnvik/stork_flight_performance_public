### Identify birds to include, segment their migration tracks, save dates out
### Written by Dr. Elham Nourani for application to Honey Buzzards
### https://github.com/mahle68/HB_ontogeny/blob/main/02_GPS_segmentation.R
### Adapted by Hester Brønnvik
### 12.10.2023
### hbronnvik@ab.mpg.de

library(sf)
library(lwgeom)
library(move)
library(data.table)
library(tidyverse)
wgs <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")

# all the GPS data downloaded in 01_access_data.R
full_data_files <- list.files("/home/hbronnvik/Documents/WS_data/GPS_data", full.names = T)

# the migration dates generated in 02_segment_tracks.R
migrants <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds") %>% 
  arrange(individual.id) %>% 
  dplyr::select(individual.id) %>% 
  deframe() %>% 
  unique()

files <- sapply(migrants, function(x){
  file <- grep(x, full_data_files, value = T)
  file
})

# the files containing migrant data
files <- files[lapply(files, length) > 0]

# the days of the migrations
migration_dates <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-11-06.rds")
# the columns worth keeping
colsToKeep <- c("individual.id","timestamp", "location.long", "location.lat", "height.above.ellipsoid",
                "gps.satellite.count","event.id", "gps.dop","eobs.status",
                "ground.speed","eobs.horizontal.accuracy.estimate",
                "study.name", "tag.local.identifier", "individual.local.identifier","individual.taxon.canonical.name")

# clean the data of GPS error and duplicate locations
# then turn into a move object with vertical speeds, ground speeds, and turning angles
lapply(1:length(files), function(n){
  ind <- readRDS(files[n])
  print(n)
  # Check for missing info in time and coords
  # These usually correspond to eobs status B to D: 
  #A = position and time within accuracy masks
  #B = only time of week and weeknumber valid
  #C = only weeknumber valid
  #D = no valid data
  
  # anyNA(df$timestamp)
  # anyNA(df$individual.local.identifier)
  # anyNA(df$location.long); anyNA(df$location.lat)
  # table(df$eobs.status)
  
  ind <- ind[ind$eobs.status %in% c("","A"),]
  ind <- ind[!is.na(ind$location.long),]
  
  if(nrow(ind) > 0){
    # remove unnecessary columns
    ind <- data.frame(ind)[,colsToKeep]
    ind$event.id <- as.character(ind$event.id) #important for solving the error: 'names' attribute [13] must be the same length as the vector [12]
    
    animalName <- unique(ind$individual.id)
    # Check for duplicates
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
          filter(speed_after_unique > 50)
        return(dd)
      }) %>% reduce(rbind)
      # filter out the ground speeds higher than reasonable
      locs_df <- locs_df %>% 
        filter(!index %in% check$index)
    }
    
    # filter out winters and breeding
    ind_dates <- migration_dates %>% 
      filter(individual.id == unique(locs_df$individual.id))%>% 
      rowwise() %>% 
      mutate(date = as.Date(str_split(id_date, "_")[[1]][2]))
    
    locs_df <- locs_df %>% 
      filter(date(timestamp) < min(ind_dates$date) | date(timestamp) %in% ind_dates$date)
    
    # calculate track variables
    mv <- move(x=locs_df$location.long, y=locs_df$location.lat, 
               time=locs_df$timestamp,
               proj=crs("+proj=longlat +ellps=WGS84"),
               animal=locs_df$individual.local.identifier,
               data=locs_df)
    mv$timelag.sec <- c(NA,timeLag(mv, units="secs"))
    mv$altitude.diff <- c(NA,(mv$height.above.ellipsoid[-1] - mv$height.above.ellipsoid[-nrow(mv)]))
    mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
    mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
    mv$step.length <- c(NA,move::distance(mv))
    mv$gr.speed <- c(NA, speed(mv))
    # save on the hard drive
    saveRDS(mv, file = paste0("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs/",animalName,"_gpsNoDup_moveObj_1024.rds"))
  }
  
})

#define variables for segmentation
min_res_tl <- 2 # 1 to max 2 sec time lag
min_burst_d <- 30 # we want bursts of at least 30 secs
sw_vs <- 2 #smoothing window of 5 seconds (> min burst duration, 2 before 2 after each loc) for vertical speed for later classification
sw_th <- 12 # HALF of the smoothing window of interest. 25 seconds for the entire thermaling behavior. 1 is row i, sw_th rows before and sw_th rows after
circl_deg <- 250 #degrees of rotation to be achieved in the time defined by sw_th*2+1 for full circular soaring
circl_deg_sh <- 145 #degrees of rotation to be achieved in the time defined by sw_th*2+1 for shallow thermal soaring
min_behav_d <- 7 #minimum duration in seconds of a specific behavior, when less than this and if in between two segments of a different behavior it will be incorporated in the previous and following segment
min_therm_d <- 20 #minimum duration for a circling event to be considered as thermalling

#calculate rotation per second required for each circling mode
ta_per_sec_circl <- circl_deg/((2*sw_th)+1)
ta_per_sec_shallow <- circl_deg_sh/((2*sw_th)+1)

fls <- list.files("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs", full.names = T)

# find 1 Hz data files
check <- lapply(list.files("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs", full.names = T), function(f){
  build <- readRDS(f)
  data.frame(file = f,
             hz = min(build$timelag.sec, na.rm = T))
}) %>% reduce(rbind)

# to remove data that have been processed:
# done <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s_2.0", full.names = T)
# 
# done_birds <- data.frame(done = done) %>% 
#   mutate(bird = sub("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s_2.0/", "",
#                     sub("_segmented_bursts.rds", "", done)))
# 
# todo <- list.files("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs", full.names = T)
# todo_birds <- data.frame(todo = todo) %>% 
#   mutate(todo_bird = sub("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs/", "",
#                          sub("_gpsNoDup_moveObj_1024.rds", "", todo))) %>% 
#   filter(!todo_bird %in% done_birds$bird)
# 
# fls <- check %>% 
#   mutate(check_bird = sub("/home/hbronnvik/Documents/WS_data/clean_geom_mv_objs/", "",
#                           sub("_gpsNoDup_moveObj_1024.rds", "", file))) %>% 
#   # filter(between(hz, .01, 1.2)) %>% 
#   filter(check_bird %in% done_birds$bird) %>% 
#   dplyr::select(file) %>% 
#   deframe()
# 
# disagree_days <- metad %>% 
#   filter(!id_date %in% old_metad$id_date)

(b <- Sys.time())
lapply(fls[1:length(fls)], function(f){
  print(paste0("Processing file ", which(fls == f), " of ", length(fls), "."))
  
  
  # STEP 1: calculate track variables -------------------------------------------------------------------------
  ind <- readRDS(f) 
  
  animalID <- ind@idData$individual.id
  ind <- ind %>% st_as_sf(coords = c("location_long", "location_lat"), crs = wgs)
  
  # select days that were missing under previous track segmentation
  ind <- ind %>% 
    mutate(id_date = paste(animalID, date(timestamp), sep = "_")) %>% 
    filter(id_date %in% disagree_days$id_date)
  
  if(nrow(ind) > 0){
    #everything is calculated from one point compared to its previous
    ind <- ind %>% 
      group_by(timestamp) %>% #remove duplicated timestamps
      slice(1) %>% 
      ungroup() %>% 
      arrange(timestamp) %>% 
      mutate(time_lag_sec = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "secs") %>% as.numeric()),
             burst_id = cumsum(time_lag_sec > min_res_tl), #increase burst id by one, every time time_lag is more than min_res_tl
             azimuth_rad = c(NA, st_geod_azimuth(.)),
             step_length = if_else(row_number() == 1, 0, st_distance(geometry, lag(geometry), by_element = TRUE) %>% as.numeric()), #meters
             azimuth =  (azimuth_rad * 180)/pi,
             azimuth_positive = if_else(azimuth >= 0, azimuth, azimuth + 360)) %>% #assign one burst id to rows of consecutive 1Hz gps
      group_by(burst_id) %>% 
      mutate(altitude_diff = if_else(row_number() == 1, NA, height.above.ellipsoid - lag(height.above.ellipsoid,1)),
             vert_speed = if_else(row_number() == 1, NA, altitude_diff/time_lag_sec), #m/s
             turn_angle = if_else(row_number() == 1 | row_number() == n(), NA,
                                  ((lag(azimuth_positive) - azimuth_positive + 180) %% 360) - 180),
             ground_speed = if_else(row_number() == 1, NA, step_length/time_lag_sec)) %>% 
      ungroup() %>% 
      # only times when the animal is fast enough for flight
      filter(ground_speed >= 3.5)
    
    #garbage collection to free up RAM
    gc()
    
    #calculate nrow of each burst
    bursts_to_keep <- ind %>% 
      st_drop_geometry() %>% 
      group_by(burst_id) %>% 
      summarize(freq = n()) %>% 
      filter(freq >= min_burst_d) %>% 
      pull(burst_id)
    
    #keep bursts with more rows than min_burst_d 
    ind <- ind %>% 
      st_drop_geometry()  %>% 
      filter(burst_id %in% bursts_to_keep)
    
    gc()
    
    if(nrow(ind) > 0){
      
      #split each individual by burst_id
      burst_ls_corr <- ind %>% 
        group_split(burst_id)  
      
      # Compute smoothed turning angle separately for each burst
      # apply a moving window of sw_th*2+1 s to calculate the absolute cumulative sum of the turning angles (hereafter cumulative turning angle) and 
      # a moving window of sw_vs*2+1 s to calculate the average vertical speed
      
      df_hr <- lapply(burst_ls_corr, function(b){ 
        b <- b %>%
          mutate(vert_speed_smooth = ifelse(row_number() <= sw_vs | row_number() > n() - sw_vs,
                                            NA,
                                            map_dbl((sw_vs + 1):(n() - sw_vs), ~ { #Within the ifelse() call, we use map_dbl() from the purrr package to loop through the row indices for which smoothing can be performed ((sw_vs + 1):(n() - sw_vs))
                                              #Inside the map_dbl() call, we define an anonymous function that calculates the mean of vert_speed values within the range of (i - sw_vs) to (i + sw_vs) for each valid row index i. The na.rm = TRUE argument ensures that NA values are removed before calculating the mean.
                                              i <- .x
                                              mean(b %>% slice((i - sw_vs):(i + sw_vs)) %>% pull(vert_speed), na.rm = T)
                                            })),
                 turn_angle_cum = ifelse(row_number() <= sw_th | row_number() > n() - sw_th,
                                         NA,
                                         map_dbl((sw_th + 1):(n() - sw_th), ~ { 
                                           i <- .x
                                           max(abs(cumsum(b %>% slice((i - sw_th):(i + sw_th)) %>% drop_na(turn_angle) %>% pull(turn_angle))))
                                         })))
        return(b)
      }) %>% 
        bind_rows()
      # classify soaring only based on vertical speed using k-means clustering
      df_hr <- df_hr %>% 
        drop_na(vert_speed_smooth) %>% 
        mutate(kmean_cluster = kmeans(vert_speed_smooth, 2)$cluster) #get two clusters
      
      #find the cluster id that matches soaring (i.e. the cluster with the higher mean vertical speed)
      soar_id <- df_hr %>%
        group_by(kmean_cluster) %>% #take the mean of smoothed vertical speed for each kmeans cluster
        summarise(mean_vert_speed_smooth = mean(vert_speed_smooth)) %>%
        ungroup() %>%
        filter(mean_vert_speed_smooth == max(mean_vert_speed_smooth)) %>%
        pull(kmean_cluster) 
      
      # classify flight based on cumulative rotation at each point AND the clustering based on vertical speed
      df_hr <- df_hr %>% 
        mutate(soar_clust = if_else(kmean_cluster == soar_id, "soar", "glide"),
               flight_clust = if_else(is.na(turn_angle_cum), NA,
                                      if_else(soar_clust == "glide", "gliding",
                                              if_else(turn_angle_cum >= circl_deg, "circular_soaring",
                                                      if_else(turn_angle_cum >= circl_deg_sh, "shallow_circular_soaring",
                                                              if_else(turn_angle_cum < circl_deg_sh, "linear_soaring",
                                                                      "other")))))
        ) 
      
      #visual check
      #mapview::mapview(df_hr, zcol = "flight_clust")
      
      # Add some steps of smoothing based on duration of behaviors:
      burst_ls_class <- df_hr %>% 
        group_split(burst_id)
      
      bursts_smooth <- lapply(burst_ls_class, function(b){ #sample:  b <- burst_ls_class[[1]]
        
        ##--------------------- SMOOTHING level 0: assign classification of the time window at the start and end of the segment (Martina does this as the third smoothing stage)
        # these first and last points can only be classified as either gliding or linear, as their classification was only based on vertical speed but not turning angle
        #calculate sum of turn_angles. If they are classified as circular soaring, but are straighter than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign linear soaring
        
        b <- b %>% 
          mutate(na_index = if_else(!is.na(flight_clust), NA, 
                                    if_else(row_number() %in% 1:(sw_th-2), "head", "tail"))) %>% #give a separate ID to the rows of NAs at the beginning and end of the segment
          group_by(na_index) %>% 
          mutate(duration = sum(time_lag_sec),
                 turn_angle_sum = abs(sum(turn_angle, na.rm = T))) %>% 
          mutate(flight_clust = if_else(!is.na(flight_clust), flight_clust, #if the row already had a flight cluster assigned to it, keep the assignment
                                        if_else(soar_clust == "glide", "gliding",
                                                if_else(soar_clust == "soar" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring",
                                                        if_else(soar_clust == "soar" & between((turn_angle_sum/duration), ta_per_sec_shallow, ta_per_sec_circl), "shallow_circular_soaring", 
                                                                if_else(soar_clust == "soar" & (turn_angle_sum/duration) < ta_per_sec_shallow, "linear_soaring", "linear_soaring")))))) %>% 
          ungroup() %>% 
          dplyr::select(-c(duration, turn_angle_sum, na_index))
        
        ##--------------------- SMOOTHING level 1: based on duration
        #Point-based smoothing and assigning flight segment ID
        b <- b %>% 
          arrange(timestamp) %>% 
          #Before getting into flight segments, do a point-wise smoothing.
          #specifically, deal with the situations where for example there are two points within a soaring bout with 2 different assignments
          #if the previous point and the point after next are circular soaring, assign circular soaring
          mutate(flight_clust_sm = if_else(row_number() %in% c(1, 2, nrow(.) - 1, nrow(.)), flight_clust, #keep the first and last points' assignments the same. otherwise they become NA
                                           if_else(lag(flight_clust, 1) == "circular_soaring" & lead(flight_clust, 2) == "circular_soaring", "circular_soaring",
                                                   #similarly, if the previous point and the point after next are shallow circular soaring, assign shallow circular soaring
                                                   if_else(lag(flight_clust, 1) == "shallow_circular_soaring" & lead(flight_clust, 2) == "shallow_circular_soaring", "shallow_circular_soaring",
                                                           #if the point before previous point and the point after are circular soaring, assign circular soaring
                                                           if_else(lag(flight_clust, 2) == "circular_soaring" & lead(flight_clust, 1) == "circular_soaring", "circular_soaring",
                                                                   #similarly, if the point before previous point and the point after are shallow circular soaring, assign shallow circular soaring
                                                                   if_else(lag(flight_clust, 2) == "shallow_circular_soaring" & lead(flight_clust, 1) == "shallow_circular_soaring", "shallow_circular_soaring",
                                                                           flight_clust)))))) %>% 
          #assign a unique flight segment id to consecutive rows with the same flight category
          mutate(flight_clust_diff = if_else(c(0, diff(as.numeric(as.factor(flight_clust_sm)))) == 0, 0, 1), 
                 flight_seg_id = cumsum(coalesce(flight_clust_diff,0))) #coalesce ignores the NAs
        
        #Segment-based smoothing. First calculate the duration and unique class of each behavioral segment   
        behav_duration <- b %>% 
          group_by(flight_seg_id) %>% 
          reframe(timestamp_start = head(timestamp,1),
                  duration = sum(time_lag_sec),
                  flight_clust_sm = unique(flight_clust_sm)) %>% 
          ungroup() %>% 
          arrange(timestamp_start) %>% 
          #First, if the assignment is linear soaring, duration is <5 sec, and the behav before and after are shallow soaring, assign shallow soaring
          mutate(flight_clust_sm1 = if_else(flight_clust_sm == "linear_soaring" & duration <= 17 & #the very shallow circles are 19 seconds long. so if only the start and end are classified as shallow soaring, the whole bout should be classified as such
                                              lag(flight_clust_sm, 1) == "shallow_circular_soaring" &
                                              lead(flight_clust_sm, 1) == "shallow_circular_soaring", "shallow_circular_soaring", flight_clust_sm)) %>% 
          #Next, more generally, unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
          #keep the assignments for the first and last two points the same. otherwise they will get NAs based on the next rule
          mutate(flight_clust_sm1 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm, 
                                            #if the previous and next flight clust are the same and the current segment is less than min_behav_d sec, assign the prev/next flight clust
                                            if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1))) %>% 
          dplyr::select(flight_seg_id, flight_clust_sm, flight_clust_sm1)
        
        #merge the new clusters with the original data
        b <- b %>% 
          left_join(behav_duration, by = join_by("flight_seg_id","flight_clust_sm")) %>% 
          #calculate flight_seg_id for the new clustering
          mutate(flight_clust_diff_sm1 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm1)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
                 flight_seg_id_sm1 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm1), NA,
                                             cumsum(coalesce(flight_clust_diff_sm1, 0))))
        #visual sanity check
        #mapview::mapview(b, zcol = "flight_clust_sm1")
        
        ##--------------------- SMOOTHING level 2: REPEAT smoothing based on segment duration
        #Segment-based smoothing. First calculate the duration and unique class of each behavioral segment   
        behav_duration2 <- b %>% 
          group_by(flight_seg_id_sm1) %>% 
          reframe(timestamp_start = head(timestamp,1),
                  duration = sum(time_lag_sec),
                  flight_clust_sm1 = unique(flight_clust_sm1)) %>% 
          ungroup() %>% 
          arrange(timestamp_start) %>% 
          #unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
          #keep the assignments for the first and last two points the same. otherwise they will get NAs based on the next rule
          mutate(flight_clust_sm2 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm1, 
                                            #if the previous and next flight clust are the same and the current segment is less than min_behav_d sec, assign the prev/next flight clust
                                            if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1)))
        
        #merge the new clusters with the original data
        b <- b %>% 
          left_join(behav_duration2 %>% dplyr::select(flight_seg_id_sm1, flight_clust_sm1, flight_clust_sm2), by = c("flight_seg_id_sm1","flight_clust_sm1")) %>% 
          #calculate flight_seg_id for the new clustering
          mutate(flight_clust_diff_sm2 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm2)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
                 flight_seg_id_sm2 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm2), NA,
                                             cumsum(coalesce(flight_clust_diff_sm2, 0))))
        # mapview::mapview(b, zcol = "flight_clust_sm2")
        
        ##--------------------- SMOOTHING level 3: based on turning angle of circular soaring bouts (Martina's code did this based on circling duration)
        #calculate sum of turn_angles. If they are classified as circular soaring, but are straighter than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign linear soaring
        #alternatively, if the assignment is linear soaring, but the rate of rotation is more than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign the appropriate soaring class
        
        behav_ta <- b %>% 
          group_by(flight_seg_id_sm2, flight_clust_sm2) %>% 
          reframe(duration = sum(time_lag_sec),
                  turn_angle_sum = sum(abs(turn_angle), na.rm = T)) %>% 
          mutate(flight_clust_sm3 = if_else(flight_clust_sm2 == "circular_soaring" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring",
                                            if_else(flight_clust_sm2 == "circular_soaring" & (turn_angle_sum/duration) < ta_per_sec_shallow, "linear_soaring", 
                                                    if_else(flight_clust_sm2 == "linear_soaring" & (turn_angle_sum/duration) > ta_per_sec_circl, "circular_soaring",
                                                            if_else(flight_clust_sm2 == "linear_soaring" & (turn_angle_sum/duration) > ta_per_sec_shallow, "shallow_circular_soaring",
                                                                    flight_clust_sm2)))))
        
        #merge the new clusters with the original data
        b <- b %>% 
          left_join(behav_ta %>% dplyr::select(flight_seg_id_sm2, flight_clust_sm2, flight_clust_sm3), by = c("flight_seg_id_sm2","flight_clust_sm2")) %>% 
          #calculate flight_seg_id for the new clustering
          mutate(flight_clust_diff_sm3 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm3)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
                 flight_seg_id_sm3 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm3), NA,
                                             cumsum(coalesce(flight_clust_diff_sm3,0))))
        
        #visual sanity check
        #mapview::mapview(b, zcol = "flight_clust_sm3")
        
        # Assign unique ID to the behavioral segment based on the final smoothest classification
        b <- b %>% 
          mutate(track_flight_seg_id = paste(animalID, "burst", unique(b$burst_id), "seg", b$flight_seg_id_sm3, sep = "_")) 
        
        #return each classified and smoothed burst to a list
        return(b) 
      }) %>%
        #rbind all bursts and save classified and smoothed dataframe per individual
        bind_rows() %>% 
        mutate(individual.id = animalID)
      
      # save locally
      saveRDS(bursts_smooth, file = paste0("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s/", 
                                           animalID,"_segmented_bursts.rds"))
      gc()
    } 
  }
})
Sys.time() - b #2.249033 hours for 13 individuals
# Time difference of 1.56451 days 410 birds


