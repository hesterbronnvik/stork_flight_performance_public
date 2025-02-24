### Measure soaring-gliding efficiency in white storks
### Hester Brønnvik
### 2025-01-08
### hbronnvik@ab.mpg.de

library(tidyverse)
library(geosphere)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 10), 
                  text = element_text(size = 15),
                  strip.background = element_rect(fill = "white", color = "white")))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

# the data with flight classifications
fls <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)
# the metadata with track IDs and dates
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
# remove eastern migrants
{
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
classified <- classified[!is.na(classified$trackID),]}

eastern_IDs <- classified %>% 
  filter(location.long >= 16) %>% 
  dplyr::select(individual.id) %>% 
  deframe() %>% 
  unique()

rm(classified);rm(status)

sg_eff <- lapply(fls, function(x){
  print(x)
  # get an ID
  bird <- readRDS(x) %>% 
    dplyr::select(individual.id, timestamp, location.lat, location.long, height.above.ellipsoid,
                  burst_id, track_flight_seg_id, flight_clust_sm3)
  if(!unique(bird$individual.id) %in% eastern_IDs){
    # reclassify soaring to match what we did in wind estimation (wind_estimation.R ln 382)
    # reclassify any shallow soaring as circular if it is brief and leads into or follows circular soaring
    tween_behav <- bird %>% 
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
    bird <- bird %>% 
      # smooth a 4th time
      mutate(flight_clust_sm3 = ifelse(track_flight_seg_id %in% tween_behav$track_flight_seg_id,
                                       "circular_soaring", flight_clust_sm3)) %>% 
      rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
      mutate(ind_burst_id = paste0(individual.id, "_", burst_id))
    # get migration
    record <- m_days %>% 
      filter(individual.id == unique(bird$individual.id))
    bird <- bird %>% 
      mutate(date = date(timestamp)) %>% 
      filter(date %in% record$date) %>% 
      left_join(record %>% dplyr::select(trackID, track_status, date))
    # for each hour, get soaring times
    hr_bird <- bird %>% 
      mutate(hr = paste0(date(timestamp), "_", hour(timestamp))) %>% 
      group_by(hr, flight_clust_sm3) %>% 
      summarize(track_id = unique(trackID),
                track_status = unique(track_status),
                mode_obs = n())
    # for each hour, get glide distances
    glides <- bird %>% 
      filter(flight_clust_sm3 == "gliding") %>% 
      # each glide
      mutate(gap = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
             glide_id = cumsum(ifelse(is.na(gap) | gap == 1, F, T))) %>% 
      # distance per glide
      group_by(glide_id) %>% 
      # start and end points
      slice(1, n()) %>% 
      # great-circle-distance on the ellipsoid
      mutate(glide_displacement = distVincentyEllipsoid(cbind(location.long[1], location.lat[1]),
                                                        cbind(location.long[2], location.lat[2]))) %>% 
      # one distance per glide (in meters)
      slice(1) %>% 
      ungroup()
    hr_glide <- glides %>% 
      mutate(hr = paste0(date(timestamp), "_", hour(timestamp))) %>% 
      group_by(hr) %>% 
      summarize(track_id = unique(trackID),
                track_status = unique(track_status),
                glide_dist = sum(glide_displacement))
    # estimate hourly efficiency
    efficiency <- hr_bird %>% 
      filter(flight_clust_sm3 == "circular_soaring") %>% 
      left_join(hr_glide) %>% 
      mutate(efficiency = glide_dist/mode_obs)
    
    # map for sanity
    # bird %>%
    #   mutate(hr = paste0(date(timestamp), "_", hour(timestamp))) %>%
    #   filter(hr == "2020-08-12_13") %>%
    #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>%
    #   mapview::mapview(zcol = "flight_clust_sm3")
    return(efficiency)
  }
}) %>% reduce(rbind)

sg_eff <- sg_eff %>% 
  rowwise() %>% 
  mutate(individual.id = str_split(track_id, "_")[[1]][1],
         season = str_split(track_id, "_")[[1]][2]) %>% 
  ungroup() %>% 
  left_join(records %>% dplyr::select(trackID, journey_number) %>% rename(track_id = trackID))

sg_eff %>% 
  filter(journey_number < 5) %>% 
  mutate(age = paste(season, journey_number)) %>% 
  ggplot(aes(log(efficiency), fill = age)) +
  geom_density() +
  scale_fill_manual(values = c("#EE5E53", "#F07268", "#F38979", "#FABBA0", "#0081A7", "#009BB1", "#24B5B8", "#B5CDB7")) +
  labs(x = "log per-hour soaring-gliding efficiency (m/s)", y = "Density", fill = "Age") +
  facet_wrap(~season)

# an extreme example:
# readRDS(fls[grep("4003322813", fls)]) %>%
#   mutate(hr = paste0(date(timestamp), "_", hour(timestamp))) %>%
#   filter(hr == "2024-08-19_17") %>%
#   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>%
#   mapview::mapview(zcol = "flight_clust_sm3")

# from 06_estimate_thermal_loss.R ln 585
hr_losses <- readRDS("/home/hbronnvik/Documents/chapter2/thermal_losses_20241130.rds") %>% 
  group_by(individual.id) %>% 
  mutate(hr = paste0(date(timestamp), "_", hour(timestamp)),
         individual.id = as.character(individual.id)) %>% 
  group_by(individual.id, hr) %>% 
  summarize(losses = sum(ridden))

sg_eff %>% 
  left_join(hr_losses) %>% 
  drop_na(losses) %>%
  mutate(lost = ifelse(losses > 0, "Loss", "No loss"),
         log_eff = log(efficiency)) %>% 
  filter(!is.infinite(log_eff) & !is.na(log_eff)) %>% 
  group_by(season, journey_number, losses) %>% 
  summarize(obs = n())

sg_eff %>% 
  left_join(hr_losses) %>% 
  drop_na(losses) %>% 
  ggplot(aes(as.factor(losses), log(efficiency))) +
  geom_boxplot() +
  facet_wrap(~season+journey_number, nrow = 2)

ridge_pls <- lapply(c("spring", "fall"), function(s){
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
  
  rp <- sg_eff %>% 
    left_join(hr_losses) %>% 
    drop_na(losses) %>%
    mutate(lost = ifelse(losses > 0, "Loss", "No loss"),
           log_eff = log(efficiency)) %>% 
    filter(!is.infinite(log_eff) & !is.na(log_eff)) %>% 
    filter(season == s) %>% 
    ggplot(aes(log_eff, as.character(journey_number), fill = lost)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                  panel_scaling = T, alpha = 0.85) +
    scale_fill_manual(values = cols, name = "") +
    labs(x = "log per-hour soaring-gliding efficiency (m/s)", y = "Age", title = "") +
    theme(legend.position = "top")
  return(rp)
})
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/hourly_efficiency_loss.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(ridge_pls[[2]], ridge_pls[[1]])
# dev.off()

sg_eff %>% 
  left_join(hr_losses) %>% 
  drop_na(losses) %>%
  mutate(lost = ifelse(losses > 0, "Loss", "No loss"),
         log_eff = log(efficiency)) %>% 
  filter(!is.infinite(log_eff) & !is.na(log_eff)) %>% 
  filter(season == "fall") %>% 
  ggplot(aes(log_eff, as.character(journey_number), fill = as.factor(losses))) +
  ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                panel_scaling = T, alpha = 0.85) +
  # scale_fill_manual(values = cols, name = "") +
  labs(x = "log per-hour soaring-gliding efficiency (m/s)", y = "Age") +
  theme(legend.position = "top")

ridge_pls2 <- lapply(c("spring", "fall"), function(s){
  lab <- ifelse(s == "spring", "Spring", "Fall")
  if(s == "spring"){
    lab2 <- c("C", "D")
  }else{
    lab2 <- c("A", "B")
  }
  if(s == "spring"){
    cols <- rev(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#becfb9"))
  }else{
    cols <- rev(c("#ee5e53", "#F07167", "#f27e71",  "#f7a58f", "#FED9B7"))
  }
 
  rp <- sg_eff %>% 
    left_join(hr_losses) %>% 
    drop_na(losses) %>%
    mutate(lost = ifelse(losses > 0, "Loss", "No loss"),
           log_eff = log(efficiency)) %>% 
    filter(!is.infinite(log_eff) & !is.na(log_eff)) %>% 
    filter(season == s) %>% 
    ggplot(aes(log_eff, as.character(journey_number), fill = as.factor(losses))) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                  panel_scaling = T, alpha = 0.85) +
    scale_fill_manual(values = cols, name = "# of losses") +
    labs(x = "log per-hour soaring-gliding efficiency (m/s)", y = "Age", title = lab) +
    theme(legend.position = "top")
  return(rp)
})

ggpubr::ggarrange(ridge_pls2[[2]], ridge_pls2[[1]])

# scale the predictors
mod_eff <- sg_eff %>% 
  left_join(hr_losses) %>% 
  drop_na(losses) %>%
  mutate(lost = losses > 0,
         log_eff = log(efficiency)) %>% 
  mutate(individual.id = as.factor(individual.id)) %>% 
  filter(!is.infinite(log_eff) & !is.na(log_eff)) %>% 
  dplyr::select(individual.id, season, lost, journey_number, log_eff, losses) %>%
  filter(between(log_eff, quantile(log_eff, 0.025), quantile(log_eff, 0.975))) %>% 
  mutate(age_z = scale(journey_number)[,1],
         log_eff_z = scale(log_eff)[,1],
         losses_z = scale(losses)[,1])

# thermals they lost had slightly lower vertical speeds and more directional changes
mod_eff %>% 
  ggplot(aes(log_eff, fill = as.character(lost))) +
  geom_density(alpha = 0.5) +
  labs(x = "Scaled value", y = "Density", fill = "Thermal loss") +
  facet_wrap(~season, ncol = 2)

effects <- lapply(c("fall", "spring"), function(s){
  ride_mod <- lme4::lmer(log_eff_z~lost+age_z+(1|individual.id), data = mod_eff[mod_eff$season == s,])
  lab <- round(performance::performance_rmse(ride_mod), 3)
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
    mutate(pred = gsub("age", "Age",
                       sub("log_eff", "Efficiency",
                           sub("lostTRUE", "Thermal loss",
                               gsub("_z", "", pred)))),
           sig = sign(lower)==sign(upper))
  
  ride_coefs <- ggplot(conf_limit_ride, aes(fixef, fct_reorder(pred, nchar(pred)))) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_linerange(aes(xmin = lower, xmax = upper), lwd = 1.75, color = col) +
    geom_point(size = 4, color = col) +
    scale_x_continuous(limits = c(-0.5, 0.01), breaks = seq(-0.5, 0.01, 0.2)) +
    labs(x = "Fixed effect estimate and 95% CI", y = "", color = "Significant", subtitle = lab)
  return(ride_coefs)
}) # RMSE: spring = 0.44, fall = 0.43

# create Figure 3
tgrob1 <- ggpubr::text_grob("Fall", size = 20)
tgrob2 <- ggpubr::text_grob("Spring", size = 20)
# Draw the text
plot_0 <- ggpubr::as_ggplot(tgrob1) + theme(plot.margin = margin(0,12,0,0, "cm"))
plot_1 <- ggpubr::as_ggplot(tgrob2) + theme(plot.margin = margin(0,12,0,0, "cm"))
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/hourly_efficiency_loss.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(ggpubr::ggarrange(plot_0, plot_1),
                  ggpubr::ggarrange(ridge_pls[[2]], ridge_pls[[1]], labels = "AUTO"),
                  ggpubr::ggarrange(effects[[1]], effects[[2]], labels = c("C", "D")),
                  nrow = 3, ncol = 1, heights = c(1, 5, 5), vjust = 6)
# dev.off()

