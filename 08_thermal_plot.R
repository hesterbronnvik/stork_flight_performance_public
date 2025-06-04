### Create Figure 1 (what thermals look like, thermal loss, 1 Hz vertical and wind speeds)
### Hester Brønnvik
### 2024-12-05
### hbronnvik@ab.mpg.de

library(tidyverse)
# two different schemes
colfunc <- colorRampPalette(c("#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53", "#EB3F33"))
colfunc2 <- colorRampPalette(c("#004E66", "#006D8F", "#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#9DD2C9"))#, "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  strip.background = element_rect(fill = "white", color = "white")))

wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}

# from 03_segment_flight.R
seg_fls <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)
# get a fledgling from the 2024 data
pre_migration <- readRDS(seg_fls[grepl(4042709187, seg_fls)])

# get wind estimates for that bird (from 04_estimate_wind.R)
wind2 <- readRDS("/home/hbronnvik/Documents/WS_data/wind_thermal_250deg_20s/4042709187_seg_wind_2024-11-08.rds") %>% 
  filter(date(timestamp) == "2024-07-23")
# select a thermal
# get those winds
wtherm2 <- wind2 %>% 
  filter(burst_id == 2963) %>% 
  mutate(wind_speed = wind_speed(windX, windY))
# and thos GPS
comp_therm2 <- pre_migration %>% 
  filter(individual.id == 4042709187 & burst_id == 2963) %>% 
  slice(55:n()) %>% 
  left_join(wtherm2 %>% dplyr::select(timestamp, wind_speed))

# also get a thermal for an adult along with its wind estimates
f <- list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s", full.names = T)[grepl("1176031140", list.files("/home/hbronnvik/Documents/WS_data/segmented_flight_250deg_20s"))]
# GPS
ind <- readRDS(f) 
therm <- ind %>% 
  filter(burst_id == 16308) %>% 
  slice(1:250)

# wind
wind <- readRDS("/home/hbronnvik/Documents/WS_data/wind_thermal_250deg_20s/1176031140_seg_wind_2024-11-07.rds")
wtherm <- wind %>% 
  filter(burst_id == 16308) %>% 
  mutate(wind_speed = wind_speed(windX, windY))

comp_therm <- therm %>% 
  left_join(wtherm %>% dplyr::select(timestamp, wind_speed))

difftime("2024-02-25", "2020-06-29")

# stick the individuals together as a list
therms <- comp_therm %>% 
  dplyr::select(individual.id, flight_clust_sm3, location.lat, height.above.ellipsoid, vert.speed, wind_speed) %>% 
  rbind(comp_therm2 %>% 
          dplyr::select(individual.id, flight_clust_sm3, location.lat, height.above.ellipsoid, vert.speed, wind_speed)) %>% 
  pivot_longer(cols = contains("speed"), names_to = "variable", values_to = "speed") %>% 
  group_by(individual.id, variable) %>% 
  group_split()

# summarize the thermals
avgs <- comp_therm %>% 
  dplyr::select(individual.id, flight_clust_sm3, location.lat, height.above.ellipsoid, vert.speed, wind_speed) %>% 
  rbind(comp_therm2 %>% 
          dplyr::select(individual.id, flight_clust_sm3, location.lat, height.above.ellipsoid, vert.speed, wind_speed)) %>% 
  filter(flight_clust_sm3 == "circular_soaring") %>% 
  pivot_longer(cols = contains("speed"), names_to = "variable", values_to = "speed") %>% 
  mutate(v = ifelse(variable == "vert.speed", "thermal strength", "wind speed"))  %>% 
  group_by(individual.id, variable) %>% 
  summarize(std = sd(speed, na.rm = T),
            avg = paste0(unique(v), " = ", 
                         round(mean(speed, na.rm = T), 2),
                         " +/- ", round(std, 2)))

# make plots of the thermals and their vertical speeds
pls <- lapply(therms, function(therm){
  if(unique(therm$individual.id) == 1176031140){
    lab <- "Thermal loss (1,336 days since fledging)"
  }else{
    lab <- "Thermal soaring (day of fledging)"
  }
  # therm <- therms[[1]]
  if(unique(therm$variable) == "wind_speed"){
    p <- therm %>% 
      mutate(flight = ifelse(flight_clust_sm3 == "circular_soaring", "Soaring", "Other"),
             flight = factor(flight, levels = c("Soaring", "Other"))) %>% 
      ggplot(aes(x = location.lat, y = height.above.ellipsoid)) +
      geom_point(aes(color = speed, size = flight)) +
      geom_path(aes(color = speed)) +
      scale_size_manual(values = c(2.5, 1)) +
      # xlim()
      scale_color_gradientn(colors = colfunc(30)) +
      labs(x = "Latitude", y = "Height above ellipsoid",
           color = "Wind\nspeed", size = "Flight", title = lab) + 
      theme(legend.position = "top",
            plot.title = element_text(hjust = 0.5),
            legend.key.size = unit(0.75, "cm"),
            legend.title = element_text(size = 11),
            axis.title.x=element_text(color = "white"),
            axis.text.x=element_text(color = "white"),
            axis.ticks.x=element_line(color = "white"),
            axis.line.x=element_line(color = "white"))
  }else{
    txt <- avgs %>% 
      filter(individual.id == unique(therm$individual.id)) %>% 
      dplyr::select(avg) %>% 
      deframe()
    xloco <- ifelse(unique(therm$individual.id) == 1176031140, 47.5868, 48.0762)
    p <- ggplot(therm, aes(x = location.lat, y = speed)) +
      geom_hline(yintercept = 0) +
      geom_segment(aes(x = location.lat, y = 0, 
                       xend = location.lat, yend = speed,
                       color = speed), 
                   arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
      scale_color_gradientn("Thermal\nstrength", colors = colfunc2(30),
                            values = scales::rescale(c(min(therm$speed),
                                                       0,
                                                       max(therm$speed)))) + 
      scale_y_continuous(limits = c(-10.2, 6.3)) +
      annotate(geom = "text", x = xloco, y = c(-5, -7, -6), label = c("Mean:", txt), hjust = 0) +
      labs(x = "Latitude", y = "Speed", size = "Soaring") +
      theme(legend.position = c(.8,.2), legend.direction="horizontal",
            axis.title.y=element_text(color = "white"),
            axis.text.y=element_text(color = "white"),
            axis.ticks.y=element_line(color = "white"),
            axis.line.y=element_line(color = "white"),
            legend.key.size = unit(0.75, "cm"),
            legend.title = element_text(size = 11))
  }
  if(unique(therm$individual.id) == 1176031140){
    p <- p + scale_x_reverse()
  }
  return(p)
})
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/november/thermals_arrows.png",
#     height = 8.2, width = 11.7, units = "in", res = 300)
ggpubr::ggarrange(ggpubr::ggarrange(pls[[4]], pls[[2]], labels = c("A", "B")),
                  ggpubr::ggarrange(pls[[3]], pls[[1]], labels = c("C", "D")),
                  nrow = 2, heights = c(2, 1.5))
# dev.off()

