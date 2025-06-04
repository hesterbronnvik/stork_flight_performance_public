# stork_flight_performance_public
Scripts and analysis for the ontogeny of white stork flight performance study

# Abstract
Movement allows animals to access resources and pursue fitness. Yet, this requires balancing movement costs against other potential benefits. These trade-offs can change over lifetimes, as optimality is contextual and depends on ability and need. For migratory soaring birds, efficient movement requires exploiting atmospheric uplift---a challenging task that negates the cost of flight. Due to this link between costs and the environment, soaring birds are limited in when and where to fly. Thus, soaring performance goes beyond movement skills and also involves coping with potentially imperfect conditions to reduce time or competition. Here we asked whether, over repeated migrations ,white storks (Ciconia ciconia) improve their ability to soar. Using high-resolution lifetime tracking data from 151 storks, we found that juveniles outperformed adults under supportive conditions but with age migration difficulty increased and adults performed well under challenging conditions. Adults traveled in less supportive conditions and spent more energy on flight, indicating a change in ability and motivation. Thus, understanding how animals improve is not simply described by a learning curve, but requires a multifaceted perspective on individual needs and skills.

# This repository consists of the R scripts:

01_access_data.R: download the data from Movebank.org

02_segment_tracks.R: identify migration days in the GPS data.

03_segment_flight.R: identify times when the storks were in flight and classify flight types.

04_estimate_wind.R: following Weinzierl et al. 2016 (https://doi.org/10.1002/ece3.2585), estimate 1 Hz wind speeds.

05_estimate_flapping.R: following A. Scharf (https://doi.org/10.3389/fevo.2019.00200), estimate instances of flapping flight in the ACC data.

06_estimate_thermal_loss.R: identify times when storks fell lost thermals and extract metrics of thermal loss.

07_estimate_efficiency.R: estimate soaring-gliding efficiency.

08_thermal_3d_plot.R: visualize the data.

09_performance_pre_migration.R: explore the data from juveniles tagged in 2024 and transmitting high-resolution data from before migration.

All raw input data are available on Movebank.org
