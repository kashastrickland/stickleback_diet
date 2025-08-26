# Load packages
library(tidyverse)
# library(lemon)
# library(cowplot)
# library(zoo)
library(sf)
library(ggspatial)



# import shape file and station coordinates
shape_file <- st_read("data/Myvatn_WSGUTM28.shp")
stations <- read_csv("data/stations.csv") %>%
  filter(station %in% c(23, 44, 135, 124, "NS" ,"GS", "HS"))

shape_file_coords <- as_tibble(st_coordinates(shape_file)) %>%
  rename(lat = Y,
         long = X)

# adjust stations
stations_adj <- stations 
stations_adj[stations_adj$station == "NS",]$lat <- 
  stations[stations$station == "NS",]$lat + 400
stations_adj[stations_adj$station == "GS",]$lat <- 
  stations[stations$station == "GS",]$lat + 400
stations_adj[stations_adj$station == "135",]$lat <- 
  stations[stations$station == "135",]$lat + 75
stations_adj[stations_adj$station == "135",]$long <- 
  stations[stations$station == "135",]$long - 30
stations_adj[stations_adj$station == "HS",]$long <- 
  stations[stations$station == "HS",]$long + 100

# specify station groupings
stations_adj$station_group <- ifelse(stations_adj$station %in% c(23, 44, 135), "South basin lake",
                                     ifelse(stations_adj$station %in% c("NS", "GS"), "North basin shore",
                                            ifelse(stations_adj$station == "124", "North basin lake", 
                                                   "Warm springs")))
# create vector of groups and group colors
station_groups <- c("North basin lake",
                    "North basin shore",
                    "South basin lake",
                    "Warm springs")
station_group_colors <- c("#cc79a7","darkgoldenrod1","#0072b2","#d55e00")


# plot
fig_map <- ggplot(data = shape_file_coords %>% filter(L1 == 1, L2 == 1),
                  aes(long,lat))+
  geom_polygon(aes(group = interaction(L1, L2)), 
               size = 0.3, color = "black", fill = "gray95")+
  geom_polygon(data = shape_file_coords %>% filter(L1 != 1),
               aes(group = interaction(L1, L2)), 
               size = 0.3, color = "black",  fill = "gray100")+
  coord_equal()+
  geom_text(data = stations_adj,
            aes(label = station,
                color = station_group),
            fontface = "bold",
            size = 3)+
  scale_color_manual("", 
                     values = station_group_colors,
                     labels = NULL)+
  guides(color = guide_legend(
    override.aes = list(label = station_groups,
                        color = station_group_colors,
                        shape = NA),
    title = NULL
  ))+
  theme_void()+
  theme(plot.margin = margin(0,0,0,0))+
  annotation_scale(location = "tl", 
                   height = unit(0.075, "in"),
                   width_hint = 0.4,
                   pad_x = unit(0.01, "in"),
                   pad_y = unit(0.15, "in"),
                   text_cex = 0.8,
                   bar_cols = c("gray60", "white")) + 
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering(fill = c("gray60", "white"),
                                                          text_size = 7),
                         pad_x = unit(0.28, "in"),
                         pad_y = unit(0.3, "in"),
                         height = unit(0.3, "in"),
                         width = unit(0.275, "in"))
fig_map

# export
ggsave("tables_figures/fig_map.pdf",
        plot = fig_map,
        width = 3.5,
        height = 3.5)

