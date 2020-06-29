## Project: Causal discovery in heavy-tailed data
## Descrip: Produce map of the upper Danube basin
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

## Import libraries ####
rm(list = ls())
library(tidyverse)
library(ggmap)
library(ggsn) # for scale bars/north arrows in ggplots

## Define functions ####
clean_river_data <- function(dfr){
  ## data_frame -> data_frame
  ## reshape each row of the dataframe

  ## helper
  clean_river_data_helper <- function(list_of_2){
    ## list -> list
    ## reshape each row of river_data

    ## helper
    split_str <- function(str, sep, invert=F){
      ## character character boolean -> character_v
      ## split string into 2 parts using separator sep and return it as a vector
      s <- unlist(strsplit(str, sep))
      if (length(s) == 1){
        s <- c("", s)

        if (invert){
          s <- s[2:1]
        }
      }
      return(s)
    }
    ## function body

    s1 <- split_str(list_of_2[[1]], "\t", F)
    s2 <- split_str(list_of_2[[2]], "\t", T)

    list_of_4 <- list(name = s1[1],
                      lat = as.numeric(s1[2]),
                      lon = as.numeric(s2[1]),
                      stat = as.numeric(s2[2]))
    return(list_of_4)
  }

  ## function body
  dfr <- transpose(dfr)
  map_dfr(.x = dfr, .f = clean_river_data_helper)
}

## Compute river coordinates ####
file_names <- c("Danube", "Iller", "Inn", "Isar",
                "Lech", "Naab", "Regen", "Saalach", "Saalzach")

rivers <- tibble(lon = double(), lat = double(),
                 group = integer(), name = character())
for (i in 1:length(file_names)){
  fname <- paste("data/river_data/river_coords/",
                 file_names[i], ".txt", sep = "")

  raw_dat <- read_csv(fname, col_names = F)
  clean_dat <- clean_river_data(raw_dat) %>%
    select(lon, lat) %>%
    mutate(group = i) %>%
    mutate(name = file_names[i])

  rivers <- rivers %>%
    bind_rows(clean_dat)
}

## Import station coordinates ####
load("data/river_data/StsInfo.RData")

# Select the 31 stations that are in the AOAS paper plus station 19 (renamed 32)
StsChos <- c(c(1:47)[-c(16,30,31,34,42,43,44,45,46,47,3,1,2,29,18,19)], 19)
NoSt <- length(StsChos)
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)

station_info <- StsInfo[StsChos,] %>%
  mutate(id_old = StsChos,
         id = 1:NoSt) %>%
  filter(id %in% station_names) %>%
  rename(name = RivNames,
         lat = Lat,
         lon = Long,
         ave_vol = AveVol)

## Download map ####
bbox = c(9.9, 47.40, 13.60, 49.40)
map_2_plot <- get_stamenmap(bbox = bbox,
                            maptype = "terrain-background",
                            zoom = 9,
                            color = "color",
                            force = TRUE
                            )

## Plot map ####
tolBlue <- "#4477AA"
tolRed <-  "#EE6677"

lbls <- tibble(name = c("Danube", "Iller", "Lech", "Danube", "Isar", "Inn", "Salzach"),
               lon = c(10.1, 10.2, 11, 12.7, 12.35, 12.15, 13),
               lat = c(48.5, 48.2, 48.2, 49, 48.7, 48.15, 48))

size_stations <- 3
size_names <- 3

output_map <- ggmap(map_2_plot) +
  labs(x="Longitude", y="Latitude") +
  geom_path(data = rivers, aes(x = lon, y = lat,
                                     group = group, fill = NULL),
            color = tolBlue, alpha = .8, size = 1.5) +
  geom_point(data = station_info, aes(x = lon, y = lat),
             pch=21, col="black", bg=tolRed, size = 5, alpha = 1) +
  geom_label(data = station_info %>% filter(!(id %in% c(1, 13))),
                                           aes(x=lon, y=lat, label=id),
             nudge_x=-0.1, nudge_y=0.1, label.size=0.1,size=size_stations,
             fontface = "bold.italic", label.r=unit(0.20, "lines")) +
  geom_label(data = station_info %>% filter(id == 1),
             aes(x=lon, y=lat, label=id),
             nudge_x=0, nudge_y=0.1, label.size=0.1,size=size_stations,
             fontface = "bold.italic", label.r=unit(0.20, "lines")) +
  geom_label(data = station_info %>% filter(id == 13),
             aes(x=lon, y=lat, label=id),
             nudge_x=-0.2, nudge_y=0, label.size=0.1, size=size_stations,
             fontface = "bold.italic", label.r=unit(0.20, "lines")) +
  geom_text(data = lbls, aes(x = lon, y = lat, label = name),
            angle = c(20, 90, 90, -30, 20, 45, -60),
            size=size_names, fontface = "bold.italic")


output_map <- output_map +
  scalebar(x.min = attr(map_2_plot, "bb")[[2]],
           y.min=attr(map_2_plot, "bb")[[1]],
           x.max =attr(map_2_plot, "bb")[[4]],
           y.max=attr(map_2_plot, "bb")[[3]],
           dist = 25, anchor = c(x=10.1, y=49.2),
           model = 'WGS84', transform = T,
           location = "topleft", st.size = 3.3, st.dist = 0.02,
           dist_unit = "km"); output_map

ggsave("output/river_map.pdf", output_map,
       width = 6, height = 5, units = c("in"))
