#Author:Peter Laurin
#Date: Aug 30, 2021
#Glucosinulate Paper - Phenotypic Mismatch v. Euclidean Distance

# R Script intended to work in parallel as part of a large dataset
# based on array num, which subsets the data into 10000-size chunks
# Haversine formula shown below, also calculated w/ geosphere distHaversine
# Katz_accessions_and_classifications is sheet "D. GSLs_data emmeans" from 
# https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjc3ODQvZWxpZmUtNjc3ODQtc3VwcDEtdjIueGxzeA--/elife-67784-supp1-v2.xlsx?_hash=DDyP7JsX8%2FBjy610ZvQD%2BgQAMe3BXtizrPYYHdByTmE%3D


library(tidyverse)
library(geosphere)

args <- commandArgs(trailingOnly = TRUE)
array_num <- args[1] %>% as.integer()
out_name <- args[2]

Katz_accession_data <- read_tsv("Katz_accessions_and_classifications.txt")

Haversine_dist <- function(CS1, CS2){
  coords1 <- Katz_accession_data %>% filter(CS == CS1) %>% select(Longitude, Latitude) %>% as.numeric()
  coords2 <- Katz_accession_data %>% filter(CS == CS2) %>% select(Longitude, Latitude) %>% as.numeric()
  earth_rad <- 6371 #km
  diff_lon <- coords1[1] - coords2[1]
  lat1 <- coords1[2]
  lat2 <- coords2[2]
  diff_lat <- lat1 - lat2
  a <- (sin(diff_lat / 2))^2 + 
    cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * 
    (sin(diff_lon / 2))^2
  c <- 2 * earth_rad * pmin(1, sqrt(a))
  return(earth_rad * c)
}

h_dist <- function(CS1, CS2){
  coords1 <- Katz_accession_data %>% filter(CS == CS1) %>% select(Longitude, Latitude) %>% as.numeric()
  coords2 <- Katz_accession_data %>% filter(CS == CS2) %>% select(Longitude, Latitude) %>% as.numeric()
  return(distHaversine(coords1, coords2))
}

reduced_data <- Katz_accession_data %>% select(CS, Latitude, Longitude, Classification_name)
all_CS_combos <- combn(Katz_accession_data$CS, 2)
subs <- seq(to = dim(all_CS_combos)[2], by = 10000)
beg <- subs[array_num]
if(array_num == 32){
  end <- dim(all_CS_combos)[2]
} else {
  end <- subs[array_num + 1] - 1
}

distances <- mapply(h_dist, CS1 = all_CS_combos[1,beg:end], CS2 = all_CS_combos[2,beg:end])

write_lines(distances, file = out_name)





