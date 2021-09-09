#Author:Peter Laurin
#Date: Aug 30, 2021
#Glucosinulate Paper - Phenotypic Mismatch v. Euclidean Distance

# Peter Laurin
# Script to plot Haversine distances by the proportion of glucosinolate 
# classification mismatches, for every combination of accession in Katz et al.
# 2021. Commented out portion is basically what would be locally run,
# but due to computational intensity it was parallized under the geo_boss.sh 
# script



library(tidyverse)

#Katz_accession_data <- read_tsv("Katz_accessions_and_classifications.txt")

#h_dist <- function(CS1, CS2){
#  coords1 <- Katz_accession_data %>% filter(CS == CS1) %>% select(Longitude, Latitude) %>% as.numeric()
#  coords2 <- Katz_accession_data %>% filter(CS == CS2) %>% select(Longitude, Latitude) %>% as.numeric()
#  return(distHaversine(coords1, coords2))
#}


all_CS_combos <- combn(Katz_accession_data$CS, 2)

#distances <- mapply(h_dist, CS1 = all_CS_combos[1,], CS2 = all_CS_combos[2,])

#write_lines(distances, file = out_name)

#distances_total.txt written in geo_boss.sh
distances <- read_lines("distances_total.txt")

#find mismatches in combinations of accessionos
plotting_data <- all_CS_combos %>% t() %>% as_tibble() %>% rename(CS = V1, sample2 = V2) %>% mutate(distance = distances %>% as.numeric())
reduced_data <- Katz_accession_data %>% select(CS, Classification_name)
plotting_data <- left_join(plotting_data, reduced_data, by = "CS") %>% rename(one_classification = Classification_name, sample1 = CS, CS = sample2)
plotting_data <- left_join(plotting_data, reduced_data, by = "CS") %>% rename(two_classification = Classification_name, sample2 = CS)
plotting_data <- plotting_data %>% transmute(distance, mismatch = (one_classification != two_classification) %>% as.integer())

#distance bins for full range of distances
plotting_full <- plotting_data %>% mutate(dist_bin = cut_interval(distance, n = 100))
plotting_full <- plotting_full %>% group_by(dist_bin) %>% summarise(n = n(), n_mismatch = sum(mismatch), prop = n_mismatch / n)

#distance bins - 
plotting_2 <- plotting_data %>% filter(distance < 2000000)
plotting_2 <- plotting_2 %>% mutate(dist_bin = cut_interval(distance, n = 50))
plotting_2 <- plotting_2 %>% group_by(dist_bin) %>% summarise(n = n(), n_mismatch = sum(mismatch), prop = n_mismatch / n)

#plot bins of distances / mismatches - full range

#ggplot(plotting_full, aes(dist_bin %>% as.numeric(), prop)) + 
#  geom_point(alpha = 0.99) + 
#  geom_smooth(color = "red", se=F) + 
#  theme_classic() +
#  xlab("bin number") + ylab("proportion of mismatches in glucosinulate classification") + 
#  labs(title = "full range of Haverstine distances (0 - 74,000 km)") + 
#  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

#plot bins of distances / mismatches - < 2000 Km
ggplot(plotting_2, aes(dist_bin %>% as.numeric() * 40, prop)) + 
  geom_point(alpha = 0.99) + 
  geom_smooth(color = "red", se=F, method = 'loess') + 
  theme_classic() +
  xlab("km") + ylab("proportion of mismatches in glucosinulate classification") + 
  labs(title = "range of Haverstine distances") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))





