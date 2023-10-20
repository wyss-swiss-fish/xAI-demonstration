#### Compare inhouse records to kanton bern records

pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# load in all spatial objects for making maps and plots
source('scripts/results/00-load-spatial.R')

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_JULY_inHouse_V1/"

# get run to make figures for
RUN <- "ubelix_SDM_RF_JULY_inHouse_V1"

sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# read in data
final_data <- paste0(dd, 'sdm-pipeline/species-records-final/fish-presenceAbsence_2010.csv')
final_data <- read_csv(final_data)

# get unique locations
final_data_unique <- st_as_sf(final_data, coords = c('X', 'Y'), crs = 'wgs84') %>% 
  select(dataset, occ, species_name) %>% 
  unique 

final_data_unique$dataset <- recode(final_data_unique$dataset, 
                                    'alte_aare_monitoring' = 'kanton', 
                                    bafi_data = 'kanton', 
                                    lanat_3_unibe = 'EAWAG', 
                                    project_fiumi = 'EAWAG', 
                                    vonlanthen = 'consultancy', 
                                    'wa_sampling' = 'EAWAG')

# map
tm_shape(final_data_unique %>% arrange(occ)) + 
  tm_dots(size = 2, col = 'occ') + 
  tm_facets(by = c('dataset'), free.coords = FALSE) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01)

for(i in sp_list){
print(tm_shape(final_data_unique %>% filter(species_name == i) %>% arrange(occ)) + 
  tm_dots(size = 2, col = 'occ') + 
  tm_facets(by = c('dataset'), free.coords = FALSE) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01))
}

# map
pdf('figures/data_distribution.pdf', width = 10, height = 10)
tm_shape(final_data_unique %>% 
           filter(species_name %in% sp_list[1:3]) %>% 
           arrange(occ)) + 
  tm_dots(size = 0.1, col = 'occ') + 
  tm_facets(by = c('dataset', 'species_name'), free.coords = FALSE) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01)
tm_shape(final_data_unique %>% 
           filter(species_name %in% sp_list[4:8]) %>% 
           arrange(occ)) + 
  tm_dots(size = 0.1, col = 'occ') + 
  tm_facets(by = c('dataset', 'species_name'), free.coords = FALSE) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01)
dev.off()




## map the different subsets of the data ----

### In house EAWAG and NAWA data
data_names <- c('fish-presenceAbsence_2010.csv', 'fish-presenceAbsence_2010_inHouse.csv', 'fish-presenceAbsence_2010_Kanton.csv', 'fish-presenceAbsence_2010_Kanton_subset.csv')
dataset_locations <- paste0(dd, 'sdm-pipeline/species-records-final/', data_names)

# load in datasets and process for mapping
all_datasubsets <- lapply(1:4, function(x){
  read.csv(dataset_locations[x]) %>% 
    st_as_sf(., coords = c('X', 'Y'), crs = 'wgs84') %>% 
    select(dataset) %>% 
    unique() %>% 
    mutate(dataset_2 = gsub('.csv', '', data_names[x]))
}) %>% 
  bind_rows()

all_datasubsets$dataset <- recode(all_datasubsets$dataset, 
                                    'alte_aare_monitoring' = 'kanton', 
                                    bafi_data = 'kanton', 
                                    lanat_3_unibe = 'EAWAG', 
                                    project_fiumi = 'EAWAG', 
                                    vonlanthen = 'consultancy', 
                                    'wa_sampling' = 'EAWAG')

# to better visualise
all_datasubsets <- all_datasubsets[sample(1:nrow(all_datasubsets), replace = F), ]


# map out difference subsets
pdf('figures/data_distribution_allsurveys.pdf', width = 10, height = 6)
tm_shape(all_datasubsets) + 
  tm_dots(size = 0.2, col = 'dataset', palette = 'Set1') +
  tm_facets(by = 'dataset_2', free.coords = FALSE)+
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.5)
dev.off()

### Excluding data from Kander, Lutschen, Simme and Upper Aare that are naturally disconnected 


