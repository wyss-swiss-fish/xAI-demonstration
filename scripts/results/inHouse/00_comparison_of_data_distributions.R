#### Compare inhouse records to kanton bern records

pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

# load in all spatial objects for making maps and plots
source('scripts/results/00-load-spatial.R')

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_JULY_inHouse_V1/"

# get run to mak figures for
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

