## script to perform shapely analysis across multiple species
library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(randomForest)

# data import locations
user    = 'cw21p621'
dd      = paste0('C:/Users/', user, '/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/')

# BASE INPUTS
## catchment areas of switzerland
## boundary area of switzerland

source('scripts/functions/spatial_data.R')

# INPUTS
## model = final random forest models fitted in the sdm-pipeline project scripts
## model_final_data = data used to fit the models in the sdm-pipeline project scripts 
## env_data = environmental data used to fit the models in the sdm-pipeline (constrained to species known distributions)

# PATHS TO INPUT FILES (each has 5 models/data sets)
sp_model <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/output/final_rf.rds'
sp_data  <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/data/model_data_final.rds'
sp_env   <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/data/env_data.TIF'
sp_pred_01  <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/output/presence_stack.TIF'

# READ IN ONE EXAMPLE SPECIES
sp_model <- readRDS(sp_model)
sp_data  <- readRDS(sp_data)
sp_env   <- rast(sp_env)

# MAKE NEW DATA OBJECT FOR SPECIES
# extract the environmental data values for the modelled variables
env_poly <- terra::extract(x = sp_env, 
                           y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script 
                           fun = function(x){mean(x, na.rm = T)})

# bind back in environmental data from extractions
subcatchments_final_env <- cbind(subcatchments_final, env_poly)

# drop geometrics and keep only catchments with data in the final objects
subcatchments_final_env_df <- st_drop_geometry(subcatchments_final_env) %>% na.omit()


#### Identify catchments with presences and absences predicted ----

# presence absence model prediction as raster
sp_pred_01 <- median(rast(sp_pred_01))

# extract catchment presence or absence
catchment_PA <- cbind(subcatchments_final_env %>% select(TEZGNR40), terra::extract(sp_pred_01, vect(subcatchments_final_env), function(x) median(x, na.rm = T)))

catchment_present = catchment_PA %>% st_drop_geometry() %>% filter(median == 1) %>% pull(TEZGNR40)
catchment_absent = catchment_PA %>% st_drop_geometry() %>% filter(median == 0) %>% pull(TEZGNR40)

#### Perform shapley function for a given set of data ----

# run shapleys over all model iteractions
sp_shapley <- lapply(1:length(sp_model), function(x){ 
  
  print(x)
  
  vars <- colnames(attr(sp_model[[x]]$terms,"factors"))
  
  get_shapleys(model = sp_model[[x]], # random forest model fitted to data
               data = sp_data[[x]],  # data used to fit the random forest
               new_data = subcatchments_final_env_df, # new data on which to predict the shapely values
               vars = vars
  )
  
  })


# bind together listed outputs
sp_shapley <- bind_rows(sp_shapley, .id = 'model_run')

# aggregate across all values
sp_shapley_summary <- sp_shapley %>% 
  select(TEZGNR40, grep('SHAP', colnames(.))) %>% 
  pivot_longer(., cols = 2:ncol(.)) %>% 
  group_by(TEZGNR40, name) %>% 
  do(value_mean = mean(.$value,na.rm = T), 
     value_sd   = sd(.$value, na.rm = T)) %>% 
  unnest()
  

#### Create maps of shapley values per variable ----


sp_tmap_shap <- lapply(unique(sp_shapley_summary$name), function(x){
  
  sp_shap_short <- sp_shapley_summary %>% filter(name == x)
  
  shap_spatial <- st_as_sf(left_join(sp_shap_short, subcatchments_final_env %>% select(TEZGNR40, Shape)), crs = target_crs)
  
  mean_shap <- tm_shape(shp = shap_spatial) + tm_polygons(col = 'value_mean', style = 'cont', palette = 'RdBu', title = "")  + tm_layout(title = paste('mean shapely of', x))
  sd_shap   <- tm_shape(shp = shap_spatial) + tm_polygons(col = 'value_sd', palette = 'viridis', style = 'cont', midpoint = NA, title = "") + tm_layout(title = paste('sd shapely of', x))
  
  return(tmap_arrange(list(mean_shap, sd_shap)))
  
})


pdf(file = 'figures/test_shaps.pdf', width = 10, height = 10)
lapply(sp_tmap_shap, print)
dev.off()


#### Rank positive and negative shapley values within catchments ----

pos_shapley <- sp_shapley_summary %>% 
  filter(value_mean >= 0.05) %>% 
  group_by(TEZGNR40) %>% 
  arrange(value_mean) %>% 
  nest() %>% 
  mutate(rank = purrr::map(data, ~order(.$value_mean, decreasing = T))) %>% 
  unnest(c(data, rank)) %>% 
  ungroup() %>% 
  filter(TEZGNR40 %in% catchment_present) %>% 
  left_join(., subcatchments_final_env %>% select(TEZGNR40, Shape)) %>% 
  st_as_sf(., crs = target_crs)

neg_shapley <- sp_shapley_summary %>% 
  filter(value_mean <= -0.05) %>% 
  group_by(TEZGNR40) %>% 
  arrange(value_mean) %>% 
  nest() %>% 
  mutate(rank = purrr::map(data, ~order(.$value_mean, decreasing = F))) %>% 
  unnest(c(data, rank)) %>% 
  ungroup() %>%  
  filter(TEZGNR40 %in% catchment_absent) %>% 
  left_join(., subcatchments_final_env %>% select(TEZGNR40, Shape)) %>% 
  st_as_sf(., crs = target_crs)


tm_shape(shp = pos_shapley %>% filter(rank < 3)) + tm_polygons(col = 'name') + tm_facets(by = 'rank') 
tm_shape(shp = neg_shapley %>% filter(rank < 3)) + tm_polygons(col = 'name') + tm_facets(by = 'rank') 


tm_shape(cropping_sf) +
  tm_polygons(col = 'gray95') + 
  tm_shape(shp = pos_shapley %>% 
             filter(name %in% c('ecoF_eco_meanSHAP', 
                                'ecoF_eco_1SHAP'))) + 
  tm_polygons(col = 'rank', palette = '-Oranges') + 
  tm_shape(ch_rivers) + 
  tm_lines() + 
  tm_shape(lakes_in_ch) + 
  tm_polygons(col = 'white') + 
  tm_layout(title = 'positive ecomorphological contributions')


tm_shape(cropping_sf) +
  tm_polygons(col = 'gray95') + 
  tm_shape(shp = neg_shapley %>% 
             filter(name %in% c('ecoF_eco_above_3SHAP', 
                                'ecoF_eco_meanSHAP'))) + 
  tm_polygons(col = 'rank', palette = '-Oranges') + 
  tm_shape(ch_rivers) + 
  tm_lines() + 
  tm_shape(lakes_in_ch) + 
  tm_polygons(col = 'white') + 
  tm_layout(title = 'negative ecomorphological contributions')







