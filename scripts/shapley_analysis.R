## script to perform shapely analysis across multiple species
library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(randomForest)

# data import locations
dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'
dd_env <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/'

# BASE INPUTS
## catchment areas of switzerland
## boundary area of switzerland

source('scripts/functions/spatial_data.R')
source('scripts/functions/shap_testing.R')

# INPUTS
## model = final random forest models fitted in the sdm-pipeline project scripts
## model_final_data = data used to fit the models in the sdm-pipeline project scripts 
## env_data = environmental data used to fit the models in the sdm-pipeline (constrained to species known distributions)

# get list of files paths list
dirs <- list.dirs(path = paste0(dd, 'sdm-pipeline/sdm-run/exAI'), recursive=T)
dirs_df <- data.frame(str_split(dirs, '/', simplify = T))
dirs_df$path <- dirs

# filter to options of interest
dirs_df <- dirs_df %>% filter(!X14 %in% c('', 'figures'), 
                              X11 != 'VSURF_2010')


# split by groupings of interest
dirs_list <- split(dirs_df, paste0(dirs_df$X11,'_', dirs_df$X13))

# define focal catchments
urtene  <- 101012
limpach <- 106961
catchment_shaps <- c(urtene = urtene, limpach = limpach)

### loop through the list of options
for(i in 1:length(dirs_list)){
  
  print(i)
  # SELECTION OF LOOPED SPECIES AND MODELS
  dir_i    <- dirs_list[[i]]
  dir_out  <- dir_i$path[2] 
  dir_data <- dir_i$path[1]
  run_type <- unique(dir_i$X11)
  sp_name  <- unique(dir_i$X13)
  
  # CREATE DIRECTORY TO OUTPUT FIGURES
  dir.create(paste0('figures/', run_type, '/', sp_name, '/'), recursive = T)
  
  # PATHS TO INPUT FILES (each has 5 models/data sets)
  sp_model <- paste0(dir_out, '/final_rf.rds')
  sp_data  <- paste0(dir_data, '/model_data_final.rds')
  sp_env   <- paste0(dir_data, '/env_data.TIF')
  sp_pred_01  <- paste0(dir_out, '/presence_stack.TIF')
  sp_suit  <- paste0(dir_out, '/suitability_stack.TIF')
  
  # READ IN ONE EXAMPLE SPECIES
  sp_model <- readRDS(sp_model)
  sp_data  <- readRDS(sp_data)
  sp_env   <- rast(sp_env)
  sp_pred_01 <- median(rast(sp_pred_01))
  sp_suit  <- mean(rast(sp_suit))
  sp_suit[sp_pred_01 == 0] <- NA
  
  
  # CREATE SUITABILITY POLYGONS
  ## extract suitability across polygons
  suit_poly <- terra::extract(x = sp_suit, 
                              y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script 
                              fun = function(x){mean(x, na.rm = T)})
  
  # bind back in environmental data from extractions
  subcatchments_final_suit <- cbind(subcatchments_final, suit_poly)
  subcatchments_final_suit_df <- st_drop_geometry(subcatchments_final_suit) %>% na.omit()
  
  # plot habitat suitability
  png(filename = paste0('figures/', run_type, '/', sp_name, '/', 'suitability_map.png'), width = 2500, height = 2500, res = 300)
  print(tm_shape(cropping_sf) +
    tm_polygons(col = 'gray95') +
    tm_shape(subcatchments_final_suit) + 
    tm_polygons(col = 'mean', palette = 'viridis', style = 'cont', breaks=c(0.5,0.6,0.8,1), 
                border.col = 'transparent', midpoint = NA, title = "", legend.reverse = T) +
    tm_shape(ch_rivers_2) +
    tm_lines(lwd = 'CUM_LEN_LOG', scale=2, col = 'gray80', legend.lwd.show = F) + 
    tm_shape(lakes_in_ch) + 
    tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
    tm_layout(title = 'habitat suitability'))
  dev.off()
  
  
  # MAKE NEW DATA OBJECT FOR SPECIES
  ## extract the environmental data values for the modelled variables
  env_poly <- terra::extract(x = sp_env, 
                             y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script 
                             fun = function(x){mean(x, na.rm = T)})
  
  
  # bind back in environmental data from extractions
  subcatchments_final_env <- cbind(subcatchments_final, env_poly)
  
  # drop geometrics and keep only catchments with data in the final objects
  subcatchments_final_env_df <- st_drop_geometry(subcatchments_final_env) %>% na.omit()
  
  
  #### IDENTIFY CATCHMENTS WITH PRESENCES AND ABSENCES PREDICTED ----
  
  # extract catchment presence or absence
  catchment_PA <- cbind(subcatchments_final_env %>% select(TEZGNR40), terra::extract(sp_pred_01, vect(subcatchments_final_env), function(x) median(x, na.rm = T)))
  
  catchment_present = catchment_PA %>% st_drop_geometry() %>% filter(median == 1) %>% pull(TEZGNR40)
  catchment_absent = catchment_PA %>% st_drop_geometry() %>% filter(median == 0) %>% pull(TEZGNR40)
  

  #### PERFORM SHAPLEY FUNCTION FOR GIVEN SET OF DATA ----

  # run shapleys over all model iterations
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
  
###  # save shapley values per species run
###  dir.create(paste0('data/', run_type, '/', sp_name, '/'), recursive = T)
###  saveRDS(sp_shapley_summary, file = paste0('data/', run_type, '/', sp_name, '/shapely_values.RDS'))
###  
  
     #### PLOT THE SHAPLEY VALUES FOR SPECIFIC CATCHMENTS ----
  
  if(!is.null(catchment_shaps)){
    
    png(filename =  paste0('figures/', run_type, '/', sp_name, '/', 'local_shapley.png'), res = 300, width = 1000, height = 2000)
    print(ggplot(data = sp_shapley_summary %>% 
             filter(TEZGNR40 %in% catchment_shaps, 
                    name %in% c('ecoF_slope_medianSHAP', 'stars_n_ch_m_cSHAP', 'ecoF_discharge_log_maxSHAP', 
                                'ecoF_eco_above_3SHAP', 'local_tcdSHAP'),) %>% 
             mutate(TEZGNR40 = plyr::revalue(.$TEZGNR40, setNames(names(catchment_shaps), catchment_shaps)), 
                    name = plyr::revalue(.$name, c(ecoF_slope_medianSHAP = 'slope', 
                                                   stars_n_ch_m_cSHAP = 'upstream nitrogen input', 
                                                   ecoF_discharge_log_maxSHAP = 'discharge', 
                                                   ecoF_eco_above_3SHAP = 'eco-morphological modification', 
                                                   local_tcdSHAP        = 'bankside tree cover')))) + 
      geom_bar(aes(x = value_mean, y = name, fill = value_mean), col = 'black', stat = 'identity') + 
      geom_vline(aes(xintercept = 0)) +
      facet_wrap(~TEZGNR40, nrow = 2) + 
      scale_fill_gradient2() + 
      theme_bw() + 
      theme(aspect.ratio = 1, 
            panel.grid = element_blank(), 
            legend.position = 'none', 
            strip.background = element_blank(), 
            strip.text = element_text(hjust = 0, size = 15)) + 
      xlab('shapley value') + 
      ylab(NULL))
      dev.off()
      
  }
  
}
  
###  #### MAP THE MEAN AND SD OF SHAPLEY VALUES ----
  
###  # perform spatial join for plotting
###  shap_spatial <- st_as_sf(left_join(sp_shapley_summary, subcatchments_final_env %>% select(TEZGNR40, Shape)), crs = target_crs)
###  
###  # map the shaps
###  sp_tmap_shap <- lapply(unique(sp_shapley_summary$name), function(x){
###    
###    sp_shap_short <- shap_spatial %>% filter(name == x)
###    
###    mean_shap <- tm_shape(cropping_sf) +
###      tm_polygons(col = 'gray95') +
###      tm_shape(shp = sp_shap_short) + 
###      tm_polygons(col = 'value_mean', style = 'cont', palette = 'RdBu', border.col = 'transparent', title = "", legend.reverse = TRUE)  + 
###      tm_layout(title = paste('mean shapely of', x)) + 
###      tm_shape(ch_rivers_2) +
###      tm_lines(lwd = 'CUM_LEN_LOG', scale=2, col = 'gray80', legend.lwd.show = F) + 
###      tm_shape(lakes_in_ch) + 
###      tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
###      tm_layout(title = x)
###    
###    sd_shap   <- tm_shape(cropping_sf) +
###      tm_polygons(col = 'gray95') +
###      tm_shape(shp = sp_shap_short) + 
###      tm_polygons(col = 'value_sd', palette = 'viridis', style = 'cont', border.col = 'transparent', midpoint = NA, title = "") + 
###      tm_layout(title = paste('sd shapely of', x)) + 
###      tm_shape(ch_rivers_2) +
###      tm_lines(lwd = 'CUM_LEN_LOG', scale=2, col = 'gray80', legend.lwd.show = F) + 
###      tm_shape(lakes_in_ch) + 
###      tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
###      tm_layout(title = x)
###    
###    return(tmap_arrange(list(mean_shap, sd_shap)))
###    
###  })
###  
###  # give output names
###  names(sp_tmap_shap) <- unique(sp_shapley_summary$name)
###  
###  # save all as pngs
###  dir.create(paste0('figures/', run_type, '/', sp_name, '/', 'shapley_values'), recursive = T)
###  lapply(names(sp_tmap_shap), function(x){
###    
###    png(file = paste0('figures/', run_type, '/', sp_name, '/', 'shapley_values', '/', x, '.png'), res = 300, height = 2000, width = 2000)
###    print(sp_tmap_shap[[x]])
###    dev.off()
###    
###  })
  
  
  #### SHAPLEY PLOTS WITH PRESENCE ABSENCE MASKS ----
  
  # perform spatial join for plotting
  shap_spatial <- st_as_sf(left_join(sp_shapley_summary, subcatchments_final_env %>% select(TEZGNR40, Shape)), crs = target_crs)
  shap_spatial_P <- filter(shap_spatial, TEZGNR40  %in% catchment_present)
  shap_spatial_A <- filter(shap_spatial, TEZGNR40  %in% catchment_absent)
  
  # map the shaps
  sp_tmap_shap_PA <- lapply(unique(sp_shapley_summary$name), function(x){
    
    sp_shap_short_P <- shap_spatial_P %>% filter(name == x, value_mean > 0.01)
    sp_shap_short_A <- shap_spatial_A %>% filter(name == x, value_mean < -0.01)
    
    mean_shap_P <- tm_shape(cropping_sf) +
      tm_polygons(col = 'white') +
      tm_shape(shp = sp_shap_short_P) + 
      tm_polygons(col = 'value_mean', style = 'cont', palette = 'Purples', border.col = 'transparent', title = "", legend.reverse = TRUE,
                  contrast = c(0.25, 1))  + 
      tm_layout(title = paste('mean shapely of', x)) + 
      tm_shape(ch_rivers_2) +
      tm_lines(lwd = 'CUM_LEN_LOG', scale=1.5, col = 'gray80', legend.lwd.show = F) + 
      tm_shape(lakes_in_ch) + 
      tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
      tm_layout(title = x)
    
    mean_shap_A <- tm_shape(cropping_sf) +
      tm_polygons(col = 'white') +
      tm_shape(shp = sp_shap_short_A) + 
      tm_polygons(col = 'value_mean', style = 'cont', palette = 'Reds', border.col = 'transparent', title = "", legend.reverse = TRUE, 
                  contrast = c(0.25, 1))  + 
      tm_layout(title = paste('mean shapely of', x)) + 
      tm_shape(ch_rivers_2) +
      tm_lines(lwd = 'CUM_LEN_LOG', scale=1.5, col = 'gray80', legend.lwd.show = F) + 
      tm_shape(lakes_in_ch) + 
      tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
      tm_layout(title = x)
    
    return(tmap_arrange(list(mean_shap_P, mean_shap_A)))
    
  })
  
  # give output names
  names(sp_tmap_shap_PA) <- unique(sp_shapley_summary$name)
  
  # remove those with no catchments
  missing_PA <- sapply(sp_tmap_shap_PA, is.null, simplify = T)
  if(sum(missing_PA) != 0){sp_tmap_shap_PA <- sp_tmap_shap_PA[-which(missing_PA)]}
  
  # save all as pngs
  dir.create(paste0('figures/', run_type, '/', sp_name, '/', 'shapley_values_PA'), recursive = T)
  lapply(names(sp_tmap_shap_PA), function(x){
    
    png(file = paste0('figures/', run_type, '/', sp_name, '/', 'shapley_values_PA', '/', x, '.png'), res = 300, height = 5000, width = 2500)
    print(sp_tmap_shap_PA[[x]])
    dev.off()
    
  })
    
  #### RANK POSITIVE AND NEGATIVE SHAPLEY VALUES WITHIN A GIVEN CATCHMENT ----
  
  pos_shapley <- sp_shapley_summary %>% 
    filter(value_mean >= 0.01) %>% 
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
    filter(value_mean <= -0.01) %>% 
    group_by(TEZGNR40) %>% 
    arrange(value_mean) %>% 
    nest() %>% 
    mutate(rank = purrr::map(data, ~order(.$value_mean, decreasing = F))) %>% 
    unnest(c(data, rank)) %>% 
    ungroup() %>%  
    filter(TEZGNR40 %in% catchment_absent) %>% 
    left_join(., subcatchments_final_env %>% select(TEZGNR40, Shape)) %>% 
    st_as_sf(., crs = target_crs)
  
  
  
  ## plot positive shapley contributions
  sp_tmap_pos_shaps <- lapply(unique(pos_shapley$name), function(x){
    
    pos_shapley_short <- pos_shapley %>% filter(name == x, rank <= 5)
    
    if(!nrow(pos_shapley_short)>1){return(NULL)}else{
    
    tm_shape(cropping_sf) +
     tm_polygons(col = 'gray95') + 
     tm_shape(shp = pos_shapley_short) + 
     tm_polygons(col = 'rank', palette = "viridis", n = 20, contrast = c(0, 0.5), border.col = 'transparent', breaks = c(1:5), labels = as.character(c(1:5))) + 
     tm_shape(ch_rivers_2) +
     tm_lines(lwd = 'CUM_LEN_LOG', scale=2, col = 'gray80', legend.lwd.show = F) + 
     tm_shape(lakes_in_ch) + 
     tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
     tm_layout(title = x)
      
    }
    }
  )
  
  # rename positive shapley values
  names(sp_tmap_pos_shaps) <- unique(pos_shapley$name)
  
  # remove those with no catchments
  missing_pos <- sapply(sp_tmap_pos_shaps, is.null, simplify = T)
  if(sum(missing_pos) != 0){sp_tmap_pos_shaps <- sp_tmap_pos_shaps[-which(missing_pos)]}
  
  # save all as pngs
  dir.create(paste0('figures/', run_type, '/', sp_name, '/', 'pos_shapley_values'), recursive = T)
  
  # loop through creation of pngs
  lapply(names(sp_tmap_pos_shaps), function(x){
    
    png(file = paste0('figures/', run_type, '/', sp_name, '/', 'pos_shapley_values', '/', x, '.png'), res = 300, height = 2000, width = 2000)
    print(sp_tmap_pos_shaps[[x]])
    dev.off()
    
  })
  
  # # save positive contributions to pdf
  # pdf(file = paste0('figures/', run_type, '/', sp_name, '/', 'pos_shapley_values.pdf'), width = 10, height = 10)
  # lapply(sp_tmap_pos_shaps, print)
  # #sp_tmap_pos_shaps[[1]]
  # dev.off()

  
  
  ## plot negative shapley contributions
  sp_tmap_neg_shaps <- lapply(unique(neg_shapley$name), function(x){
    
    neg_shapley_short <- neg_shapley %>% filter(name == x, rank <= 5)
    
    if(!nrow(neg_shapley_short)>1){return(NULL)}else{
        
    tm_shape(cropping_sf) +
      tm_polygons(col = 'gray95') + 
      tm_shape(shp = neg_shapley_short) + 
      tm_polygons(col = 'rank',palette = "viridis", n = 20, contrast = c(1, 0.5), border.col = 'transparent', breaks = c(1:5), labels = as.character(c(1:5))) + 
      tm_shape(ch_rivers_2) +
      tm_lines(lwd = 'CUM_LEN_LOG', scale=2, col = 'gray80', legend.lwd.show = F) + 
      tm_shape(lakes_in_ch) + 
      tm_polygons(border.col = 'gray80', col = 'gray95', legend.show = F) + 
      tm_layout(title = x)
      
    }
  }
  )
  
  # rename negitive shapley values
  names(sp_tmap_neg_shaps) <- unique(neg_shapley$name)
  
  # remove those with no catchments
  missing_negs <- sapply(sp_tmap_neg_shaps, is.null, simplify = T)
  if(sum(missing_negs) != 0){sp_tmap_neg_shaps <- sp_tmap_neg_shaps[-which(missing_negs)]}
  
  # save all as pngs
  dir.create(paste0('figures/', run_type, '/', sp_name, '/', 'neg_shapley_values'), recursive = T)
  
  # loop through creation of pngs
  lapply(names(sp_tmap_neg_shaps), function(x){
    
    png(file = paste0('figures/', run_type, '/', sp_name, '/', 'neg_shapley_values', '/', x, '.png'), res = 300, height = 2000, width = 2000)
    print(sp_tmap_neg_shaps[[x]])
    dev.off()
    
  })
  
}

## end of loop creating shapley values across species

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







