#### Read in final random forest models ----

# read in required packages
pacman::p_load(tidyverse, gridExtra, randomForest)

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"

# figure directory
fig_dir <- paste0("figures/", RUN, '/')
dir.create(fig_dir, recursive = T)
dir.create(paste0(fig_dir, 'perm_importance/'))

# get species
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# get directories for random forest objects
sdm_dirs <- list.files(paste0("D:/sdm-pipeline/sdm-run/", RUN_SDM), full.names = T)
sdm_dirs <- sdm_dirs[grepl(paste0(sp_list, collapse = "|"), sdm_dirs)]
all_rfs  <- paste0(sdm_dirs, '/output/final_models_rf_pa.RDS')

# get rasterized suitability scores
sp_raster_suit_pa <- paste0(sdm_dirs, "/output/raster_maps/suitability_stack_pa_rf.TIF")

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
                        full.names = T
)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")

# get data inputs across all species
all_data <- paste0(sdm_dirs, '/data/sdm_input_data.rds')




#### QUALITATITIVE SHADOW DISTRIBUTION ----

# load spatial objects
source('scripts/results/00-load-spatial.R')

# custom function to estimate quanlitative shadow distribution properties
source('scripts/functions/generate_SD_outputs.R')

# custom function to estimate quantitative shadow distribution properties 
source('scripts/functions/generate_quant_SD.R')

# natural niche factors
nn_factors <- c('ecoF_discharge_max_log10_SHAP',
                'ecoF_slope_min_log10_SHAP', 
                'stars_t_mn_m_c_SHAP', 
                'stars_t_mx_m_c_SHAP',
                'ecoF_flow_velocity_mean_SHAP', 
                'local_dis2lake_SHAP')
hab_factors <- c('local_wet_SHAP', 
                 'local_flood_SHAP', 
                 'local_imd_log10_ele_residual_SHAP', 
                 'ecoF_eco_mean_ele_residual_SHAP')
con_factors <- c('local_asym_cl_log10_SHAP')
all_factors <- c(nn_factors, hab_factors, con_factors)

# run shadow distributions across all species
all_shadow <- lapply(1:length(sp_list), function(x){
  
  generate_SD_outputs(sdm_input_data = sdm_dirs[x], 
                      raster_data = sp_raster_suit_pa[x], 
                      shap = shap_pa[x], 
                      natural_niche_factors = nn_factors,
                      habitat_factors = hab_factors,
                      conn_factors = con_factors,
                      species = sp_list[x],
                      output_folder = paste0(fig_dir, '/check_scenarios/statistics'))
  
})

names(all_shadow) <- sp_list



#### QUANTITATIVE SHADOW DISTRIBUTION ----

# output list
all_dd_method <- list()

# run loop
for(method in 1:4){
  
  # run quantitative shadow distribution
  all_dd <- lapply(1:length(all_shadow), function(x) generate_quant_SD(all_shadow[[x]]))
  
  # get names
  names(all_dd) <- sp_list
  
  # assign to output list
  all_dd_method[[method]] <- all_dd
  
}

# check structure
all_dd_method[[1]]


#### MODEL SCENARIOS IN FEATURE SPACE ----

# set example
x = 1

# read in feature space
env_data <- terra::rast('scripts/workflow example/data example/env_data.tif')


# adjust feature space by multiple "scenarios"
hum_factor <- gsub('_SHAP', '', c(hab_factors, con_factors))
list_rast <- list()
for(i in hum_factors){
  
  # improve wetland, floodplains and connectivity
  # upper clamp or lower clamp (depends on direction of threat)
  if(i %in% c("local_wet", "local_flood", "local_asym_cl_log10")){
    clamp_lwr <- as.numeric(global(env_data[[i]], fun = function(x) quantile(x, 0.99, na.rm = T)))
    clamp_upr <- as.numeric(global(env_data[[i]], fun = function(x) quantile(x, 1, na.rm = T)))
    list_rast[[which(hum_factors %in% i)]] <- clamp(env_data[[i]], lower = clamp_lwr, upper = clamp_upr)
    terra::varnames(list_rast[[which(hum_factors %in% i)]]) <- i
  }
  
  
  # upper clamp or lower clamp (depends on direction of threat)
  if(!i %in% c("local_wet", "local_flood", "local_asym_cl_log10")){
    clamp_lwr <- as.numeric(global(env_data[[i]], fun = function(x) quantile(x, 0, na.rm = T)))
    clamp_upr <- as.numeric(global(env_data[[i]], fun = function(x) quantile(x, 0.01, na.rm = T)))
    list_rast[[which(hum_factors %in% i)]] <- clamp(env_data[[i]], lower = clamp_lwr, upper = clamp_upr)
    terra::varnames(list_rast[[which(hum_factors %in% i)]]) <- i
  }
  
}

# get scenario rasters
scen_rast <- rast(list_rast)

# adjust layers
scen_data2 <- env_data
scen_data2[[hum_factor]] <- scen_rast

tm_shape(scen_data2) + 
  tm_raster(style = 'cont', n = 5, midpoint = NA) + 
  tm_facets(free.scales = T)



# run through all species making predictions over adjusted feature space
for(sp in sp_list){

  # read in random forests
  rf1 <- readRDS(all_rfs[[grep(sp,all_rfs)]])
  
  # get predictions based on adjustment scenarios
  exp_pred <- terra::predict(scen_data2, rf1, type = "prob")[[2]]
  obs_pred <- terra::predict(env_data,   rf1, type = "prob")[[2]]
  
  shadow <- (obs_pred / exp_pred)
  
  # aggregate to sub-catchment values
  ext_shadow <- terra::extract(shadow, 
                               vect(subcatchments_final %>% select(TEILEZGNR)), 
                               fun = function(x) {mean(x, na.rm = T)})
  
  ext_shadow2 <- cbind(subcatchments_final %>% select(TEILEZGNR), ext_shadow)
  ext_shadow2 <- ext_shadow2 %>% select(-ID) %>% rename(SCEN_SD_OratioE = X1)
  
  # join to subcatchment object
  all_scen_shadow_sp <- lapply(1:length(all_dd_method), function(x){
    z <- left_join(all_dd_method[[x]][[sp]], ext_shadow2 %>% st_drop_geometry()) %>% 
      mutate(method = paste0('method', x), 
             species_name = sp) %>% 
      st_drop_geometry()
    return(z)
  })
  
  # create output directory and save files for creating plots all together
  dir.create(paste0(fig_dir, 'check_scenarios/shadow_scenario/'), recursive = T)
  saveRDS(all_scen_shadow_sp, file = paste0(fig_dir, 'check_scenarios/shadow_scenario/', sp, '.RDS'))
  
}


shad_files <- list.files(paste0(fig_dir, 'check_scenarios/shadow_scenario/'), full.names = T)

shad_data <- lapply(shad_files, readRDS)

shad_data <- bind_rows(shad_data)

# add in R2
r2 <- shad_data %>%
  filter(SD_OratioE < 1, 
         method != 'method4') %>% 
  rename(method_shad = method) %>% 
  group_by(species_name, method_shad) %>% 
  do(cor = tidy(cor.test(.$SCEN_SD_OratioE, .$SD_OratioE, method = 'pearson', na.rm = T))) %>% 
  unnest(cor) %>% 
  select(-method) %>% 
  rename(method = method_shad)

r2 %>% group_by(method) %>% do(test = summary(.$estimate)) %>% pull(test)
r2 %>% group_by(method) %>% do(test = sd(.$estimate)) %>% pull(test)

# plot correspondence between approach based on SHAP adjustment and scenarios in the feature space
png(filename = paste0(fig_dir, 'check_scenarios/shadow_scenario.png'), res = 300, height = 3000, width = 2000)
ggplot(data = shad_data %>% 
         filter(method != 'method4')) + 
  # add in r2 values
  geom_text(data = r2, aes(x = 0.2, y = 0.9, label = round(estimate, digits = 2))) + 
  # add in hexagons
  geom_hex(aes(x = SD_OratioE, y = SCEN_SD_OratioE, fill = stat(log10(count)), alpha = stat(log10(count))), bins = 20) + 
  #geom_point(aes(x = SD_OratioE, y = SCEN_SD_OratioE)) + 
  viridis::scale_fill_viridis(option = 'inferno', 
                              begin = 1, end = 0,
                              na.value = viridis::viridis(10,option = 1, end = 0.9)[10], 
                              name = 'log10(no. obs.)') +
  scale_alpha(range = c(0.1, 1)) + 
  stat_smooth(aes(x = SD_OratioE, y = SCEN_SD_OratioE), method = 'loess', col = 'black', size = 0.5) + 
  stat_smooth(aes(x = SD_OratioE, y = SCEN_SD_OratioE), method = 'lm', lty = 2, col = 'black', size = 0.5) + 
  geom_abline(col = 'black') + 
  facet_grid(species_name ~ method) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.text.y = element_text(size = 6)) + 
  ylim(c(0,1)) + 
  xlim(c(0,1)) + 
  xlab('SHAP shadow distribution') + 
  ylab('Feature scenario shadow distribution')
dev.off()

# check correspondence
ggplot(data = test_shad) + 
  geom_point(aes(x = SD_OratioE, y = SCEN_SD_OratioE)) + 
  geom_abline() + 
  stat_smooth(aes(x = SD_OratioE, y = SCEN_SD_OratioE), method = 'lm') + 
  theme_bw() + 
  theme(aspect.ratio = 1)

test_shad <- test_shad %>% filter(SD_OratioE != 0)

stats::cor.test(x = test_shad$SCEN_SD_OratioE, y = test_shad$SD_OratioE, na.rm = T, method = 'pearson')
stats::cor.test(x = test_shad$SCEN_SD_OratioE, y = test_shad$SD_OratioE, na.rm = T, method = 'spearman')

tm_shape(test_shad %>% filter(!is.na(SD_OratioE), SCEN_SD_OratioE <= 1)) +
  tm_fill(col = 'SCEN_SD_OratioE', 
          style = 'cont',
          palette = 'Spectral',
          title = 'Quantitative shadow distribution', 
          colorNA = 'transparent', 
          #breaks = c(0, 0.5, 1), 
          legend.is.portrait = F) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.5) + 
  tm_layout(frame = F,
            bg.color = "transparent")


#### ----

data1 <- readRDS(all_data[[x]])
shap1 <- readRDS(shap_pa[[x]])
train1 <- data1@pa_data$full_data
occ <- as.numeric(train1[,1])
train1 <- train1[,-c(1:3)]





