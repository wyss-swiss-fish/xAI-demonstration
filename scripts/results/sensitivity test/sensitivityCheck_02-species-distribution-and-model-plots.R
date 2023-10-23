#### Aim of script is to generate series of figures summarise SMDs and shapley values for each species ----

# This script generates a series of figures for exploration including:
# - maps of species shapley values for presence-only and presence-absence SDMs (although we further consider only PA models)
# - response curves of shapley values for presence-only and presence-absence SDMs (although we further consider only PA models)

#### 1. Load in packages data ----

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

#### 2. Set directories for loading and saving objects ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# read in environmental data
env_data <- rast(paste0(dd_env, '/ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif'))

# load in all spatial objects for making maps and plots
source('scripts/results/00-load-spatial.R')


### Set up variable names
# set vector of focal variables for assigning new names, this must be consistent with the models fitted
vars <- c('ecoF_discharge_max_log10', 
          'ecoF_slope_min_log10', 
          'ecoF_flow_velocity_mean', 
          'stars_t_mx_m_c', 
          'stars_t_mn_m_c',
          'local_asym_cl_log10',
          'local_dis2lake',
          'ecoF_eco_mean_ele_residual', 
          'local_imd_log10_ele_residual',
          'local_wet',
          'local_flood')

# set vars
vars_shap <- paste0(vars, "_SHAP")

# rename variables
vars_renamed = c('discharge', 
                 'slope', 
                 'flow velocity', 
                 'temperature max', 
                 'temperature min', 
                 'connectivity', 
                 'distance to lake', 
                 'ecomorphology', 
                 'urbanisation',
                 'wetland', 
                 'floodplains')

vars_renamed <- cbind(vars_shap, vars_renamed)


#### For each run save response curves across sensitivity tests and investigate variation in results ----

MODEL_list <- c('ubelix_SDM_RF_APRIL_V1',
              'ubelix_SDM_RF_JULY_inHouse_V2', 
              'ubelix_SDM_RF_JULY_Kanton_subset_V1',
              'ubelix_SDM_RF_JULY_Kanton_V1')

SHAP_list <-  c('ubelix_SDM_RF_APRIL_V1_02',
                'ubelix_SDM_RF_JULY_inHouse_V2_01', 
                'ubelix_SDM_RF_JULY_Kanton_subset_V1_01',
                'ubelix_SDM_RF_JULY_Kanton_V1_01')

# output response curves
all_rc <- list()

for(run in 1:length(MODEL_list)){

# figure directory
fig_dir <- paste0("figures/", MODEL_list[run], '/')

# get run to mak figures for
RUN <- SHAP_list[run]

# get species of interest
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
  full.names = T
)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]

# create shapley folders per species
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")


#### 6. Make shapley based response curves for interpretations ----

dir.create(paste0(fig_dir, "all_sp_plots/"), recursive = T)

# make plots of response curves for presence-absence data
pdf(paste0(fig_dir, "all_sp_plots/shap_rc_all_pa", ".pdf"), width = 14, height = 14)
for (i in 1:length(sp_list)) {
  
    # read in raw shapleys
    shap_i <- readRDS(shap_pa[i])
    # join together shapley values to environmental data for each TEILEZGNR
    shap_i <- left_join(shap_i %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
    
    # get the shapley values
    grouped_rc_shap <- shap_i %>% 
      pivot_longer(cols = matches(vars)) %>% 
      select(name, value) %>% 
      mutate(shap = ifelse(grepl('_SHAP', name), 'shapley', 'env'),
             name = gsub('_SHAP', '', name)) %>% 
      pivot_wider(names_from = shap, values_from = value) %>% 
      unnest(c(shapley, env)) %>% 
      filter(!is.na(shapley)) %>% 
      merge(., vars_renamed %>% 
              data.frame %>% 
              mutate(vars_shap = gsub('_SHAP', '', vars_shap)), 
            by.x = 'name', 
            by.y = 'vars_shap') %>% 
      data.frame
    
    # get the minimum and maximum
    min_y <- min(grouped_rc_shap$shapley, na.rm = T)
    max_y <- max(grouped_rc_shap$shapley, na.rm = T)
    
    print(head(grouped_rc_shap))
    
    # make shapley plot for all variables for a single species
    rc_shap_po_plot <- ggplot(data = grouped_rc_shap %>% group_by(name) %>% sample_n(1000), aes(x=env, y = shapley)) + 
      facet_wrap(~vars_renamed, scales = 'free') +
      geom_point(col = 'gray50') +
      stat_smooth(col = 'black', se = F) +
      geom_hline(aes(yintercept = 0)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            aspect.ratio = 0.5) +
      scale_y_continuous(limits = c(min_y, max_y)) +
      ggtitle(sp_list[i])
    
    print(rc_shap_po_plot)
    
  }
dev.off()



## get together all hapes per subset
all_shap_env <- lapply(1:length(sp_list), function(x){ 
  
  if (file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    # join together shapley values to environmental data for each TEILEZGNR
    shap_x <- left_join(shap_x %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
    
    # get the shapley values
    grouped_rc_x <- shap_x %>% 
      pivot_longer(cols = matches(vars)) %>% 
      select(name, value) %>% 
      mutate(shap = ifelse(grepl('_SHAP', name), 'shapley', 'env'),
             name = gsub('_SHAP', '', name)) %>% 
      pivot_wider(names_from = shap, values_from = value) %>% 
      unnest(c(shapley, env)) %>% 
      filter(!is.na(shapley)) %>% 
      merge(., vars_renamed %>% 
              data.frame %>% 
              mutate(vars_shap = gsub('_SHAP', '', vars_shap)), 
            by.x = 'name', 
            by.y = 'vars_shap') %>% 
      data.frame
    
    grouped_rc_x$species_name <- sp_list[x]
    
    print(x)
    return(grouped_rc_x)
  }else{NULL}
}
)

# bind all together
all_shap_env_pa <- bind_rows(all_shap_env)

# estimate variable importance (mean(abs(shapley))) for species and variable
all_shap_env_pa_simple <- all_shap_env_pa %>% 
  group_by(vars_renamed, species_name) %>% 
  do(mean_var_imp = mean(abs(.$shapley), na.rm=T), 
     range_var_imp = sd(.$shapley), na.rm = T) %>% 
  unnest(c(mean_var_imp, range_var_imp))

# plot the response curves as a large matrix across all species
pdf(paste0(fig_dir, "all_sp_plots/allsp_rc", ".pdf"), width = 13, height = 13)
ggplot(data = all_shap_env_pa %>% 
         left_join(., all_shap_env_pa_simple) %>% 
         group_by(vars_renamed, species_name) %>% 
         sample_n(5000)) +
  geom_point(aes(x = env, y = shapley, col = mean_var_imp)) +
  stat_smooth(aes(x = env, y = shapley), col = 'black', se = F) +
  facet_grid(species_name ~ vars_renamed, scales = 'free_x') +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        legend.position = 'bottom') +
  scale_colour_viridis_c(name = 'variable importance') + 
  xlab('environmental variable') + 
  ylab('Shapley value')
dev.off()


## output all shapley across all variables and all data
all_rc[[run]] <- all_shap_env_pa %>% 
  left_join(., all_shap_env_pa_simple) %>% 
  group_by(vars_renamed, species_name) %>% 
  mutate(data_subset = MODEL_list[run])

} # end of run loop 


# bind all the rows together and combine into one dataset and subset
all_rc_sub <- bind_rows(all_rc) %>% 
  group_by(vars_renamed, species_name, data_subset) %>% 
  sample_n(5000)

# response curve sensitivity checks
dir.create('figures/sensitivity-check/', recursive = T)
pdf(paste0("figures/sensitivity-check/", "allsp_rc", ".pdf"), width = 13, height = 13)
ggplot(data = all_rc_sub) +
  stat_smooth(aes(x = env, y = shapley, col = data_subset), se = F) +
  facet_grid(species_name ~ vars_renamed, scales = 'free_x') +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        legend.position = 'bottom') +
  xlab('environmental variable') + 
  ylab('Shapley value')
dev.off()




#### for each of the datasets have a look at the distribution of the environmental data that they cover ----

# set the vocal variables
DATA_list <- c('fish-presenceAbsence_2010.csv',
               'fish-presenceAbsence_2010_inHouse.csv', 
               'fish-presenceAbsence_2010_Kanton.csv', 
               'fish-presenceAbsence_2010_Kanton_subset.csv')


values <- list()

for(subset in 1:length(DATA_list)){
  
  # read in raw records
sdm_data <- read_csv(paste0(dd, 'sdm-pipeline/species-records-final/', DATA_list[subset]))

# subset to unique records
sdm_xy <- sdm_data %>% filter(year > 2010) %>% select(X, Y) %>% unique() %>% 
  st_as_sf(., coords = c('X', 'Y'), crs = 'wgs84') %>% 
  st_transform(., st_crs(env_data))

# extract the environmental values of each dataset
extr_xy <- extract(env_data, sdm_xy)

# process extracted data
values[[subset]] <- extr_xy %>% 
  select(any_of(vars)) %>% 
  pivot_longer(cols = vars) %>% 
  mutate(data_subset = gsub('.csv', '', DATA_list[subset]))

}
values_combined <- bind_rows(values)

# revalue names
values_combined$data_subset <- recode(values_combined$data_subset, 
                                      'fish-presenceAbsence_2010' = 'all data', 
                                      'fish-presenceAbsence_2010_inHouse' = 'in house', 
                                      'fish-presenceAbsence_2010_Kanton' = 'kanton only', 
                                      'fish-presenceAbsence_2010_Kanton_subset' = 'remove alpine')

# clean variable names
values_combined <- left_join(values_combined, data.frame(vars_renamed) %>% mutate(vars_shap = gsub('_SHAP', '', .$vars_shap)), by = c('name' = 'vars_shap'))


# across all variables view environmental distributions of the data

png('figures/env_data_across_subsets.png', res = 300, width = 3000, height = 2000)
ggplot(data = values_combined) + 
  geom_boxplot(aes(x = data_subset, y = value, col = data_subset)) +
  theme_bw() + 
  facet_wrap(~vars_renamed, scales = 'free') + 
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(angle = 0, hjust = 0, vjust = 1)) + 
  labs(color = " ")
dev.off()
