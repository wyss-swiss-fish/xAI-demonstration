### Script to calculate areas within the natural niche for all species and then quantify the areas of dark distributions and dark diversity ----


#### 1. Load packages----

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)


#### 2. Set directories for loading and saving objects and load spatial ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_MARCH_v6/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_MARCH_v6"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS('figures/ubelix_SDM_RF_MARCH_v6/evaluations/subset_sp.RDS')

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
                        full.names = T
)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]

# create shapley folders per species: subcatchments
shap_po <- paste0(shap_dirs, "/shapley_rf_po.RDS")
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")

# create shapley folders per species: rasters
shap_po_rast <- paste0(shap_dirs, "/shapley_rf_po.TIF")
shap_pa_rast <- paste0(shap_dirs, "/shapley_rf_pa.TIF")


# get directories for response curve objects
sdm_dirs <- list.files(paste0("D:/sdm-pipeline/sdm-run/", RUN), full.names = T)
sdm_dirs <- sdm_dirs[grepl(paste0(sp_list, collapse = "|"), sdm_dirs)]

# response curve paths
rc_po <- paste0(sdm_dirs, "/output/response_curves/response_curves_rf_po.rds")
rc_pa <- paste0(sdm_dirs, "/output/response_curves/response_curves_rf_pa.rds")

# range maps
sp_raster_pres_po <- paste0(sdm_dirs, "/output/raster_maps/presence_stack_po_rf.TIF")
sp_raster_pres_pa <- paste0(sdm_dirs, "/output/raster_maps/presence_stack_pa_rf.TIF")

# suitability maps
sp_raster_suit_po <- paste0(sdm_dirs, "/output/raster_maps/suitability_stack_po_rf.TIF")
sp_raster_suit_pa <- paste0(sdm_dirs, "/output/raster_maps/suitability_stack_pa_rf.TIF")

# read in environmental data
env_data <- rast(paste0(dd_env, '/ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif'))

# load spatial objects
source('scripts/results/00-load-spatial.R')


#### 3. Set variable names ----

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


#### 4. Get raster maps of shapley values ----

# natural niche factors
nn_factors <- c('ecoF_discharge_max_log10_SHAP','ecoF_slope_min_log10', 'stars_t_mn_m_c_SHAP', 
                'stars_t_mx_m_c_SHAP','ecoF_flow_velocity_mean_SHAP', 'local_dis2lake')
hab_factors <- c('local_flood_SHAP', 'local_imd_log10_ele_residual_SHAP', 'ecoF_eco_mean_ele_residual_SHAP')
con_factors <- c('local_asym_cl_log10_SHAP')
all_factors <- c(nn_factors, hab_factors, con_factors)

all_dd <- lapply(1:length(sp_list), function(i){
  
   ## loop for i 
   do_pa <- file.exists(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
   if(do_pa != T){next()}
   
   ## ORGANISE SHAPELY VALUES
   # read in shapley values for a given species
   shap_rast_pa_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
   # take only the focal natural factors
   shap_rast_pa_i <- shap_rast_pa_i[[which(names(shap_rast_pa_i) %in% all_factors)]]
   # define within niche area
   nn_shap_rast <- shap_rast_pa_i[[which(names(shap_rast_pa_i) %in% nn_factors)]]
   # remove non_na
   shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]
   # order by name
   shap_rast_pa_i <- shap_rast_pa_i[[sort(names(shap_rast_pa_i))]]
   
   # ORGANISE SPECIES PRESENCE OR ABSENCE
   model_data   <- readRDS(paste0(sdm_dirs[i], '/data/sdm_input_data.rds'))@pa_data$full_data
   sp_rast <- rast(sp_raster_suit_pa[i])
   pred_occ <- terra::extract(sp_rast, model_data[c('X', 'Y')])
   threshold <- ecospat::ecospat.max.tss(pred_occ[,2], model_data$occ)$max.threshold
   presence_map <- sp_rast
   absence_map  <- sp_rast
   # define presence map
   presence_map[sp_rast<threshold] <- 0
   presence_map[sp_rast>threshold] <- 1
   # define absence map
   absence_map[sp_rast<threshold] <- 1
   absence_map[sp_rast>threshold] <- 0
   
   # align maps
   shap_rast_pa_i <- resample(shap_rast_pa_i, sp_rast)
   absence_map[is.na(shap_rast_pa_i[[1]])]  <- NA
   presence_map[is.na(shap_rast_pa_i[[1]])] <- NA
   
   # DEFINE NICHES, EXPECTED DISTRIBUTIONS AND DARK DISTRIBUTIONS
   # define within niche area by summation
   nn_pos_sum <- app(nn_shap_rast, sum, na.rm = T)
   nn_pos_sum <- resample(nn_pos_sum, absence_map)
   
   # get areas that are absent that could be present
   suit_nn_pos_sum <- sp_rast + nn_pos_sum
   suit_nn_pos_sum[suit_nn_pos_sum < threshold] <- 0
   suit_nn_pos_sum[suit_nn_pos_sum > threshold] <- 1
   
   
   # define within niche area by mean
   nn_pos_mean <- app(nn_shap_rast, mean, na.rm = T)
   nn_pos_mean <- resample(nn_pos_mean, absence_map)
   
   # make expected map based on summation
   expected_dist_sum <- nn_pos_sum
   expected_dist_sum[expected_dist_sum[]<0] <- 0
   expected_dist_sum[expected_dist_sum[]>0] <- 1
   
   # make expected map based on mean
   expected_dist_mean <- nn_pos_mean
   expected_dist_mean[expected_dist_mean[]<0] <- 0
   expected_dist_mean[expected_dist_mean[]>0] <- 1
   
   # make expected but absent map based on summation
   dark_dist_sum <- absence_map
   dark_dist_sum[nn_pos_sum[]<0] <- 0 # exclude areas from absences where natural factors do not support
   
   # make expected but absent map based on summation
   dark_dist_mean <- absence_map
   dark_dist_mean[nn_pos_mean[]<0] <- 0 # exclude areas from absences where natural factors do not support
   
   # make names cleaner
   names(presence_map)   <- sp_list[i]
   names(expected_dist_sum)  <- sp_list[i]
   names(expected_dist_mean) <- sp_list[i]
   names(dark_dist_sum)  <- sp_list[i]
   names(dark_dist_mean) <- sp_list[i]
   
   return(list(observed_dist       = presence_map, 
               expected_dist_sum   = expected_dist_sum,
               expected_dist_mean  = expected_dist_mean, 
               dark_dist_sum       = dark_dist_sum, 
               dark_dist_mean      = dark_dist_mean))
})

# get names
names(all_dd) <- sp_list

## REMOVE NON-NATIVE SPECIES
all_dd <- all_dd[-which(names(all_dd) == 'Oncorhynchus mykiss')]

## SUM RESULTING DISTRIBUTIONS OVER ALL SPECIES
# compile observed distributions
observed_dist <- app(rast(lapply(all_dd, function(x) x$observed_dist)), sum, na.rm = T)
# compile expected distributions summing shapley values of natural niche
sum_expected_dist<- app(rast(lapply(all_dd, function(x) x$expected_dist_sum)), sum, na.rm = T)
# compile expected distributions mean of shapley values of natural niche
mean_expected_dist<- app(rast(lapply(all_dd, function(x) x$expected_dist_mean)), sum, na.rm = T)
# compile dark distributions based on summation of shapley values of natural niche
sum_dark_dist <- app(rast(lapply(all_dd, function(x) x$dark_dist_sum)), sum, na.rm = T)
# compile dark distributions based on mean of shapley values of natural niche
mean_dark_dist <- app(rast(lapply(all_dd, function(x) x$dark_dist_mean)), sum, na.rm = T)

# look at rough plots: OBSERVED DISTRIBUTION
plot(observed_dist)
# look at rough plots: OBSERVED DISTRIBUTION
plot(sum_expected_dist)
plot(sum_dark_dist)
# look at rough plots: EXPECTED AND DARK FOR MEAN SUMMATION OF SHAPLEY VALUES
plot(mean_expected_dist)
plot(mean_dark_dist)
# look at rough plots of community completeness
completeness <- observed_dist / sum_expected_dist
completeness[completeness[] > 1] = 1
plot(1/completeness)
# look at plots of the ratio of dark to observed diversity
obs_unobs_ratio <- observed_dist / sum_dark_dist
plot(obs_unobs_ratio)
plot(1/obs_unobs_ratio)

plot(completeness - obs_unobs_ratio)
plot(sum_expected_dist - sum_dark_dist)

# the approach for the sum and the mean give identical results. 
plot(sum_expected_dist-mean_expected_dist)

# testing
plot(rast_dd_sum)
hist(rast_dd_sum[rast_dd_sum>0])
