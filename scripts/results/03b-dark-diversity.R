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

# get run to make figures for
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
   
   # get areas that are absent which could be present based on increasing suitability
   expected_dist_sum <- sp_rast + nn_pos_sum
   expected_dist_sum[expected_dist_sum < threshold] <- 0
   expected_dist_sum[expected_dist_sum > threshold] <- 1
   
   # make expected but absent map based on niche effect that are predict presence
   dark_dist_sum <- absence_map
   dark_dist_sum[expected_dist_sum[]==0] <- 0 # exclude areas from absences where natural factors do not support

   # make names cleaner
   names(presence_map)   <- sp_list[i]
   names(expected_dist_sum)  <- sp_list[i]
   names(dark_dist_sum)  <- sp_list[i]

   return(list(observed_dist       = presence_map, 
               expected_dist_sum   = expected_dist_sum,
               dark_dist_sum       = dark_dist_sum))
})

# get names
names(all_dd) <- sp_list

## REMOVE NON-NATIVE SPECIES
all_dd <- all_dd[-which(names(all_dd) == 'Oncorhynchus mykiss')]

## SUM RESULTING DISTRIBUTIONS OVER ALL SPECIES
# compile observed distributions
sum_observed_dist <- app(rast(lapply(all_dd, function(x) x$observed_dist)), sum, na.rm = T)
# compile expected distributions summing shapley values of natural niche
sum_expected_dist<- app(rast(lapply(all_dd, function(x) x$expected_dist_sum)), sum, na.rm = T)
# compile dark distributions based on summation of shapley values of natural niche
sum_dark_dist <- app(rast(lapply(all_dd, function(x) x$dark_dist_sum)), sum, na.rm = T)

# look at rough plots: OBSERVED DISTRIBUTION
plot(sum_observed_dist)
# look at rough plots: OBSERVED DISTRIBUTION
plot(sum_expected_dist)
plot(sum_dark_dist)
# look at rough plots of community completeness
completeness <- sum_observed_dist / sum_expected_dist
completeness[completeness[] > 1] = 1
# clip to env_data
completeness <- mask(completeness, vect(cropping_sf))
plot(completeness)

# look at plots of the ratio of dark to observed diversity
obs_unobs_ratio <- sum_observed_dist - sum_dark_dist
obs_unobs_ratio[obs_unobs_ratio[] > 0] = 0
plot(obs_unobs_ratio)

## make baseplot for rivers
river_base <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F) + 
  tm_layout(frame = F,
            bg.color = "transparent")


## make tmap plots
observed_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(sum_observed_dist, raster.downsample = F) +
  tm_raster(style = 'cont', legend.is.portrait = F, n = 9) +
  river_base
observed_diversity_plot

expected_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(sum_expected_dist, raster.downsample = F) +
  tm_raster(style = 'cont', legend.is.portrait = F, n = 9) +
  river_base
expected_diversity_plot

dark_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(sum_dark_dist, raster.downsample = F) +
  tm_raster(style = 'cont', legend.is.portrait = F, n = 9) +
  river_base
dark_diversity_plot

completeness_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(completeness, raster.downsample = F) +
  tm_raster(legend.is.portrait = F, 
            palette = "Spectral", n = 100, contrast = c(0, 1)) +
  river_base
completeness_plot 

# get legends
observed_diversity_plot_legend <- observed_diversity_plot + tm_layout(legend.only = T)
completeness_plot_legend       <- completeness_plot + tm_layout(legend.only = T)

dir.create(paste0(fig_dir, '/dark_diversity'), recursive = T)
pdf(paste0(fig_dir, '/dark_diversity/combined_maps.pdf'), width = 16, height = 8, 
    bg = 'transparent')
tmap_arrange(observed_diversity_plot + tm_layout(legend.show = F), 
             expected_diversity_plot + tm_layout(legend.show = F), 
             dark_diversity_plot + tm_layout(legend.show = F), 
             completeness_plot + tm_layout(legend.show = F),
             nrow = 2, ncol = 2)
dev.off()

# save pdf of legends
pdf(paste0(fig_dir, '/dark_diversity/diversity_legend.pdf'), width = 2, height = 2, 
    bg = 'transparent')
observed_diversity_plot_legend
dev.off()

# save pdf of legends
pdf(paste0(fig_dir, '/dark_diversity/completeness_legend.pdf'), width = 2, height = 2, 
    bg = 'transparent')
completeness_plot_legend
dev.off()

hist(completeness, breaks = 5)

quantile(completeness[], na.rm  = T)
table(completeness[] < 0.5)

#### 5. Extract values of diversity across subcatchment units and make summaries ----

## here take the diversity estimates and summarise, but remove all areas where we do not expect species to occur anyway (and also NA areas in our models)
subcatchments_dark_div <- subcatchments_final

# extractions for expected diversity
catchment_expected_div <- extract(sum_expected_dist, 
                                  subcatchments_dark_div %>% select(TEILEZGNR), 
                                  fun = mean, 
                                  na.rm = T)


# extractions for observed diversity
catchment_observed_div <- extract(sum_observed_dist, 
                                  subcatchments_dark_div %>% select(TEILEZGNR), 
                                  fun = mean, 
                                  na.rm = T)


# extractions for dark diversity
catchment_dark_div <- extract(sum_dark_dist, 
                              subcatchments_dark_div %>% select(TEILEZGNR), 
                              fun = mean, 
                              na.rm = T)

# extractions for completeness
catchment_completeness <- extract(completeness, 
                                  subcatchments_dark_div %>% select(TEILEZGNR), 
                              fun = mean, 
                              na.rm = T)

# assign as new columns
subcatchments_dark_div$catchment_expected_div <- catchment_expected_div$sum
subcatchments_dark_div$catchment_observed_div <- catchment_observed_div$sum
subcatchments_dark_div$catchment_dark_div <- catchment_dark_div$sum
subcatchments_dark_div$catchment_completeness <-  catchment_completeness$sum

# remove all catchments where expected diversity is 0
subcatchments_dark_div <- subcatchments_dark_div %>% filter(catchment_expected_div != 0, !is.na(catchment_expected_div))

# get the summaries
summary(subcatchments_dark_div$catchment_observed_div)
summary(subcatchments_dark_div$catchment_expected_div)
summary(subcatchments_dark_div$catchment_dark_div)
summary(subcatchments_dark_div$catchment_completeness)

# get histograms
hist(subcatchments_dark_div$catchment_observed_div)
hist(subcatchments_dark_div$catchment_expected_div)
hist(subcatchments_dark_div$catchment_dark_div)
hist(subcatchments_dark_div$catchment_completeness)

# histograms of observed diversity
pdf(paste0(fig_dir, '/dark_diversity/observed_div_hist.pdf'), width = 2, height = 2)
ggplot(data = subcatchments_dark_div) + 
  geom_histogram(aes(x = catchment_observed_div), bins = 8, fill = '#FE9F2F') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 1.3) + 
  xlab('no. species') + 
  scale_x_continuous(breaks = seq(0,8,2))
dev.off()

# histogram of expected diversity
pdf(paste0(fig_dir, '/dark_diversity/expected_div_hist.pdf'), width = 2, height = 2)
ggplot(data = subcatchments_dark_div) + 
  geom_histogram(aes(x = catchment_expected_div), bins = 8, fill = '#FE9F2F') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 1.3) + 
  xlab('no. species') + 
  scale_x_continuous(breaks = seq(0,8,2))
dev.off()

# histogram of dark diversity
pdf(paste0(fig_dir, '/dark_diversity/dark_div_hist.pdf'), width = 2, height = 2)
ggplot(data = subcatchments_dark_div) + 
  geom_histogram(aes(x = catchment_dark_div), bins = 8, fill = '#FE9F2F') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 1.3) + 
  xlab('no. species') + 
  scale_x_continuous(breaks = seq(0,8,2))
dev.off()


# histogram of completeness
subcatchments_dark_div$catchment_completeness_fact <- round(subcatchments_dark_div$catchment_completeness, 1)*100
hist_this <- subcatchments_dark_div %>% st_drop_geometry() %>% group_by(catchment_completeness_fact) %>% count()
pdf(paste0(fig_dir, '/dark_diversity/completeness_hist.pdf'), width = 2, height = 2)
ggplot(data = hist_this) + 
  geom_bar(aes(x = catchment_completeness_fact, y = n, fill = catchment_completeness_fact), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 1.3, 
        legend.position = 'none') + 
  xlab('completeness') + 
  ylab('count') + 
  scale_x_continuous(breaks = seq(0,100,50), 
                     labels = paste0(seq(0,100,50), '%')) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral"))
dev.off()

