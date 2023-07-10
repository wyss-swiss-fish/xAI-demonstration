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

# load actual random forests
rf_pa <- paste0(sdm_dirs, "/output/final_models_rf_pa.RDS")

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


#### 4. Shadow distribution summaries for all species ----

# custom function to estimate shadow distribution properties and output as a shapefile
source('scripts/functions/shadow_distribution.R')

# natural niche factors
nn_factors <- c('ecoF_discharge_max_log10_SHAP','ecoF_slope_min_log10_SHAP', 'stars_t_mn_m_c_SHAP', 
                'stars_t_mx_m_c_SHAP','ecoF_flow_velocity_mean_SHAP', 'local_dis2lake_SHAP')
hab_factors <- c('local_wet_SHAP', 'local_flood_SHAP', 'local_imd_log10_ele_residual_SHAP', 'ecoF_eco_mean_ele_residual_SHAP')
con_factors <- c('local_asym_cl_log10_SHAP')
all_factors <- c(nn_factors, hab_factors, con_factors)

# run shadow distributions across all species
all_shadow <- lapply(1:length(sp_list), function(x){
  
  shadow_distribution(sdm_input_data = sdm_dirs[x], 
                      raster_data = sp_raster_suit_pa[x], 
                      shap = shap_pa[x], 
                      natural_niche_factors = nn_factors,
                      habitat_factors = hab_factors,
                      conn_factors = con_factors,
                      species = sp_list[x],
                      output_folder = paste0(fig_dir, '/shadow_dist_summaries/statistics'))
  
})

names(all_shadow) <- sp_list

# summarise shadow distribution summaries across all species
summary_shadow <- bind_rows(lapply(list.files(paste0(fig_dir, '/shadow_dist_summaries/statistics'), full.names = T), read_csv))

# group and aggregate
summary_shadow_mean <- summary_shadow %>% 
  group_by(property) %>% 
  do(mean = round(mean(.$value, na.rm = T), 2),
     sd   = round(sd(.$value, na.rm = T),2)) %>% 
  unnest(cols = c(mean, sd))

write_csv(summary_shadow_mean, file = paste0(fig_dir, '/shadow_dist_summaries/summary_across_species.csv'))

#### 5. Biplot of threatened sub-cathcments for each species ----

## biplots 
dir.create(paste0(fig_dir, '/shadow_dist_summaries/biplot/'), recursive = T)
lapply(all_shadow, function(x){
 
   # create scatter plot of shapley values against suitability classified by threat
   pdf(paste0(fig_dir, '/shadow_dist_summaries/biplot/', unique(x$species_name), '_biplot.pdf'), 
       width = 10, height = 5, bg = 'transparent')
  
   print(ggplot() + 
     geom_point(data = x %>% arrange(niche_categories), 
                aes(y = suitability, 
                    x = natural_niche_value, 
                    col = niche_categories, 
                    size = niche_categories,
                    alpha = is.na(presence)), 
                stroke = 0) + 
     xlab(expression('abiotic niche' ~phi~ 'values')) + 
     ylab('predicted habitat suitability') + 
     # scale_colour_manual('', values = c('#E8E8E8','#06C4C4','#D1D111', '#D16111', '#F70000')) + 
     theme_bw() + 
     theme(panel.grid = element_blank(),
           legend.position = 'left', 
           panel.border = element_blank(),
           panel.background = element_blank(),
           plot.background = element_blank(),
           axis.line = element_line()) + 
     geom_vline(aes(xintercept = 0), lty = 2) + 
     scale_size_manual(values = c(1,2,2,2,2)) + 
     scale_alpha_manual(values = c(1,0.4)))

   dev.off()
   
})


## threat maps within ecological niche
dir.create(paste0(fig_dir, '/shadow_dist_summaries/threat_niche/'), recursive = T)
levels = c('1. outside ecological niche',
           '2. inside ecological niche',
           '3. poor connectivity',
           '4. poor habitat',
           '5. poor connectivity and habitat')
# map the threat maps
lapply(all_shadow, function(x) {
  # create scatter plot of shapley values against suitability classified by threat
  pdf(paste0(fig_dir, "/shadow_dist_summaries/threat_niche/", unique(x$species_name), "_threat_niche.pdf"),
    width = 10, height = 5, bg = "transparent"
  )
  print(tm_shape(x) +
    tm_fill(
      col = "niche_categories",
      style = "fixed",
      palette = c("#E8E8E8", "#06C4C4", "#D1D111", "#D16111", "#F70000"),
      breaks = levels,
      title = ""
    ) +
    tm_shape(river_intersect_lakes) +
    tm_lines(legend.show = F, col = "black") +
    tm_shape(lakes) +
    tm_polygons(border.col = "black", col = "white", legend.show = F) +
    tm_layout(frame = F, bg.color = "transparent", title = unique(x$species_name)))
  dev.off()
})


#### 6. Create the maps of shadow distributions from spatial objects shapley values ----

# this creates additional columns on our spatial objects that summarise the expected, observed and shadow distributions

all_dd <- lapply(1:length(all_shadow), function(i){
  
  ####
  ## READ IN DATA
  # get species i 
  sp <- all_shadow[[i]]$species_name[1]
  print(sp)
  
  # get shapley shape files distribution
  shap_sp_i <- all_shadow[[i]]
  
  
  ####
  ## DEFINE THE SUITABILITY PREDICTIONS USING THE BASELINE PREDICTION AND THE SUMMED SHAPLEY VALUES
  
  # get the baseline value as the mean of all habitat suitability predictions
  baseline_value <- mean(shap_sp_i$suitability, na.rm = T)
    
  # sum all the shapley values to get the net prediction as a deviaiton from the baseline mean prediction
  shap_sp_i$shap_all_sum <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% all_factors]), na.rm = T)
  
  # get the baseline suitability values estimated as the sum of the shapely + baseline
  # this ensure internal consistency and that all values are estimated from the shapely values
  shap_sp_i$shap_suit_baseline <- shap_sp_i$shap_all_sum + baseline_value
  # the shap_suit_baseline object is also the predicted distribution
  
  # check the correlation is very high between methods, some error is expected from randomisation proceedure in the shapley estimations.
  suit_method_correlation <- cor(shap_sp_i$shap_suit_baseline, shap_sp_i$suitability, method = 'pearson')
  print(suit_method_correlation)
  
  
  ####
  ## DEFINE THE EXPECTED DISTRIBUTION USING THE BASELINE PREDICTION AND THE POSITIVE SHAPLEY VALUES, RELEASING THREATS
  ## here we require 1. the areas with positive values for natural factors, and within 
  ## this set, the 2. areas with positive values for threats (unthreatened), 3. areas that would be suitable if
  
  ## PARTITION
  # partition out the different shapley components in different groupings
  shap_sp_i$nn_sum      <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% nn_factors]), na.rm = T)
  shap_sp_i$nn_sum_mask <- ifelse(shap_sp_i$nn_sum < 0, NA, shap_sp_i$nn_sum)
  shap_sp_i$hab_sum     <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% hab_factors]), na.rm = T)
  shap_sp_i$con_sum     <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% con_factors]), na.rm = T)
  shap_sp_i$threat_continuous <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)]), na.rm = T)
    
  ## EXPECTED DISTRIBUTION
  # get positive shapley values for threats (i.e., where habitat are positive contributions) 
  pos_threat_shaps <- rowSums(abs(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])), na.rm = T)

  # add back in positive effects of threats to get the expected distribution
  shap_sp_i$expected_distribution <- shap_sp_i$nn_sum + pos_threat_shaps + baseline_value
  
  # mask the expected distribution 
  shap_sp_i$expected_distribution <- ifelse(shap_sp_i$nn_sum < 0, NA, shap_sp_i$expected_distribution)
  
  ## OBSERVED DISTRIBUTION
  # take the habitat suitability predictions within the natural niche area
  shap_sp_i$observed_distribution <- ifelse(is.na(shap_sp_i$nn_sum_mask), NA, shap_sp_i$shap_suit_baseline)
  
  
  ## CONTRIBUTIONS TO SHADOW DISTRIBUTIONS
  # get matrix of negative shapley values
  any_neg_threat_shaps <- st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])
  any_neg_threat_shaps[any_neg_threat_shaps >= 0] <- NA
  
  # get the areas where any threat is negative
  shap_sp_i$sum_presence_negative_threat <- as.numeric(apply(any_neg_threat_shaps < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_con <- as.numeric(apply(any_neg_threat_shaps[names(any_neg_threat_shaps) %in% c(con_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_habitat <- as.numeric(apply(any_neg_threat_shaps[names(any_neg_threat_shaps) %in% c(hab_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  
  # define the shadow distribution quantitatively as the observed / expected
  shap_sp_i$SD_OratioE <- shap_sp_i$observed_distribution / shap_sp_i$expected_distribution

  # define shadow distribution as the observed as a proportion of E
  shap_sp_i$SD_OpercentOfE <- (shap_sp_i$observed_distribution - shap_sp_i$expected_distribution) / shap_sp_i$expected_distribution
  
  # define the shadow distribution quantitatively as the observed - expected
  shap_sp_i$SD_OminusE <- shap_sp_i$observed_distribution - shap_sp_i$expected_distribution
  
  # define the shadow distribution as the areas in the expected distribution that have negative net threat effects
  shap_sp_i$SD_expected_threat_continuous <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$threat_continuous)
  
  # define the shadow distribution as the areas in the expected where any threat is negative, and the magnitude of this threat
  shap_sp_i$SD_expected_sum_presence_negative_threat <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_threat)
  shap_sp_i$SD_expected_sum_presence_negative_con <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_con)
  shap_sp_i$SD_expected_sum_presence_negative_habitat <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_habitat)
  
  # test maps of species different distribution types
  # tm_shape(shap_sp_i) + tm_fill(col = 'shap_suit_baseline')
  # tm_shape(shap_sp_i) + tm_fill(col = 'observed_distribution')
  # tm_shape(shap_sp_i) + tm_fill(col = 'expected_distribution')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_OminusE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_ThreatInE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_anyNegThreatInE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'expected_but_absent')
  # tm_shape(shap_sp_i) + tm_fill(col = 'observed_and_present')
  
  # NOTE:
  # the observed_distribution and expected_distribution have the same masking, 
  # such that the observed is simply the expected but with the effects of threat factors retained in the
  # habitat suitability score. In this way, we can directly compare these layers as they are both built from the 
  # partitioned suitability scores using the shapley values.
  # Note that we avoided metrics that use the threshold defined by the original suitability model. This value is inherently 
  # biased towards representing all the processes that affect a species distribution together, whereas our 
  # comparison of observed_distribution and expected_distribution are built from partitioning the contributions
  # based on different categories so is less biased towards the currently observed distribution, but we lack a clear way 
  # to define a 'presence' or 'absence' based on thresholds as is often reported in the literature. 

  # return the objects of interest
  return(shap_sp_i)
  
})

# get names
names(all_dd) <- sp_list


#### 7. Plot each distribution type across all species ----

tmap_mode('plot')

river_lake_tm <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F)

# create output directory
dir.create(paste0(fig_dir, '/shadow_dist_summaries/shadow_distribution_maps/'))

# create maps across all species
lapply(1:length(all_dd), function(x){
  
  sp_dd <- all_dd[[x]]
  
  pdf(paste0(fig_dir, '/shadow_dist_summaries/shadow_distribution_maps/', unique(sp_dd$species_name), '.pdf'),
      width = 10, height = 15,
      bg = "transparent"
  )
  
  print(tmap_arrange(
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              breaks = c(0, 0.5, 1),
              col = 'observed_distribution', 
              title = "observed suitability inside niche",
        legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              breaks = c(0, 0.5, 1),
              col = 'expected_distribution', 
              title = "expected suitability inside niche",
              legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              palette = 'Spectral',
              breaks = round(seq(min(sp_dd$SD_OpercentOfE, na.rm = T), 
                                 0,
                                 length.out = 4),2),
              col = 'SD_OpercentOfE', 
              title = "Shadow distribution",
              legend.is.portrait = F,
              midpoint = NA
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              col = 'SD_expected_threat_continuous', 
              title = "",
              legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              col = 'SD_expected_sum_presence_negative_threat', 
              title = "",
              legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              col = 'SD_expected_sum_presence_negative_con', 
              title = "",
              legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    tm_shape(sp_dd) +
      tm_fill(style = 'cont',
              col = 'SD_expected_sum_presence_negative_habitat', 
              title = "",
              legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = unique(sp_dd$species_name)) + 
      river_lake_tm,
    
    
    
    
    ncol = 2 ))
  
  dev.off()
  
  })



#### 8. Aggregate together shadow distribution values ----

# bind together all outputs
all_dd_bind <- bind_rows(all_dd)

# remove areas from distributions where they are not naturally predicted to occur
all_dd_bind <- all_dd_bind %>% filter(!is.na(expected_distribution))

# summarise the distribution properties per subcatchment
dist_properties <- c('observed_distribution', 
                     'expected_distribution',
                     'SD_expected_threat_continuous', 
                     'SD_expected_sum_presence_negative_con', 
                     'SD_expected_sum_presence_negative_habitat', 
                     'SD_OratioE',
                     'SD_OpercentOfE')

# filter to relevant columns, 
# sum the number of species in each subcatchment, 
# filter to be more than 2 species per subcatchment.
all_dd_toSum <- all_dd_bind %>% 
  select(TEILEZGNR, species_name, 
         any_of(dist_properties)) %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  nest() %>% 
  mutate(n_sp = purrr::map(data, ~sum(!is.na(.$expected_distribution), na.rm = T))) %>% 
  unnest(c(data, n_sp)) %>% 
  filter(n_sp > 2)

# create summaries of all focal columns (mean, min, max, sd)
all_dd_sum <- all_dd_toSum %>% 
    summarise_at(c(dist_properties, 'n_sp'),
               list(mean = mean, min = min, max = max, sd = sd), 
               na.rm = T) %>% 
  left_join(., all_dd_bind %>% select(TEILEZGNR) %>% unique())

# convert back to sf object for plotting
all_dd_sum_sf <- all_dd_sum %>% st_as_sf()

# make baseplot for rivers
river_base <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F) + 
  tm_layout(frame = F,
            bg.color = "transparent")

# function to standardize plotting
map_values <- function(x, 
                       breaks = NULL,
                       palette, 
                       legend.reverse = F){tm_shape(subcatchments_rhine_union_2) + 
    tm_polygons(col = 'gray90') + 
    #tm_shape(all_dd_sum_sf %>% sample_n(., 1000)) + 
    tm_shape(all_dd_sum_sf) + 
    tm_fill(col = x, 
            style = 'cont', 
            title = '', 
            legend.is.portrait = T,
            legend.reverse = legend.reverse,
            palette = palette, 
            n = 10, 
            contrast = c(0, 1), 
            colorNA = 'transparent', 
            textNA = "", 
            breaks = breaks, 
            midpoint = NA) + 
    river_base + 
    tm_shape(subcatchments_rhine_union_2) + 
    tm_borders(col = 'black') + 
    tm_layout(legend.text.size = 1,
              legend.outside = F,
              legend.position = c('left', 'top'))}


# plots of percentage difference between observed and expected distributions 
mean_completeness_plot <- map_values('SD_OpercentOfE_mean', breaks = signif(seq(min(all_dd_sum[,c('SD_OpercentOfE_mean')], na.rm = T), 
                                                                       0, 
                                                                       length.out = 4), 2), 
                                     palette = 'Spectral', legend.reverse = T)

min_completeness_plot <- map_values('SD_OpercentOfE_min', breaks = signif(seq(min(all_dd_sum[,c('SD_OpercentOfE_min')], na.rm = T), 
                                                                          max(all_dd_sum[,c('SD_OpercentOfE_min')], na.rm = T), 
                                                                          length.out = 4), 2), 
           palette = 'Spectral', legend.reverse = T)

sd_completeness_plot <- map_values('SD_OpercentOfE_sd', breaks = signif(seq(min(all_dd_sum[,c('SD_OpercentOfE_sd')], na.rm = T), 
                                                                        max(all_dd_sum[,c('SD_OpercentOfE_sd')], na.rm = T), 
                                                                        length.out = 4), 2), 
                                    palette = 'Spectral', legend.reverse = T)


range <- signif(seq(min(all_dd_sum[,c('observed_distribution_mean', 'expected_distribution_mean')], na.rm = T), 
                    max(all_dd_sum[,c('observed_distribution_mean', 'expected_distribution_mean')], na.rm = T), 
                    length.out = 4), 2)

png(paste0(fig_dir, '/dark_diversity/all_maps.png'),
    width = 4000, height = 2250, res = 300,
    bg = "transparent")
tmap_arrange(map_values('observed_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             map_values('expected_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             map_values('expected_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             
             mean_completeness_plot,
             min_completeness_plot, 
             sd_completeness_plot,
             
             map_values('SD_expected_sum_presence_negative_con_mean', breaks = round(seq(0, 1, length.out = 4),2), palette = '-Spectral', legend.reverse = T),
             map_values('SD_expected_sum_presence_negative_habitat_mean', breaks = round(seq(0, 4, length.out = 4), 2), palette = '-Spectral', legend.reverse = T), 
             map_values('SD_expected_sum_presence_negative_habitat_mean', breaks = round(seq(0, 4, length.out = 4), 2), palette = '-Spectral', legend.reverse = T),
            
             nrow = 3, ncol = 3)
dev.off()


#### 9. Data summaries of each shadow distribution property ----

## Mean suitability
# Summary of mean observed distribution
round(summary(all_dd_sum_sf$observed_distribution_mean), 2)
mean_suit <- mean(all_dd_sum_sf$observed_distribution_mean, na.rm = T)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.21    0.47    0.56    0.56    0.66    0.96 

# Summary of mean expected distribution
round(summary(all_dd_sum_sf$expected_distribution_mean), 2)
mean_exp <- mean(all_dd_sum_sf$expected_distribution_mean, na.rm = T)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.35    0.54    0.62    0.62    0.70    0.99 

# Mean reduction in suitability, calculate directly from species-metric
1 - mean(all_dd_sum_sf$SD_OratioE_mean) # 0.1078135% reduction in suitability 
t.test(all_dd_sum_sf$observed_distribution_mean, all_dd_sum_sf$expected_distribution_mean)


## Minimum suitability
# Summary of minimum observed distribution
round(summary(all_dd_sum_sf$observed_distribution_min), 2)
mean_suit <- mean(all_dd_sum_sf$observed_distribution_min,na.rm = T)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.13    0.28    0.34    0.35    0.42    0.74 

# Summary of minimum expected distribution
round(summary(all_dd_sum_sf$expected_distribution_min), 2)
mean_exp <- mean(mean(all_dd_sum_sf$expected_distribution_min,na.rm = T))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.30    0.40    0.45    0.45    0.50    0.74 


## Percentage difference between observed and expected
# Mean percentage reduction in suitability comparing observed to expected 
1-mean(all_dd_sum_sf$SD_OratioE_min) # 0.2911247% reduction in suitability 
t.test(all_dd_sum_sf$observed_distribution_min, all_dd_sum_sf$expected_distribution_min)

# Summary of mean percentage difference between observed and expected suitability
summary(all_dd_sum_sf$SD_OpercentOfE_mean)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.53072 -0.13891 -0.08651 -0.10781 -0.05085  0.00000 

# Summary of minmum percentage difference between observed and expected suitability
summary(all_dd_sum_sf$SD_OpercentOfE_min)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.7390 -0.4778 -0.2367 -0.2911 -0.1201  0.0000 


## Percentage difference in worst 10% of catchments 
# What is the mean reduction in suitability in the lowest 10th quantile 
mean(all_dd_sum_sf$SD_OpercentOfE_mean[all_dd_sum_sf$SD_OpercentOfE_mean<quantile(all_dd_sum_sf$SD_OpercentOfE_mean, 0.1)])
quantile(all_dd_sum_sf$SD_OpercentOfE_min, 0.1)


## Maximum reduction in habitat suitability
# across-catchment minimum of mean reduction
min(all_dd_sum_sf$SD_OratioE_mean)
# across-catchment minimum of minimum reduction
min(all_dd_sum_sf$SD_OratioE_min)


## Comparative correlations between threat numbers and percentage reductions
# spearmans rank between shadow distribution and threat number
cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_con_mean, all_dd_sum_sf$SD_OratioE_mean, method = 'spearman')
cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_habitat_mean, all_dd_sum_sf$SD_OratioE_mean, method = 'spearman')

# spearmans rank between shadow distribution and threat number
cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_con_mean, all_dd_sum_sf$SD_OratioE_min, method = 'spearman')
cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_habitat_mean, all_dd_sum_sf$SD_OratioE_min, method = 'spearman')


png(paste0(fig_dir, '/dark_diversity/biplot_threats.png'),
    width = 750, height = 1500, res = 300,
    bg = "transparent")
grid.arrange(
  ggplot(all_dd_sum_sf %>% arrange(-SD_OratioE_mean)) + 
  geom_point(aes(x = expected_distribution_mean, y = observed_distribution_mean, col = SD_OratioE_mean)) + 
  geom_abline(col = 'black', lty = 2) + 
  theme_bw() + 
  ylab('observed suitability') + 
  xlab('expected suitability') + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        plot.background = element_blank(),
        axis.line = element_line(), 
        axis.ticks = element_blank(), 
        legend.position = 'none', 
        aspect.ratio = 0.5) + 
  scale_colour_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral")),
  
  ggplot(all_dd_sum_sf %>% arrange(-SD_OratioE_mean)) + 
    geom_point(aes(y = SD_expected_sum_presence_negative_habitat_mean, 
                   x = SD_expected_sum_presence_negative_con_mean, 
                   col = SD_OratioE_mean)) + 
    theme_bw() + 
    ylab('net habitat threat') + 
    xlab('net connectivity threat') + 
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          plot.background = element_blank(),
          axis.line = element_line(), 
          axis.ticks = element_blank(), 
          legend.position = 'none', 
          aspect.ratio = 0.5) + 
    scale_colour_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral")), 
  
  ncol = 1, nrow = 2
  
)
dev.off()

  

#### 10. Summarise environmental properties catchments with strong vs. weak shadow distribution ----

# identify high completeness regions
subcatch_complete100 <- all_dd_sum %>% 
  filter(SD_OratioE_mean > 0.95) %>% 
  pull(TEILEZGNR) %>% 
  unique() 

# identify low completeness regions
subcatch_complete_others <- all_dd_sum %>% 
  filter(SD_OratioE_mean < 0.95) %>% 
  pull(TEILEZGNR) %>% 
  unique()

# what are the properties of these streams
data.frame(high_completeness = all_env_subcatchments %>% 
        filter(TEILEZGNR %in% subcatch_complete100) %>% 
        select(vars) %>% 
        summarise_all(., function(x) round(mean(x),2)) %>% t,
      
      low_completeness = all_env_subcatchments %>% 
        filter(TEILEZGNR %in% subcatch_complete_others) %>% 
        select(vars) %>% 
        summarise_all(., function(x) round(mean(x),2)) %>% t) %>% 
  view()



  

