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
nn_factors <- c('ecoF_discharge_max_log10_SHAP','ecoF_slope_min_log10', 'stars_t_mn_m_c_SHAP', 
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

#### 5. Map shadow distributions for each species ----

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


#### 6. ----


#### 4A GET RASTER MAPS OF SHAPLEY VALUES ----

all_dd <- lapply(1:length(sp_list), function(i){
  
  # get species
  sp <- sp_list[i]
  print(sp)
  
  # Check file exists for shapley raster object and baseline random forest
  do_pa <- file.exists(paste0(shap_dirs[grep(sp, shap_dirs)], "/shap_raster_pa.TIF"))
  do_rf <- file.exists(paste0(rf_pa[grep(sp, rf_pa)]))
  if(do_pa != T | do_rf != T){next()}
  
  ## READING IN SHAPELY VALUES AND MODEL OBJECTS FOR FURTHER ANALYSIS
  # read in shapley values for a given species
  shap_rast_pa_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
  # take all factors
  shap_rast_pa_i <- shap_rast_pa_i[[which(names(shap_rast_pa_i) %in% all_factors)]]
  # get shapley values for natural factors only
  nn_shap_rast <- shap_rast_pa_i[[which(names(shap_rast_pa_i) %in% nn_factors)]]
  # get shapley values for threats
  threat_shap_rast <- shap_rast_pa_i[[which(names(shap_rast_pa_i) %in% c(con_factors, hab_factors))]]
  # remove non_na
  shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]
  # order by name
  shap_rast_pa_i <- shap_rast_pa_i[[sort(names(shap_rast_pa_i))]]
  
  ## GET BASELINE VALUES TO CONVERT SHAPLEY VALUES TO SUITABILITIES THROUGH ADDITION TO REFERENCE BASELINE
  # read in baseline model and obtain baseline values for prediction
  # this value indicates the deviation from which the shapley values predict
  baseline_model <- readRDS(rf_pa[grep(sp, rf_pa)])
  baseline_value <- mean(as.numeric(predict(baseline_model, type = 'prob')[,2]), na.rm = T)
  
 #  # ORGANISE SPECIES PRESENCE OR ABSENCE
 model_data   <- readRDS(paste0(sdm_dirs[i], '/data/sdm_input_data.rds'))@pa_data$full_data
 sp_rast <- rast(sp_raster_suit_pa[i])
 sp_rast <- resample(sp_rast, shap_rast_pa_i)
 #  pred_occ <- terra::extract(sp_rast, model_data[c('X', 'Y')])
 #  threshold <- ecospat::ecospat.max.tss(pred_occ[,2], model_data$occ)$max.threshold
 #  presence_map <- sp_rast
 #  absence_map  <- sp_rast
 #  # define presence map
 #  presence_map[sp_rast<threshold] <- 0
 #  presence_map[sp_rast>threshold] <- 1
 #  # define absence map
 #  absence_map[sp_rast<threshold] <- 1
 #  absence_map[sp_rast>threshold] <- 0
 #  # align maps
 #  shap_rast_pa_i <- resample(shap_rast_pa_i, sp_rast)
 #  absence_map[is.na(shap_rast_pa_i[[1]])]  <- NA
 #  presence_map[is.na(shap_rast_pa_i[[1]])] <- NA
  
  
  ## DEFINE THE OBSERVED DISTRIBUTION USING THE BASELINE PREDICTION AND THE SHAPLEY VALUES
  # sum all the shapley values
  shap_all_sum <- app(shap_rast_pa_i, sum, na.rm = T)
  
  # add the baseline value to the sum of the shapley values to obtain a reconstruction of the habitat suitability
  shap_suit_baseline <- shap_all_sum + baseline_value
  # the shap_suit_baseline object is also the observed distribution
  
  # check the correlation is very high between methods, some error is expected from randomisation proceedure in the shapley estimations.
  suit_method_correlation <- layerCor(c(shap_suit_baseline, sp_rast), fun = 'pearson', na.rm = T)
  print(suit_method_correlation[[1]])
  # define threshold for identifying range limits
  pred_occ <- terra::extract(shap_suit_baseline, model_data[c('X', 'Y')])
  threshold <- ecospat::ecospat.max.tss(pred_occ[,2], model_data$occ)$max.threshold
  
  # shapley defined absences
  absence <- shap_suit_baseline
  absence[shap_suit_baseline<threshold]  <- 1
  absence[shap_suit_baseline>=threshold] <- 0
  
  ## DEFINE THE EXPECTED DISTRIBUTION USING THE BASELINE PREDICTION AND THE POSITIVE SHAPLEY VALUES, RELEASING THREATS
  ## here we require 1. the areas with positive values for natural factors, and within 
  ## this set, the 2. areas with positive values for threats (unthreatened), 3. areas that would be suitable if
  
  # expected distribution by masking suitability to only positive shapley values
  nn_sum <- app(nn_shap_rast, sum, na.rm = T)
  nn_sum_mask <- nn_sum
  nn_sum_mask[nn_sum_mask<0] <- NA
  
  # take the distribution within the natural niche
  nn_dist <- shap_suit_baseline
  nn_dist[is.na(nn_sum_mask)] <- NA
  
  # add back in positive and negative threat effects (represent unthreatened areas and removal of threats respectively)
  nn_dist_corrected <-  app(nn_shap_rast, sum, na.rm = T) + app(abs(threat_shap_rast), sum, na.rm = T)
  nn_dist_corrected[is.na(nn_sum_mask)] <- NA
  
  # add to natural distribution and corrected distribution the baseline
  nn_dist_corrected <- nn_dist_corrected + baseline_value
  
  # define expected distribution object
  expected_distribution <- nn_dist_corrected
  
  
  #### DEFINE A STRICTER VERSION OF THE EXPECTED DISTRIBUTION WHERE ALL NICHE VALUES MUST BE POSITIVE
  # get number of niche layers

  threshold_shapley <- nn_sum[absence==0] %>% mean(., na.rm=T)
  nn_sum_mask_STRICT <- nn_sum
  nn_sum_mask_STRICT[nn_sum<threshold_shapley] <- NA
   
  # get the shapley values in the strict natural niche
  # add back in positive and negative threat effects (represent unthreatened areas and removal of threats respectively)
  nn_dist_corrected_STRICT <-  app(nn_shap_rast, sum, na.rm = T) + app(abs(threat_shap_rast), sum, na.rm = T)
  nn_dist_corrected_STRICT[is.na(nn_sum_mask_STRICT)] <- NA
  nn_dist_corrected_STRICT <- nn_dist_corrected_STRICT + baseline_value
  expected_distribution_STRICT <- nn_dist_corrected_STRICT
  
  # defining a threshold here is problematic because the observational data are influenced by mutliple processes, 
  # so is biased and does not then match up conceptually with any corrected expected distributions. 
  
  # define observed distribution
  observed_distribution <- shap_suit_baseline
  observed_distribution[is.na(expected_distribution)] <- NA
  
  # define strict version of observed distribution
  observed_distribution_STRICT <- shap_suit_baseline
  observed_distribution_STRICT[is.na(expected_distribution_STRICT)] <- NA
  
  # define observed and present by removing locations predicted as absences
  observed_and_present <- shap_suit_baseline
  observed_and_present[absence!=0] <- NA
  
  # find areas of the expected distribution that are absences in the modelled distribution
  expected_but_absent <- expected_distribution
  expected_but_absent[absence!=1] <- NA
  
  # define strict version of expected but absent distribution
  expected_but_absent_STRICT <- expected_distribution_STRICT
  expected_but_absent_STRICT[absence!=1] <- NA
  
  # NOTE:
  # the observed_distribution and expected_distribution have the same masking, 
  # such that the observed is simply the expected but with the effects of threat factors retained in the
  # habitat suitability score. In this way, we can directly compare these layers as they are both build from the 
  # partitioned suitability scores using the shapley values.
  # In contrast, the observed_and_present is the observed_distribution but removing all locations that are absences
  # and the expected_but_absence are the expected_distribution but removing all locations that are presences. The 
  # removal of locaitons is based on a thresholded suitability score from the original data, so is inherently 
  # biased towards representing all the processes that affect a species distribution together, whereas our 
  # comparison of observed_distribution and expected_distribution are built from partitioning the contributions
  # based on different categories so is unbiased, but we lack a clear way to define a 'presence' or 'absence' 
  # as is often reported in the literature. 
  
  # use of these objects in community stacking:
  # the stack of of observed and present can give the species richness
  # the stack of expected but absent can give the dark diversity
  # the stack of expected distribution can dive the expected diversity
  # the difference between expected and observed distributions represents the shadow distribution
  
  # defining the shadow distribution as the different between the expected distribution and the observed distribution
  shadow_distribution = expected_distribution - observed_distribution
  # higher values means a greater shadow and anthropic influence, 
  # negative values indicate where human influences are supporting species where they would normally 
  # be less suitable 
  
  # define strict version of above
  shadow_distribution_STRICT = expected_distribution_STRICT - observed_distribution_STRICT
  
  
  # make names cleaner for processing afterwards
  names(observed_distribution) <- sp; names(observed_distribution_STRICT) <- sp
  names(expected_distribution) <- sp; names(expected_distribution_STRICT) <- sp
  names(expected_but_absent)   <- sp; names(expected_but_absent_STRICT)   <- sp
  names(observed_and_present)  <- sp
  names(shadow_distribution)   <- sp; names(shadow_distribution_STRICT)   <- sp
  
  ## return the objects of interest
  return(list(observed_distribution = observed_distribution, 
              expected_distribution = expected_distribution, 
              expected_but_absent   = expected_but_absent, 
              observed_and_present  = observed_and_present, 
              shadow_distribution   = shadow_distribution, 
              observed_distribution_STRICT = observed_distribution_STRICT, 
              expected_distribution_STRICT = expected_distribution_STRICT, 
              expected_but_absent_STRICT = expected_but_absent_STRICT, 
              shadow_distribution_STRICT = shadow_distribution_STRICT))
  
})

# get names
names(all_dd) <- sp_list

#### 4Bi MAKE PLOT OF ALL DISTRIBUTION TYPES FOR A GIVEN SPECIES ----

tmap_mode('plot')

river_lake_tm <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F)

pdf(paste0(fig_dir, "/dark_diversity/species_maps.pdf"),
  width = 10, height = 15,
  bg = "transparent"
)

for (i in 1:length(names(all_dd))) {
  sp_dist_all <- rast(all_dd[[i]])

  print(tmap_arrange(
    tm_shape(sp_dist_all$observed_distribution) +
      tm_raster(
        style = "cont", title = "observed suitability inside niche",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = names(all_dd)[i]) + 
      river_lake_tm,
    tm_shape(sp_dist_all$expected_distribution) +
      tm_raster(
        style = "cont", title = "expected suitability inside niche",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$expected_but_absent) +
      tm_raster(
        style = "cont", title = "inside niche, but absent",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$observed_and_present) +
      tm_raster(
        style = "cont", title = "inside niche, and present",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$shadow_distribution) +
      tm_raster(
        style = "cont", title = "shadow distribution",
        legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_shape(sp_dist_all$expected_but_absent) +
      tm_raster(palette = "red", legend.show = F) + 
      river_lake_tm,
    nrow = 3, ncol = 2
  ))
}
dev.off()

#### 4Bii MAKE PLOT OF ALL DISTRIBUTION TYPES FOR A GIVEN SPECIES BASED ON STRICT NICHE CRITERIA ----

tmap_mode('plot')

river_lake_tm <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F)

pdf(paste0(fig_dir, "/dark_diversity/species_maps_STRICT_NICHE.pdf"),
    width = 10, height = 15,
    bg = "transparent"
)

for (i in 1:length(names(all_dd))) {
  sp_dist_all <- rast(all_dd[[i]])
  
  print(tmap_arrange(
    tm_shape(sp_dist_all$observed_distribution_STRICT) +
      tm_raster(
        style = "cont", title = "observed suitability inside niche",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_layout(main.title = names(all_dd)[i]) + 
      river_lake_tm,
    tm_shape(sp_dist_all$expected_distribution_STRICT) +
      tm_raster(
        style = "cont", title = "expected suitability inside niche",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$expected_but_absent_STRICT) +
      tm_raster(
        style = "cont", title = "inside niche, but absent",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$observed_and_present) +
      tm_raster(
        style = "cont", title = "inside niche, and present",
        breaks = c(0, 0.5, 1), legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) + 
      river_lake_tm,
    tm_shape(sp_dist_all$shadow_distribution_STRICT) +
      tm_raster(
        style = "cont", title = "shadow distribution",
        legend.is.portrait = F
      ) +
      tm_legend(scale = 1.5) +
      tm_shape(sp_dist_all$expected_but_absent_STRICT) +
      tm_raster(palette = "red", legend.show = F) + 
      river_lake_tm,
    nrow = 3, ncol = 2
  ))
}
dev.off()


#### 4Ci AGGREGATE RESULTING DISTRIBUTIONS AS AN AVERAGE ACROSS ALL SPECIES ----

sp_dist_summaries <- lapply(1:length(all_dd), function(x){
  
  # the difference between expected and observed distributions
  cont_shadow <- all_dd[[x]]$shadow_distribution
  
  # get the average reduction in suitability in the shadow distribution relative to expected distribution (i.e., expected - observed) 
  shadow_prop_expected <- mean(values(cont_shadow) / values(all_dd[[x]]$expected_distribution), na.rm = T)
  
  # get the area comparisons of expected but absent and in the expected distribution
  # get values for expected but absent
  e_b_a <- values(all_dd[[x]]$expected_but_absent)

  # whats the proportion across all expected but absent locations
  prop_absence_from_expected <- sum(e_b_a, na.rm = T) / sum(values(all_dd[[x]]$expected_distribution), na.rm = T)
  
  # whats the proportion in areas only absences
  shadow_prop_in_absence <- (sum(e_b_a, na.rm = T) - sum(values(all_dd[[x]]$observed_distribution)[!is.na(e_b_a)], na.rm = T)) / sum(e_b_a, na.rm = T)
  
  return(data.frame(species = names(all_dd)[x], 
                    shadow_prop_expected = shadow_prop_expected, 
                    prop_absence_from_expected = prop_absence_from_expected, 
                    shadow_prop_in_absence = shadow_prop_in_absence)) 
  })

# bind together outputs
sp_dist_summaries <- bind_rows(sp_dist_summaries)

# remove non-native species 
sp_dist_summaries <- sp_dist_summaries %>% filter(species != 'Oncorhynchus mykiss')

summary(sp_dist_summaries$shadow_prop_expected)
sd(sp_dist_summaries$shadow_prop_expected)
summary(sp_dist_summaries$shadow_prop_in_absence)
sd(sp_dist_summaries$shadow_prop_in_absence)
summary(sp_dist_summaries$prop_absence_from_expected)
sd(sp_dist_summaries$prop_absence_from_expected)




#### 4Cii AGGREGATE RESULTING DISTRIBUTIONS AS AN AVERAGE ACROSS ALL SPECIES ----

sp_dist_summaries_STRICT <- lapply(1:length(all_dd), function(x){
  
  # the difference between expected and observed distributions
  cont_shadow <- all_dd[[x]]$shadow_distribution_STRICT
  
  # get the average reduction in suitability in the shadow distribution relative to expected distribution (i.e., expected - observed) 
  shadow_prop_expected <- mean(values(cont_shadow) / values(all_dd[[x]]$expected_distribution_STRICT), na.rm = T)
  
  # get the area comparisons of expected but absent and in the expected distribution
  # get values for expected but absent
  e_b_a <- values(all_dd[[x]]$expected_but_absent_STRICT)
  
  # whats the proportion across all expected but absent locations
  prop_absence_from_expected <- sum(e_b_a, na.rm = T) / sum(values(all_dd[[x]]$expected_distribution_STRICT), na.rm = T)
  
  # whats the proportion in areas only absences
  shadow_prop_in_absence <- (sum(e_b_a, na.rm = T) - sum(values(all_dd[[x]]$observed_distribution_STRICT)[!is.na(e_b_a)], na.rm = T)) / sum(e_b_a, na.rm = T)
  
  return(data.frame(species = names(all_dd)[x], 
                    shadow_prop_expected = shadow_prop_expected, 
                    prop_absence_from_expected = prop_absence_from_expected, 
                    shadow_prop_in_absence = shadow_prop_in_absence)) 
})

# bind together outputs
sp_dist_summaries_STRICT <- bind_rows(sp_dist_summaries_STRICT)

# remove non-native species 
sp_dist_summaries_STRICT <- sp_dist_summaries_STRICT %>% filter(species != 'Oncorhynchus mykiss')

summary(sp_dist_summaries_STRICT$shadow_prop_expected)
sd(sp_dist_summaries_STRICT$shadow_prop_expected)
summary(sp_dist_summaries_STRICT$shadow_prop_in_absence)
sd(sp_dist_summaries_STRICT$shadow_prop_in_absence)
summary(sp_dist_summaries_STRICT$prop_absence_from_expected)
sd(sp_dist_summaries_STRICT$prop_absence_from_expected)



#### COMMUNITY LEVEL DISTRIBUTION SUMMARIES (STACKING DISTRIBUTIONS) ----

## REMOVE NON-NATIVE SPECIES
all_dd <- all_dd[-which(names(all_dd) == 'Oncorhynchus mykiss')]

## SUM RESULTING DISTRIBUTIONS OVER ALL SPECIES
# compile observed distributions
sum_observed_dist <- app(rast(lapply(all_dd, function(x) x$observed_distribution)), sum, na.rm = T)
# compile expected distributions summing shapley values of natural niche
sum_expected_dist<- app(rast(lapply(all_dd, function(x) x$expected_distribution)), sum, na.rm = T)
# compile dark distributions based on summation of shapley values of natural niche
sum_expected_but_absent <- app(rast(lapply(all_dd, function(x) x$expected_but_absent)), sum, na.rm = T)
# compile sum of observed and present
sum_observed_and_present <- app(rast(lapply(all_dd, function(x) x$observed_and_present)), sum, na.rm = T)
# compile sum of shadow distributions
mean_shadow_distribution <- app(rast(lapply(all_dd, function(x) x$shadow_distribution)), mean, na.rm = T)
sum_shadow_distribution <- app(rast(lapply(all_dd, function(x) x$shadow_distribution)), sum, na.rm = T)

# load in tmap to view
library(tmap)
tmap_mode('view')
tm_shape(sum_shadow_distribution) + 
  tm_raster()

# look at rough plots: OBSERVED DISTRIBUTION
plot(sum_observed_dist)
# look at rough plots: OBSERVED DISTRIBUTION
plot(sum_expected_dist)
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
  tm_raster(palette = 'Spectral', legend.is.portrait = F, n = 5) +
  river_base
expected_diversity_plot

dark_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(sum_shadow_distribution, raster.downsample = F) +
  tm_raster(palette = 'Spectral', legend.is.portrait = F, n = 5) +
  river_base
dark_diversity_plot

completeness_v2 <- completeness
completeness_v2[sum_expected_dist<=2] <- NA
completeness_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(completeness_v2, raster.downsample = F) +
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

