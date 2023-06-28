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
  
  # get species
  sp <- all_shadow[[i]]$species_name[1]
  print(sp)
  
  # get shapley shape files distribution
  shap_sp_i <- all_shadow[[i]]
  
  ## DEFINE THE OBSERVED DISTRIBUTION USING THE BASELINE PREDICTION AND THE SHAPLEY VALUES
  
  # get the baseline value as the mean of all habitat suitability predictions
  baseline_value <- mean(shap_sp_i$suitability, na.rm = T)
    
  # sum all the shapley values
  shap_sp_i$shap_all_sum <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% all_factors]), na.rm = T)
  
  # get the baseline suitability values estimated as the sum of the shapely + baseline
  # this ensure internal consistency and that all values are estimated from the shapely values
  shap_sp_i$shap_suit_baseline <- shap_sp_i$shap_all_sum + baseline_value
  # the shap_suit_baseline object is also the observed distribution
  
  # check the correlation is very high between methods, some error is expected from randomisation proceedure in the shapley estimations.
  suit_method_correlation <- cor(shap_sp_i$shap_suit_baseline, shap_sp_i$suitability, method = 'pearson')
  print(suit_method_correlation)
  
  # shapley defined absences
  shap_sp_i$shap_presence_baseline <- ifelse(shap_sp_i$shap_suit_baseline >= shap_sp_i$threshold, 1, 0)
    
  ## DEFINE THE EXPECTED DISTRIBUTION USING THE BASELINE PREDICTION AND THE POSITIVE SHAPLEY VALUES, RELEASING THREATS
  ## here we require 1. the areas with positive values for natural factors, and within 
  ## this set, the 2. areas with positive values for threats (unthreatened), 3. areas that would be suitable if
  
  # expected distribution by masking suitability to only positive shapley values
  shap_sp_i$nn_sum      <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% nn_factors]), na.rm = T)
  shap_sp_i$nn_sum_mask <- ifelse(shap_sp_i$nn_sum < 0, NA, shap_sp_i$nn_sum)
  
  # get positive shapley values for threats (i.e., where habitat is good also) and add to expected distribution
  pos_threat_shaps <- st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])
  pos_threat_shaps[pos_threat_shaps < 0] <- NA
  
  # add back in positive effects of threats
  shap_sp_i$nn_sum <- shap_sp_i$nn_sum + rowSums(pos_threat_shaps, na.rm = T)
  shap_sp_i$nn_sum_mask <- shap_sp_i$nn_sum_mask + rowSums(pos_threat_shaps, na.rm = T)
  
  # take the distribution within the natural niche
  shap_sp_i$nn_dist <- ifelse(is.na(shap_sp_i$nn_sum_mask), NA, shap_sp_i$shap_suit_baseline)
  
  # get sum of threat effects 
  # add back in positive and negative threat effects (represent unthreatened areas and removal of threats respectively)
  neg_threat_shaps <- st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])
  neg_threat_shaps[neg_threat_shaps >= 0] <- NA
  # get sum of negative threats
  sum_negative_threats <- rowSums(neg_threat_shaps, na.rm = T)
  shap_sp_i$nn_sum_mask_minus_threat <- shap_sp_i$nn_sum_mask + sum_negative_threats

  # get sum of the threats as a single layer
  shap_sp_i$threat_continuous <- sum_negative_threats
  
  # get the areas where any threat is negative
  shap_sp_i$sum_presence_negative_threat <- as.numeric(apply(neg_threat_shaps < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_con <- as.numeric(apply(neg_threat_shaps[names(neg_threat_shaps) %in% c(con_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_habitat <- as.numeric(apply(neg_threat_shaps[names(neg_threat_shaps) %in% c(hab_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  
  # add to natural distribution and corrected distribution the baseline
  shap_sp_i$expected_distribution <- shap_sp_i$nn_sum_mask + baseline_value
  
  # define observed distribution
  shap_sp_i$observed_distribution <- shap_sp_i$nn_sum_mask_minus_threat + baseline_value

  # define observed and present by removing locations predicted as absences
  shap_sp_i$observed_and_present <- ifelse(shap_sp_i$shap_presence_baseline == 0, NA, shap_sp_i$shap_suit_baseline)

  # find areas of the expected distribution that are absences in the modelled distribution
  shap_sp_i$expected_but_absent <- ifelse(shap_sp_i$shap_presence_baseline == 0, shap_sp_i$expected_distribution, NA)

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
  # In contrast, the observed_and_present is the observed_distribution but removing all locations that are absences
  # and the expected_but_absence are the expected_distribution but removing all locations that are presences. The 
  # removal of locaitons is based on a thresholded suitability score from the original data, so is inherently 
  # biased towards representing all the processes that affect a species distribution together, whereas our 
  # comparison of observed_distribution and expected_distribution are built from partitioning the contributions
  # based on different categories so is unbiased, but we lack a clear way to define a 'presence' or 'absence' 
  # as is often reported in the literature. 
  # Probably the most important layer here is the anyNegThreatInE (any negative threat in expected distribution) 
  # as this defines whether a species is exposed, and responding, negatively to any threat within its ecological niche.
  # The threat in E is also interesting as it approximates the net effects of threats within the expected distribution. 
  
  # use of these objects are used in community stacking:
  # the stack of of observed and present can give the species richness
  # the stack of expected but absent can give the predicted dark diversity
  # the stack of expected distribution can dive the expected diversity
  # the sum of the different shadow distribution components gives the total effect of threats on diversity
  
  ## return the objects of interest
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

# threshold the expected distribution
all_dd_bind$shap_suitability_threshold <- ifelse(all_dd_bind$shap_suit_baseline < threshold, NA, all_dd_bind$shap_suit_baseline)
all_dd_bind$expected_distribution_threshold <- ifelse(all_dd_bind$expected_distribution < threshold, NA, all_dd_bind$expected_distribution)

# summarise the distribution properties per subcatchment
dist_properties <- c('suitability', 
                     'expected_distribution',
                     'observed_distribution', 
                     'SD_expected_threat_continuous', 
                     'SD_expected_sum_presence_negative_threat',
                     'SD_expected_sum_presence_negative_con', 
                     'SD_expected_sum_presence_negative_habitat')
all_dd_sum <- all_dd_bind %>% 
  select(TEILEZGNR, species_name, 
         any_of(dist_properties)) %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  summarise_at(dist_properties, 
               sum, 
               na.rm = T) %>% 
  left_join(., all_dd_bind %>% select(TEILEZGNR) %>% unique()) %>% 
  # calculate community completeness here
  mutate(community_completeness = (.$observed_distribution / .$expected_distribution)*100, 
         community_completeness_threshold = (.$suitability / expected_distribution)*100) %>% 
  mutate(community_completeness = ifelse(.$community_completeness > 100, 100, .$community_completeness), 
         community_completeness_threshold  = ifelse(.$community_completeness_threshold > 100, 100, .$community_completeness_threshold)) %>% 
  # filter out areas outside the ecological niche of all species
  filter(!is.na(expected_distribution))

# calculate the weighted average shapley value
all_dd_sum <- all_dd_bind %>% 
  select(TEILEZGNR, species_name, 
         expected_distribution, shadow_distribution_OminusE) %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  do(shadow_distribution_OminusE_weighted_mean = weighted.mean(.$shadow_distribution_OminusE, .$expected_distribution, na.rm = T)) %>% 
  unnest(c(shadow_distribution_OminusE_weighted_mean)) %>% 
  ungroup() %>% 
  left_join(all_dd_sum %>% unique(), .)
  
# turn to NA all areas where we do not expect species
all_dd_sum[all_dd_sum$expected_distribution == 0,dist_properties] <- NA

# get summaries of each distribution type
all_dd_sum %>% 
   filter(expected_distribution != 0, !is.na(expected_distribution)) %>% 
   summarise_at(c(dist_properties,'community_completeness'), mean, na.rm = T)

# plot the summarised distributions
all_dd_sum_sf <- all_dd_sum %>% st_as_sf()

## make baseplot for rivers
river_base <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F) + 
  tm_layout(frame = F,
            bg.color = "transparent")

map_values <- function(x){print(tm_shape(subcatchments_rhine_union_2) + 
    tm_polygons(col = 'gray90') + 
    tm_shape(all_dd_sum_sf) + 
    tm_fill(col = x, 
            style = 'cont', 
            title = '', 
            legend.is.portrait = T,
            legend.reverse = T,
            palette = "Spectral", 
            n = 10, 
            contrast = c(0, 1), 
            colorNA = 'transparent', 
            textNA = "") + 
    river_base + 
    tm_shape(subcatchments_rhine_union_2) + 
    tm_borders(col = 'black') + 
    tm_layout(legend.text.size = 1,
              legend.outside = F,
              legend.position = c('right', 'top')))}

dist_properties
map_values('suitability')
map_values('expected_distribution')
map_values('observed_distribution')
map_values('SD_expected_threat_continuous')
map_values('SD_expected_sum_presence_negative_threat')
map_values('SD_expected_sum_presence_negative_con')
map_values('SD_expected_sum_presence_negative_habitat')


# SUM OF SUITABILITIES
suitability_plot <-  tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'suitability', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(floor(min(all_dd_sum_sf$expected_distribution , na.rm = T)), 
                             ceiling(max(all_dd_sum_sf$expected_distribution , na.rm = T)), 
                             length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))

# EXPECTED DIVERSITY
expected_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'expected_distribution', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(floor(min(all_dd_sum_sf$expected_distribution, na.rm = T)), 
                             ceiling(max(all_dd_sum_sf$expected_distribution, na.rm = T)), 
                             length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))

# DARK DIVERSITY
dark_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'shadow_distribution_OminusE', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = '-Spectral', 
          n = 10, 
          contrast = c(0,1), 
          colorNA = 'transparent',
          midpoint = NA, 
          textNA = "", 
          breaks = signif(seq(0, 
                    max(all_dd_sum_sf$shadow_distribution_OminusE, na.rm = T), 
                    length.out = 5),2)) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))

# DARK DIVERSITY AS A WEIGHTED MEAN
dark_diversity_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'shadow_distribution_OminusE_weighted_mean', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = '-Spectral', 
          n = 10, 
          contrast = c(0,1), 
          colorNA = 'transparent',
          midpoint = NA, 
          textNA = "", 
          breaks = signif(seq(0, 
                              max(all_dd_sum_sf$shadow_distribution_OminusE_weighted_mean, na.rm = T), 
                              length.out = 5),2)) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))


# SUM OF NEGATIVE THREATS
sumAllThreatInE <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'shadow_distribution_sumAllThreatInE', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = 'Spectral', 
          n = 10, 
          contrast = c(0,1), 
          colorNA = 'transparent',
          midpoint = NA, 
          textNA = "") + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))


# SUM OF ANY NEGATIVE THREAT (E.G., THREAT RICHNESS)
sumAnyNegThreatInE <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'shadow_distribution_anyNegThreatInE', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = '-Spectral', 
          n = 10, 
          contrast = c(0,1), 
          colorNA = 'transparent',
          midpoint = NA, 
          textNA = "") + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))


# COMMUNITY COMPLETENESS
community_completeness_plot <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'community_completeness', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(floor(min(all_dd_sum_sf$community_completeness, na.rm = T)), 
                         ceiling(max(all_dd_sum_sf$community_completeness, na.rm = T)), 
                         length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))

# COMMUNITY COMPLETENESS THRESHOLD
completeness_threshold <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf %>% filter(!is.na(shap_suit_baseline))) + 
  tm_fill(col = 'community_completeness_threshold', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(50, 
                             ceiling(max(all_dd_sum_sf$community_completeness_threshold, na.rm = T)), 
                             length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('right', 'top'))



# make outputs fro the above plots
png(paste0(fig_dir, '/dark_diversity/combined_maps_v2.png'), 
    width = 1600*3, height = 800*3, res = 300,
    bg = 'transparent')
tmap_arrange(suitability_plot, 
             expected_diversity_plot, 
             dark_diversity_plot,
             community_completeness_plot,
             nrow = 2, ncol = 2)
dev.off()


#### MAKE BIG PAIRS PLOT OF FOUR PROPERTIES
pairs(st_drop_geometry(all_dd_sum_sf[, c("community_completeness", "shadow_distribution_OminusE_weighted_mean", 
                                   "expected_distribution", "observed_distribution")]))



## MAKE HISTOGRAMS OF PROPERTIES
all_dd_sum_sf$observed_distribution_fact <- round(all_dd_sum_sf$observed_distribution)
hist_suit <- all_dd_sum_sf %>% st_drop_geometry() %>% group_by(observed_distribution_fact) %>% count() %>% na.omit()
pdf(paste0(fig_dir, '/dark_diversity/hist_suitability.pdf'), width = 1.5, height = 1.5)
ggplot(data = hist_suit) + 
  geom_bar(aes(x = observed_distribution_fact, y = n, fill = observed_distribution_fact), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 0.75, 
        axis.text.y = element_blank(), 
        legend.position = 'none') + 
  xlab('') + ylab('') + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral")) + 
  scale_x_continuous(breaks = seq(0,8,2))
dev.off()

## MAKE HISTOGRAMS OF PROPERTIES
all_dd_sum_sf$expected_distribution_fact <- round(all_dd_sum_sf$expected_distribution)
hist_exp <- all_dd_sum_sf %>% st_drop_geometry() %>% group_by(expected_distribution_fact) %>% count() %>% na.omit()
pdf(paste0(fig_dir, '/dark_diversity/hist_expected_diversity.pdf'), width = 1.5, height = 1.5)
ggplot(data = hist_exp) + 
  geom_bar(aes(x = expected_distribution_fact, y = n, fill = expected_distribution_fact), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 0.75, 
        axis.text.y = element_blank(), 
        legend.position = 'none') + 
  xlab('') + ylab('') + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral")) + 
  scale_x_continuous(breaks = seq(0,8,2))
dev.off()



all_dd_sum_sf$shadow_distribution_OminusE_weighted_mean_fact <- round(all_dd_sum_sf$shadow_distribution_OminusE_weighted_mean,2)
hist_OminusE <- all_dd_sum_sf %>% st_drop_geometry() %>% group_by(shadow_distribution_OminusE_weighted_mean_fact) %>% count() %>% na.omit()
pdf(paste0(fig_dir, '/dark_diversity/hist_dark_diversity.pdf'), width = 1.5, height = 1.5)
ggplot(data = hist_OminusE) + 
  geom_bar(aes(x = shadow_distribution_OminusE_weighted_mean_fact, y = n, fill = shadow_distribution_OminusE_weighted_mean_fact), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 0.75, 
        axis.text.y = element_blank(), 
        legend.position = 'none') + 
  xlab('') + ylab('') + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral"))) 
dev.off()



all_dd_sum_sf$community_completeness_fact <- round(all_dd_sum_sf$community_completeness/100,1)
hist_complete <- all_dd_sum_sf %>% st_drop_geometry() %>% group_by(community_completeness_fact) %>% count() %>% na.omit()
hist_complete$community_completeness_fact <- hist_complete$community_completeness_fact*100
pdf(paste0(fig_dir, '/dark_diversity/hist_completeness.pdf'), width = 1.5, height = 1.5)
ggplot(data = hist_complete) + 
  geom_bar(aes(x = community_completeness_fact, y = n, fill = community_completeness_fact), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(), 
        aspect.ratio = 0.75, 
        axis.text.y = element_blank(), 
        legend.position = 'none') + 
  xlab('') + ylab('') + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10, "Spectral")) + 
  scale_x_continuous(breaks = hist_complete$community_completeness_fact[c(1,3,5,7)])
dev.off()


#### Sections of completeness plots ----

# summarise the distribution properties per subcatchment
dist_properties <- c('suitability', 
                     'expected_distribution', 
                     c(hab_factors, con_factors, nn_factors))
all_dd_sum <- all_dd_bind %>% 
  select(TEILEZGNR, species_name, 
         any_of(dist_properties)) %>% 
  filter(species_name != 'Oncorhynchus mykiss') %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  summarise_at(dist_properties, 
               sum, 
               na.rm = T) %>% 
  left_join(., all_dd_bind %>% select(TEILEZGNR) %>% unique()) %>% 
  # calculate community completeness here
  mutate(community_completeness_threshold = (.$suitability / expected_distribution)*100) %>% 
  mutate(community_completeness_threshold = ifelse(.$community_completeness_threshold > 100, 100, .$community_completeness_threshold)) %>% 
  # filter out areas outside the ecological niche of all species
  filter(!is.na(expected_distribution), 
         expected_distribution != 0) # remove areas where we do not expect species

all_dd_sum$habitat_sum <- rowMeans(all_dd_sum[,hab_factors],na.rm = T)
all_dd_sum$nn_sum <- rowMeans(all_dd_sum[,nn_factors],na.rm = T)

# plot the summarised distributions
all_dd_sum_sf <- all_dd_sum %>% st_as_sf()

### MAP OF COMMUNITY COMPLETENESS THRESHOLD
completeness_threshold <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'community_completeness_threshold', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(50,#floor(min(all_dd_sum_sf$community_completeness_threshold, na.rm = T)), 
                             ceiling(max(all_dd_sum_sf$community_completeness_threshold, na.rm = T)), 
                             length.out = 5)), 
          values = c(0, 0.6, 0.65,  0.7, 0.8, 0.9, 0.95,1)) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('left', 'top'))


### MAP OF COMPLETENESS 
png(paste0(fig_dir, '/dark_diversity/completeness_maps.png'), 
    width = 1600, height = 800, res = 300,
    bg = 'transparent')
completeness_threshold
dev.off()


### MAP SUMMED EXPECTED DISTRIBUTIONS AS EXPECTED DIVERSITY
expected_diversity_map <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'expected_distribution', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(floor(min(all_dd_sum_sf$expected_distribution, na.rm = T)), 
                             ceiling(max(all_dd_sum_sf$expected_distribution, na.rm = T)), 
                             length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('left', 'top'))

png(paste0(fig_dir, '/dark_diversity/expected_diversity_maps.png'), 
    width = 1600, height = 800, res = 300,
    bg = 'transparent')
expected_diversity_map
dev.off()


### MAP SUMMED SUITABILITIES INSIDE ECOLOGICAL NICHE
richness_map <- tm_shape(subcatchments_rhine_union_2) + 
  tm_polygons(col = 'gray90') + 
  tm_shape(all_dd_sum_sf) + 
  tm_fill(col = 'suitability', 
          style = 'cont', 
          title = '', 
          legend.is.portrait = T,
          legend.reverse = T,
          palette = "Spectral", 
          n = 10, 
          contrast = c(0, 1), 
          colorNA = 'transparent', 
          textNA = "", 
          breaks = round(seq(floor(min(all_dd_sum_sf$suitability, na.rm = T)), 
                             ceiling(max(all_dd_sum_sf$suitability, na.rm = T)), 
                             length.out = 5))) + 
  river_base + 
  tm_shape(subcatchments_rhine_union_2) + 
  tm_borders(col = 'black') + 
  tm_layout(legend.text.size = 1,
            legend.outside = F,
            legend.position = c('left', 'top'))

png(paste0(fig_dir, '/dark_diversity/richness_maps.png'), 
    width = 1600, height = 800, res = 300,
    bg = 'transparent')
richness_map
dev.off()


### MAKE PLOT OF COMPLETNESS VS. SHAPLEY VALUES
all_dd_long <- all_dd_sum %>% 
  st_drop_geometry() %>% 
  select(-Shape) %>% 
  pivot_longer(., cols = c(habitat_sum, con_factors, nn_sum))

# input data for plot
shapley_vs_completeness <- ggplot(data = all_dd_long %>% 
         filter(community_completeness_threshold != 100)) +
  # points and smoothers
  geom_vline(aes(xintercept = 0), lty = 2) + 
  geom_point(aes(x = value, y = community_completeness_threshold, col = name), 
             alpha = 0.2, 
             pch = 19, 
             stroke = 0, 
             size = 2) + 
  stat_smooth(aes(x = value, y = community_completeness_threshold, group = name, 
                  col = name), 
              fill = 'black', 
              alpha = 0.5, 
              linewidth = 0.1, 
              show.legend = FALSE) + 
  # aesthetics
  theme_bw() + 
  theme(legend.position = c(0.75, 0.25), 
        panel.grid = element_blank(), 
        aspect.ratio = 1, 
        panel.border = element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.background = element_blank()) + 
  # labels
  ylab('Community completeness (%)') + 
  xlab('Shapley value') + 
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5)) + 
    scale_colour_discrete(labels = c('habitat factors', 'connectivity factors', 'natural niche factors'), 
                          name = '') + 
  guides(color = guide_legend(override.aes = list(size = 5, 
                                                  alpha = 1)))

# print figure to pdf
pdf(paste0(fig_dir, '/dark_diversity/shapley_vs_completeness.pdf'), width = 3, height = 3)
print(shapley_vs_completeness)
dev.off()



## MAKE PLOT OF EXPECTED DIVERSITY AGAINST COMPLETENESS
library(ggridges)
pdf(paste0(fig_dir, '/dark_diversity/expected_vs_completeness.pdf'), width = 3, height = 3)
ggplot(data = all_dd_sum %>% filter(expected_distribution >= 1)) +
  geom_point(aes(x = community_completeness_threshold, 
                  y = as.character(round(expected_distribution)), 
                 col = community_completeness_threshold), 
             position = position_jitter(height = 0.1, width = 0)) + 
  geom_violin(aes(x = community_completeness_threshold, 
                  y = as.character(round(expected_distribution))), 
              scale = 'width', 
              trim = T, 
              fill = 'transparent') + 
  # aesthetics
  theme_bw() + 
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        aspect.ratio = 1, 
        panel.border = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_line()) + 
  scale_colour_gradientn(
    colours = RColorBrewer::brewer.pal(10, "Spectral"), 
    values = c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) + 
  xlab('Community completeness (%)') + 
  ylab('Expected species richness')
dev.off()
  
  

# ### split histograms based on river size
# # read in rivers
# rivers_v2 <- paste0(dd_env, "CH_HYDRO.shp")
# rivers_v2 <- st_read(rivers_v2)
# rivers_large <- rivers_v2 %>%
#   filter(STRAHLE > 2) %>% 
#   select(geometry)
# 
# # get intersection of rivers
# all_dd_large_rivers <- st_join(rivers_large, all_dd_sum_sf) %>% 
#   st_drop_geometry() %>% 
#   unique() %>% 
#   na.omit()
# 
# all_dd_small_rivers <- all_dd_sum_sf %>% 
#   filter(!TEILEZGNR %in% unique(rivers_large$TEILEZGNR)) %>% 
#   st_drop_geometry() %>% 
#   unique() %>% 
#   na.omit()
# 
# # combine small and large rivers
# all_small_large <- bind_rows(all_dd_large_rivers %>% mutate(size = 'large'), 
#           all_dd_small_rivers %>% mutate(size = 'small')) %>% 
#   mutate(round_ed = as.character(round(expected_distribution))) %>% 
#   filter(round_ed != '0')
# 
# # summarise mean values
# all_small_large %>% 
#   group_by(round_ed, size) %>% 
#   summarise_at(., 'community_completeness_threshold', 
#                list(quantile_5 = function(x) round(quantile(x, 0.05, na.rm = T)), 
#                     mean = function(x) round(mean(x))))
# 
# ggplot() +
#   geom_jitter(data = all_small_large,
#              aes(x = community_completeness_threshold, 
#                  y = round_ed, 
#                  col = community_completeness_threshold, 
#                  group = paste0(size, round_ed)), 
#              position = position_jitterdodge(jitter.width = 0.1, 
#                                              jitter.height = 0, 
#                                              dodge.width = 0.75)) + 
#   geom_boxplot(data = all_small_large, 
#                aes(x = community_completeness_threshold, 
#                   y = round_ed, 
#                   group = paste0(size, round_ed), 
#                   fill = size), 
#               trim = T, 
#               outlier.colour =  'transparent') + 
#   # aesthetics
#   theme_bw() + 
#   theme(legend.position = 'none', 
#         panel.grid = element_blank(), 
#         aspect.ratio = 1, 
#         panel.border = element_blank(), 
#         axis.ticks = element_blank(), 
#         axis.line = element_line()) + 
#   scale_colour_gradientn(
#     colours = RColorBrewer::brewer.pal(10, "Spectral"), 
#     values = c(0, 0.6, 0.65,  0.7, 0.8, 0.9, 0.95,1)) + 
#   xlab('Community completeness (%)') + 
#   ylab('Expected species richness') + 
#   xlim(c(NA, 101))



#### 9. Investigate properties of complete and incomplete catchments ----

# see what species are present in the high completeness regions
subcatch_complete100 <- all_dd_sum %>% 
  filter(community_completeness_threshold == 100) %>% 
  pull(TEILEZGNR) %>% 
  unique() 
  
subcatch_complete_others <- all_dd_sum %>% 
  filter(community_completeness_threshold != 100, 
         expected_distribution >= 1) %>% 
  pull(TEILEZGNR) %>% 
  unique()

# what are the properties of these streams
cbind(all_env_subcatchments %>% 
  filter(TEILEZGNR %in% subcatch_complete100) %>% 
  summarise_all(., function(x) round(mean(x),2)) %>% t,
  
  all_env_subcatchments %>% 
  filter(TEILEZGNR %in% subcatch_complete_others) %>% 
  summarise_all(., function(x) round(mean(x),2)) %>% 
    t) %>% 
  view()


## take the full species dataset, subset to those 
# subcatchments with 100% completeness
# and table up modelling presence of species
sp_id_in_100_present <- all_dd_bind %>% 
  filter(TEILEZGNR %in% subcatch_complete100,
         !is.na(presence)) %>% 
  select(TEILEZGNR, species_name) %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  summarize(type = paste(sort(unique(species_name)),collapse=", "), 
            species_n = length(unique(species_name)))

table(sp_id_in_100$species_n)
sort(table(sp_id_in_100$type)) / length(sp_id_in_100$type)



## take the full species dataset, subset to those 
# subcatchments with 100% completeness
# and table up the expected distribution
sp_id_in_100_expected <- all_dd_bind %>% 
  filter(TEILEZGNR %in% subcatch_complete100,
         !is.na(expected_distribution)) %>% 
  select(TEILEZGNR, species_name) %>% 
  st_drop_geometry() %>% 
  group_by(TEILEZGNR) %>% 
  summarize(type = paste(sort(unique(species_name)),collapse=", "), 
            species_n = length(unique(species_name)))

table(sp_id_in_100_expected$species_n)
sort(table(sp_id_in_100_expected$type)) / length(sp_id_in_100_expected$type)

  
  

