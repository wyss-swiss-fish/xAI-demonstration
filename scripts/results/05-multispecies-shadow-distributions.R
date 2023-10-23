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
fig_dir <- "figures/ubelix_SDM_RF_APRIL_V1_02/"

# get run to make figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"


# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

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
sdm_dirs <- list.files(paste0("D:/sdm-pipeline/sdm-run/", RUN_SDM), full.names = T)
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
source('scripts/functions/generate_SD_outputs.R')

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

#### 5. Boxplots of shadow distribution summaries across species ----

# colour palette for species
library(seecolor)
pal <- c('#B436FF', '#E83197', '#FF6242', '#E88C31', '#FFCC35', 
         '#E2EB9E', '#65EBA5', '#4DE3FF', '#6772FF' )
print_color(pal)

# create plot of summary properties for positive contributions of niche effects
pos_boxplot <- ggplot(data = summary_shadow %>% 
         filter(grepl('% all subcatchments with a positive contribution', property)) %>% 
         mutate(property = recode(.$property, 
                                  '% all subcatchments with a positive contribution of ecoF_discharge_max_log10_SHAP' = 'discharge', 
                                  '% all subcatchments with a positive contribution of ecoF_slope_min_log10_SHAP' = 'slope', 
                                  '% all subcatchments with a positive contribution of stars_t_mn_m_c_SHAP' = 'min. temp',
                                  '% all subcatchments with a positive contribution of stars_t_mx_m_c_SHAP' = 'max. temp.',
                                  '% all subcatchments with a positive contribution of ecoF_flow_velocity_mean_SHAP' = 'flow velocity',
                                  '% all subcatchments with a positive contribution of local_dis2lake_SHAP' = 'distance to lake'))) + 
  geom_boxplot(aes(x = property, y = value), outlier.colour = NA, width = 0.25, position= position_nudge(x=-0.2)) + 
  geom_jitter(aes(x = property, y = value, fill = species), size = 3, pch = 21, 
              position = position_nudge(x=0.2)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank(),
        plot.title = element_text(size = 10), 
        legend.position = 'none') + 
  ylab('% subcatchments') + 
  xlab(NULL) + 
  scale_fill_manual(values = pal)
pos_boxplot

# create plot of summary properties for negative contributions of threat effects
neg_boxplot <- ggplot(data = summary_shadow %>% 
         filter(grepl('% subcatchments inside niche with a negative contribution of ', property)) %>% 
         mutate(property = recode(.$property, 
                                  '% subcatchments inside niche with a negative contribution of local_flood_SHAP' = 'floodplain', 
                                  '% subcatchments inside niche with a negative contribution of local_imd_log10_ele_residual_SHAP' = 'urbanisation', 
                                  '% subcatchments inside niche with a negative contribution of ecoF_eco_mean_ele_residual_SHAP' = 'morph. mod.',
                                  '% subcatchments inside niche with a negative contribution of local_asym_cl_log10_SHAP' = 'connectivity', 
                                  '% subcatchments inside niche with a negative contribution of local_wet_SHAP' = 'wetland'))) + 
  geom_boxplot(aes(x = property, y = value), outlier.colour = NA, width = 0.25, position= position_nudge(x=-0.2)) + 
  geom_jitter(aes(x = property, y = value, fill = species), size = 3, pch = 21, 
              position = position_nudge(x=0.2)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank(),
        plot.title = element_text(size = 10), 
        legend.position = 'none') + 
  ylab('% expected distribution') + 
  xlab(NULL) + 
  scale_fill_manual(values=pal)
neg_boxplot

# percentages of subcatchments with 'n' number of threats
threat_boxplot <- ggplot(data = summary_shadow %>% 
         filter(grepl('% of subcatchments inside niche with', property)) %>% 
         mutate(property = recode(.$property, 
                                  '% of subcatchments inside niche with 0 threats negative' = '0', 
                                  '% of subcatchments inside niche with 1 threats negative' = '1', 
                                  '% of subcatchments inside niche with 2 threats negative' = '2', 
                                  '% of subcatchments inside niche with 3 threats negative' = '3', 
                                  '% of subcatchments inside niche with 4 threats negative' = '4', 
                                  '% of subcatchments inside niche with 5 threats negative' = '5'))) + 
  geom_boxplot(aes(x = property, y = value), outlier.colour = NA, width = 0.25, position= position_nudge(x=-0.2)) + 
  geom_jitter(aes(x = property, y = value, fill = species), size = 3, pch = 21, 
              position = position_nudge(x=0.2)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank(),
        plot.title = element_text(size = 10), 
        legend.position = 'none') + 
  ylab('% expected distribution') + 
  xlab(NULL) + 
  scale_fill_manual(values = pal)
threat_boxplot

png(filename = paste0(fig_dir, 'shadow_dist_summaries/perc_across_species.png'),
    res = 300, width = 2500, height = 750)
cowplot::plot_grid(pos_boxplot, neg_boxplot, threat_boxplot, align = 'h', nrow = 1)
dev.off()

png(filename = paste0(fig_dir, 'shadow_dist_summaries/species_legend.png'),
    res = 300, width = 1000, height = 1000)
plot(cowplot::get_legend(ggplot(data = summary_shadow) + 
                            geom_point(aes(x=property, y=value, fill = species), pch = 21, size = 5) + 
                           theme(
                                 legend.key = element_rect(colour = 'white', fill = 'white'),
                                 legend.text = element_text(size = 15, face = 'italic')
                                 ) +
                            scale_fill_manual(values = pal)))
dev.off()

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

# run the quantitative shadow distribution estimation
# which creates additional columns on our spatial objects 
# that summarise the expected, observed and shadow distributions
source('scripts/functions/generate_quant_SD.R')

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


#### 7. Plot each distribution type across all species ----

tmap_mode('plot')

river_lake_tm <- tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F)

for(method in 1:length(all_dd_method)){

# create output directory
shad_plot_dir <- paste0(fig_dir, '/shadow_dist_summaries/shadow_distribution_maps', '_method_', method, '/')
dir.create(shad_plot_dir)

# subset to method in loop
all_dd <- all_dd_method[[method]]

# create maps across all species
lapply(1:length(all_dd), function(x){
  
  sp_dd <- all_dd[[x]]
  
  pdf(paste0(shad_plot_dir, unique(sp_dd$species_name), '.pdf'),
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

}



#### 8. Aggregate together shadow distribution values ----

lapply(1:4, function(x) dir.create(paste0(fig_dir, '/dark_diversity/', 'method_', x), recursive = T))

for(method in 1:length(all_dd_method)){
  
# bind together all outputs
all_dd_bind <- bind_rows(all_dd_method[[method]])

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

png(paste0(fig_dir, '/dark_diversity/', 'method_', method, '/', 'all_maps_method_', method, '.png'),
    width = 4000, height = 2250, res = 300,
    bg = "transparent")
print(tmap_arrange(map_values('observed_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             map_values('expected_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             map_values('expected_distribution_mean', breaks = range, palette = 'Spectral', legend.reverse = T),
             
             mean_completeness_plot,
             min_completeness_plot, 
             sd_completeness_plot,
             
             map_values('SD_expected_sum_presence_negative_con_mean', breaks = round(seq(0, 1, length.out = 4),2), palette = '-Spectral', legend.reverse = T),
             map_values('SD_expected_sum_presence_negative_habitat_mean', breaks = round(seq(0, 4, length.out = 4), 2), palette = '-Spectral', legend.reverse = T), 
             map_values('SD_expected_sum_presence_negative_habitat_mean', breaks = round(seq(0, 4, length.out = 4), 2), palette = '-Spectral', legend.reverse = T),
            
             nrow = 3, ncol = 3))
dev.off()


#### 9. Data summaries of each shadow distribution property ----


sink(paste0(fig_dir, '/dark_diversity/', 'method_', method, '/', 'summaries_method_', method, '.txt'))

print('MEAN SUITABILITY OF THE OBSERVED DISTRIBTUION')
## Mean suitability
# Summary of mean observed distribution
print(round(summary(all_dd_sum_sf$observed_distribution_mean), 2))
mean_suit <- mean(all_dd_sum_sf$observed_distribution_mean, na.rm = T)

print('MEAN SUITABILITY OF THE EXPECTED DISTRIBTUION')
# Summary of mean expected distribution
print(round(summary(all_dd_sum_sf$expected_distribution_mean), 2))
mean_exp <- mean(all_dd_sum_sf$expected_distribution_mean, na.rm = T)

print('MEAN PERCENTAGE REDUCTION IN SUITABILITY BETWEEN OBSERVED AND EXPECTED')
# Mean reduction in suitability, calculate directly from species-metric
print(1 - mean(all_dd_sum_sf$SD_OratioE_mean)) 

print('T-TEST OF REDUCTION IN SUITABILITIES BETWEEN OBSERVED AND EXPECTED')
print(t.test(all_dd_sum_sf$observed_distribution_mean, all_dd_sum_sf$expected_distribution_mean))


print('MINIMUM SUITABILITY OF OBSERVED DISTRIBUTION WITHIN SET OF SPECIES')
## Minimum suitability
# Summary of minimum observed distribution
print(round(summary(all_dd_sum_sf$observed_distribution_min), 2))
mean_suit <- mean(all_dd_sum_sf$observed_distribution_min,na.rm = T)

print('MINIMUM SUITABILITY OF EXPECTED DISTRIBUTION WITHIN SET OF SPECIES')
# Summary of minimum expected distribution
print(round(summary(all_dd_sum_sf$expected_distribution_min), 2))
mean_exp <- mean(mean(all_dd_sum_sf$expected_distribution_min,na.rm = T))

print('PERCENTAGE DIFFERENCE BETWEEN MINIMUM EXPECTED AND OBSERVED SUITABILITY')
## Percentage difference between observed and expected
# Mean percentage reduction in suitability comparing observed to expected 
print(1-mean(all_dd_sum_sf$SD_OratioE_min)) # 0.2911247% reduction in suitability 
t_test_exp_obs <- t.test(all_dd_sum_sf$observed_distribution_min, all_dd_sum_sf$expected_distribution_min)
print(t_test_exp_obs)


print('SUMMARY OF MEAN ACROSS SPECIES PERCENTAGE DIFFERENCE BETWEEN OBSERVED AND EXPECTED')
# Summary of mean percentage difference between observed and expected suitability
print(summary(all_dd_sum_sf$SD_OpercentOfE_mean))

print('SUMMARY OF MINIMUM ACROSS SPECIES PERCENTAGE DIFFERENCE BETWEEN OBSERVED AND EXPECTED')
# Summary of minmum percentage difference between observed and expected suitability
print(summary(all_dd_sum_sf$SD_OpercentOfE_min))
mean_of_min_shad_per_species <- signif(mean(all_dd_sum_sf$SD_OpercentOfE_min, na.rm = T),2)

print('PERCENTAGE DIFFERENCE IN WORST 10% OF CATCHMENTS')
## Percentage difference in worst 10% of catchments 
# What is the mean reduction in suitability in the lowest 10th quantile 
worst_10perc <- mean(all_dd_sum_sf$SD_OpercentOfE_mean[all_dd_sum_sf$SD_OpercentOfE_mean<quantile(all_dd_sum_sf$SD_OpercentOfE_mean, 0.1)])
print(worst_10perc)


print('LARGEST REDUCTION IN HABITAT SUITABILITY (MINIMUM OF MEAN REDUCTION)')
## Maximum reduction in habitat suitability
# across-catchment minimum of mean reduction
print(min(all_dd_sum_sf$SD_OratioE_mean))
print('LARGEST REDUCTION IN HABITAT SUITABILITY (MINIMUM OF MINIMUM REDUCTION)')
# across-catchment minimum of minimum reduction
print(min(all_dd_sum_sf$SD_OratioE_min))


## Comparative correlations between threat numbers and percentage reductions
print('spearmans rank between shadow distribution and threat number - mean shadow in community')
# spearmans rank between shadow distribution and threat number
cor_con_shadmean <- cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_con_mean, all_dd_sum_sf$SD_OratioE_mean, method = 'spearman')
print(cor_con_shadmean)
cor_hab_shadmean <- cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_habitat_mean, all_dd_sum_sf$SD_OratioE_mean, method = 'spearman')
print(cor_hab_shadmean)

# spearmans rank between shadow distribution and threat number
print('spearmans rank between shadow distribution and threat number - minimum shadow in community')
cor_con_shadmin <- cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_con_mean, all_dd_sum_sf$SD_OratioE_min, method = 'spearman')
print(cor_con_shadmin)
cor_hab_shadmin <- cor.test(all_dd_sum_sf$SD_expected_sum_presence_negative_habitat_mean, all_dd_sum_sf$SD_OratioE_min, method = 'spearman')
print(cor_hab_shadmin)

sink()



### save object as named list that contains information for summaries across method
grab_stats <- function(x, cor_or_stat){
  glance_t <- broom::glance(x)
  p <- x$p.value
  if(cor_or_stat == 'stat'){  out <- paste0(signif(x$statistic,2),'; ', ifelse(p > 0.05, 'p>0.05',
                                                                    ifelse(p < 0.001, 'p<0.001', 
                                                                           ifelse(p<0.01, 'p<0.01', 
                                                                                  ifelse(p<0.05, 'p<0.05', 'error')))))
  }
  if(cor_or_stat == 'cor'){    out <- paste0(signif(x$estimate,2),'; ', ifelse(p > 0.05, 'p>0.05',
                                                     ifelse(p < 0.001, 'p<0.001', 
                                                            ifelse(p<0.01, 'p<0.01', 
                                                                   ifelse(p<0.05, 'p<0.05', 'error')))))
    
  }
  return(out)
}


### get the objects and values that are cited directly in the manuscript text for saving externally
summary_sd_per_method <- list(
  `mean suitability of expected distribution` = signif(mean_exp, 2),
  `mean suitability of observed distribution` = signif(mean_suit,2),
  `t-test summary between expected and observed` = grab_stats(t_test_exp_obs, cor_or_stat = 'stat'), 
  `perc. reduction in observed vs. expected` = signif(1 - mean(all_dd_sum_sf$SD_OratioE_mean, na.rm = T), 2),
  `suitability in lowest 10% catchments` = signif(worst_10perc,2), 
  `minimum species shadow per catchment` = mean_of_min_shad_per_species, 
  `spearman: connectivity vs expected mean` = grab_stats(cor_con_shadmean, cor_or_stat = 'cor'), 
  `spearman: habitat vs expected mean` = grab_stats(cor_hab_shadmean, cor_or_stat = 'cor'), 
  `spearman: connectivity vs expected min` = grab_stats(cor_con_shadmin, cor_or_stat = 'cor'), 
  `spearman: habitat vs expected mean` = grab_stats(cor_hab_shadmin, cor_or_stat = 'cor')
)
saveRDS(summary_sd_per_method, 
        paste0(fig_dir, '/dark_diversity/', 'method_', method, '/summary_statistics.rds'))


### plot biplots of expected vs. observed distributions
png(paste0(fig_dir, '/dark_diversity/', 'method_', method, '/biplot_threats.png'),
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

}



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


#### write results to text ----

m1_shad <- readRDS(paste0(fig_dir, '/dark_diversity/', 'method_', 1, '/summary_statistics.rds'))
m2_shad <- readRDS(paste0(fig_dir, '/dark_diversity/', 'method_', 2, '/summary_statistics.rds'))
m3_shad <- readRDS(paste0(fig_dir, '/dark_diversity/', 'method_', 3, '/summary_statistics.rds'))

save(summary_shadow_mean, 
     m1_shad, 
     m2_shad,
     m3_shad,
     file = 'data/05-results-text.RData')

