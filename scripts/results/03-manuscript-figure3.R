### Figure 3. Look at ways how disconnection occurs between species ----

### Pick interesting examples of pairs of species and variables to investigate ----

### Barbu or Abramis and Cyprinus carpio for ecomorphology (divergent responses)
### Cottus gobio and Anguilla anguilla for connectivity (similar responses but different magnitude)
### Little response.. 

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

#### 2. Set directories for loading and saving objects ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures-march2023-v2"

# get run to mak figures for
RUN <- "ubelix_test_PARA_RF_v3"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- c("Barbus barbus", "Cyprinus carpio", "Cottus gobio", "Anguilla anguilla", "Lampetra planeri", "Telestes souffia")
# sp_list <- unique(records_table$species_name)


# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
                        full.names = T
)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]

# create shapley folders per species
shap_po <- paste0(shap_dirs, "/shapley_rf_po.RDS")
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")

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

#### Read in spatial data and assign simpler names to variables ----

# load spatial objects
source('scripts/results/00-load-spatial.R')

# set vector of focal variables for assigning new names
vars <- c('ecoF_discharge_max_log10', 
          'ecoF_slope_min_log10', 
          'ecoF_flow_velocity_mean', 
          'stars_t_mx_m_c', 
          'stars_t_mn_m_c',
          'local_asym_cl_log10',
          'ecoF_eco_mean', 
          'stars_lu_crp_p_e_ele_residual',
          'local_tcd_0.1_log10', 
          'local_imd_log10_ele_residual',
          'stars_lud_m_m_e_ele_residual',
          'local_wet',
          'local_flood',
          'stars_n_ch_m_e_ele_residual',
          'stars_p_ch_m_e_ele_residual',
          'stars_isct_m_e_ele_residual')
# set vars
vars_shap <- paste0(vars, "_SHAP")

# rename variables
vars_renamed = c('discharge', 'slope', 'flow velocity', 'temperature max', 'temperature min', 
                 'connectivity', 'ecomorphology', 'cropland', 'tree cover', 'urbanisation',
                 'livestock', 'wetland', 'floodplains', 'nitrogen', 'phosphorous', 
                 'insecticide')

vars_renamed <- cbind(vars_shap, vars_renamed)

#### Make rasters of shapley values for faster future use ----

for (i in 1:length(sp_list)) {
  # check if PO files exist
  if (file.exists(shap_po[i])) {
    # read in the RDS
    sp_shap_i <- readRDS(shap_po[i])
    
    # join shapley values to the teilenzugsgebeit
    sp_shap_i <- left_join(
      subcatchments_final %>%
        select(TEILEZGNR),
      sp_shap_i
    ) %>%
      filter(!is.na(suitability))
    
    # turn to vector in terra for rasterization
    sp_shap_vect_i <- vect(sp_shap_i)
    
    # rasterize
    shap_rast_i <- lapply(vars_shap, function(x) {
      print(x)
      if (x %in% names(sp_shap_vect_i)) {
        rasterize(
          x = sp_shap_vect_i,
          y = base_rast,
          field = x,
          touches = T,
          fun = function(x) mean(x, na.rm = T)
        )
      }
    })
    
    shap_rast_i <- rast(shap_rast_i)
    
    writeRaster(shap_rast_i, filename = paste0(shap_dirs[i], "/shap_raster_po.TIF"), overwrite = T)
  }
  
  # check if PA files exist
  if (file.exists(shap_pa[i])) {
    # read in the RDS
    sp_shap_i <- readRDS(shap_pa[i])
    
    # join shapley values to the teilenzugsgebeit
    sp_shap_i <- left_join(
      subcatchments_final %>%
        select(TEILEZGNR),
      sp_shap_i
    ) %>%
      filter(!is.na(suitability))
    
    # turn to vector in terra for rasterization
    sp_shap_vect_i <- vect(sp_shap_i)
    
    # rasterize
    shap_rast_i <- lapply(vars_shap, function(x) {
      print(x)
      if (x %in% names(sp_shap_vect_i)) {
        rasterize(
          x = sp_shap_vect_i,
          y = base_rast,
          field = x,
          touches = T,
          fun = function(x) mean(x, na.rm = T)
        )
      }
    })
    
    shap_rast_i <- rast(shap_rast_i)
    
    writeRaster(shap_rast_i, filename = paste0(shap_dirs[i], "/shap_raster_pa.TIF"), overwrite = T)
    
  }
}

#### Example 1. Contrast response curves for ecomorphology across species -----

# read in rasters for focal variables
shap_rast_Bb <- rast(paste0(shap_dirs[grepl('Barbus', shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_Cc <- rast(paste0(shap_dirs[grepl('Cyprinus', shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_Bb <- readRDS(shap_pa[grepl('Barbus', shap_pa)])
shap_Cc <- readRDS(shap_pa[grepl('Cyprinus', shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_Bb <- left_join(shap_Bb %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_Cc <- left_join(shap_Cc %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_Bb_Cc <- bind_rows(shap_Bb, shap_Cc) %>% 
  select(species_name, TEILEZGNR, matches('ecoF_eco_mean'))

ggplot(data = shap_Bb_Cc, aes(x = ecoF_eco_mean, y = ecoF_eco_mean_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = T, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5) + 
  ylab('Shapley value') + 
  xlab('Environmental value') + 
  scale_colour_manual(values = c('gray50', 'gray75'))


# make wide for plot of shapley values for each species against eachother
shap_Bb_Cc_wide <- shap_Bb_Cc %>% 
  pivot_wider(names_from = 'species_name', values_from = 'ecoF_eco_mean_SHAP')

shap_Bb_Cc_wide$col <- ifelse(shap_Bb_Cc_wide$`Barbus barbus` >= shap_Bb_Cc_wide$`Cyprinus carpio` & 
                                shap_Bb_Cc_wide$`Barbus barbus` > 0 & shap_Bb_Cc_wide$`Cyprinus carpio` < 0,
                              1, 
                              ifelse(shap_Bb_Cc_wide$`Barbus barbus` <= shap_Bb_Cc_wide$`Cyprinus carpio` & 
                                       shap_Bb_Cc_wide$`Cyprinus carpio` > 0 & shap_Bb_Cc_wide$`Barbus barbus` < 0, 
                                     2, 3))
ggplot(data = shap_Bb_Cc_wide, 
       aes(x = `Barbus barbus`, y = `Cyprinus carpio`, col = as.factor(col))) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) + 
  scale_colour_manual(values = c('#bd0f06', '#2200c9', 'gray90')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', aspect.ratio = 0.5) + 
  xlab('Shapley values: Barbus barbus') + 
  ylab('Shapley values: Cyprinus carpio')
  

tm_shape(shap_rast_Bb['ecoF_eco_mean_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = F,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Barbus barbus')

tm_shape(shap_rast_Cc['ecoF_eco_mean_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = T,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Cyprinus carpio')


# read in presence predictions for both species
pres_Bb <- rast(sp_raster_pres_pa[grepl('Barbus', sp_raster_pres_pa)])
pres_Cc <- rast(sp_raster_pres_pa[grepl('Cyprinus', sp_raster_pres_pa)])
pres_Bb_Cc <- mosaic(pres_Bb, pres_Cc)

shap_div_Bb <- shap_rast_Bb['ecoF_eco_mean_SHAP'] > shap_rast_Cc['ecoF_eco_mean_SHAP'] & shap_rast_Bb['ecoF_eco_mean_SHAP'] > 0
shap_div_Cc <- shap_rast_Cc['ecoF_eco_mean_SHAP'] > shap_rast_Bb['ecoF_eco_mean_SHAP'] & shap_rast_Cc['ecoF_eco_mean_SHAP'] > 0

shap_div_Bb[shap_div_Bb[]==F] <- NA
shap_div_Bb[shap_div_Bb[]==T] <- 1
shap_div_Cc[shap_div_Cc[]==F] <- NA
shap_div_Cc[shap_div_Cc[]==T] <- 2

shap_div_Bb_Cc <- mosaic(shap_div_Cc, shap_div_Bb)
pres_Bb_Cc <- resample(pres_Bb_Cc, shap_div_Bb_Cc)
shap_div_Bb_Cc[is.na(pres_Bb_Cc)] <- NA

tm_shape(shap_div_Bb_Cc) +
  tm_raster(n = 2,
            palette = c('#bd0f06', '#2200c9'), 
            labels = c('Ecomorphology favours B. barbus', 'Ecomorphology favours C. carpio'), 
            title = "Contrasted Ecomorphology effects") + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent")

# plot the environmental data
tm_shape(env_data['ecoF_eco_mean']) + 
  tm_raster(style = 'cont', 
            palette = '-viridis', 
            legend.is.portrait = F, 
            title = 'Ecomorphological modification') + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent")


#### Example 2. Contrast response curves for connnectivity across species ----

# read in rasters for focal variables
shap_rast_Aa <- rast(paste0(shap_dirs[grepl('Anguilla', shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_Cg <- rast(paste0(shap_dirs[grepl('Cottus gobio', shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_Aa <- readRDS(shap_pa[grepl('Anguilla', shap_pa)])
shap_Cg <- readRDS(shap_pa[grepl('Cottus gobio', shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_Aa <- left_join(shap_Aa %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_Cg <- left_join(shap_Cg %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_Aa_Cg <- bind_rows(shap_Aa, shap_Cg) %>% 
  select(species_name, TEILEZGNR, matches('local_asym_cl_log10'))

ggplot(data = shap_Aa_Cg, aes(x = local_asym_cl_log10, y = local_asym_cl_log10_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black',se = F, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4)) + 
  geom_hline(aes(yintercept = 0), lty = 2) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5) + 
  ylab('Shapley value') + 
  xlab('Environmental value')

# make wide for plot of shapley values for each species against eachother
shap_Aa_Cg_wide <- shap_Aa_Cg %>% 
  pivot_wider(names_from = 'species_name', values_from = 'local_asym_cl_log10_SHAP')

shap_Aa_Cg_wide$col <- ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` > 0 & shap_Aa_Cg_wide$`Cottus gobio` > 0 & shap_Aa_Cg_wide$`Anguilla anguilla` >= shap_Aa_Cg_wide$`Cottus gobio`,
                              1, 
                              ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` > 0 & shap_Aa_Cg_wide$`Cottus gobio` > 0 & shap_Aa_Cg_wide$`Anguilla anguilla` <= shap_Aa_Cg_wide$`Cottus gobio`, 
                                     2, 
                                     ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` < 0 & shap_Aa_Cg_wide$`Cottus gobio` < 0 & shap_Aa_Cg_wide$`Anguilla anguilla` <= shap_Aa_Cg_wide$`Cottus gobio`,
                                            3, 
                                            ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` < 0 & shap_Aa_Cg_wide$`Cottus gobio` < 0 & shap_Aa_Cg_wide$`Anguilla anguilla` >= shap_Aa_Cg_wide$`Cottus gobio`,
                                                   4, 5))))
                                            
ggplot(data = shap_Aa_Cg_wide, 
       aes(x = `Anguilla anguilla`, y = `Cottus gobio`, col = as.factor(col))) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) + 
  scale_colour_manual(values =  c('#2200c9', '#00b5c9', '#bd0f06', '#fa8eec', 'gray75')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', aspect.ratio = 0.5) + 
  xlab('Shapley values: Anguilla anguilla') + 
  ylab('Shapley values: Cottus gobio') + 
  geom_abline()


tm_shape(shap_rast_Aa['local_asym_cl_log10_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = F,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Anguilla anguilla')

tm_shape(shap_rast_Cg['local_asym_cl_log10_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = F,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Cottus gobio')


# read in presence predictions for both species
pres_Aa <- rast(sp_raster_pres_pa[grepl('Anguilla', sp_raster_pres_pa)])
pres_Cg <- rast(sp_raster_pres_pa[grepl('Cottus', sp_raster_pres_pa)])
pres_Aa_Cg <- mosaic(pres_Aa, pres_Cg)

shap_morepositive_Aa <- shap_rast_Aa['local_asym_cl_log10_SHAP'] > 0 & shap_rast_Cg['local_asym_cl_log10_SHAP'] > 0 & shap_rast_Aa['local_asym_cl_log10_SHAP'] >= shap_rast_Cg['local_asym_cl_log10_SHAP']
shap_morepositive_Cg <- shap_rast_Aa['local_asym_cl_log10_SHAP'] > 0 & shap_rast_Cg['local_asym_cl_log10_SHAP'] > 0 & shap_rast_Aa['local_asym_cl_log10_SHAP'] <= shap_rast_Cg['local_asym_cl_log10_SHAP']
shap_morenegative_Aa <- shap_rast_Aa['local_asym_cl_log10_SHAP'] < 0 & shap_rast_Cg['local_asym_cl_log10_SHAP'] < 0 & shap_rast_Aa['local_asym_cl_log10_SHAP'] <= shap_rast_Cg['local_asym_cl_log10_SHAP']
shap_morenegative_Cg <- shap_rast_Aa['local_asym_cl_log10_SHAP'] < 0 & shap_rast_Cg['local_asym_cl_log10_SHAP'] < 0 & shap_rast_Aa['local_asym_cl_log10_SHAP'] >= shap_rast_Cg['local_asym_cl_log10_SHAP']

shap_morepositive_Aa[shap_morepositive_Aa[]==F] <- NA
shap_morepositive_Aa[shap_morepositive_Aa[]==T] <- 1
shap_morepositive_Cg[shap_morepositive_Cg[]==F] <- NA
shap_morepositive_Cg[shap_morepositive_Cg[]==T] <- 2
shap_morenegative_Aa[shap_morenegative_Aa[]==F] <- NA
shap_morenegative_Aa[shap_morenegative_Aa[]==T] <- 3
shap_morenegative_Cg[shap_morenegative_Cg[]==F] <- NA
shap_morenegative_Cg[shap_morenegative_Cg[]==T] <- 4

shap_pos_neg_Aa_Cg <- mosaic(shap_morepositive_Aa, shap_morepositive_Cg, shap_morenegative_Aa, shap_morenegative_Cg)

levels(shap_pos_neg_Aa_Cg) <- data.frame(c(1,2,3,4), c('A. anguilla more positive', 'C. gobio more positive', 
                                                       'A. anguilla more negative', 'C. gobio more negative'))

tm_shape(shap_pos_neg_Aa_Cg) +
  tm_raster(palette = c('#2200c9', '#00b5c9', '#bd0f06', '#fa8eec'), 
            #labels = c('shap_morepositive_Aa', 'shap_morepositive_Cg', 'shap_morenegative_Aa', 'shap_morenegative_Cg'), 
            title = "Similar Connectivity effects") + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent")


# plot the environmental data
tm_shape(env_data['local_asym_cl_log10']) + 
  tm_raster(style = 'cont', 
            palette = '-viridis', 
            legend.is.portrait = F, 
            title = 'Connectivity') + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            legend.title.size = 1)



#### Example 3. Contrast response curves for low effect of cropland across target species ----

# read in rasters for focal variables
shap_rast_Lp <- rast(paste0(shap_dirs[grepl('Lampetra', shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_Ts <- rast(paste0(shap_dirs[grepl('Telestes', shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_Lp <- readRDS(shap_pa[grepl('Lampetra', shap_pa)])
shap_Ts <- readRDS(shap_pa[grepl('Telestes', shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_Lp <- left_join(shap_Lp %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_Ts <- left_join(shap_Ts %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_Lp_Ts <- bind_rows(shap_Lp, shap_Ts) %>% 
  select(species_name, TEILEZGNR, matches('stars_lu_crp_p_e_ele_residual'))

ggplot(data = shap_Lp_Ts, aes(x = stars_lu_crp_p_e_ele_residual, 
                              y = stars_lu_crp_p_e_ele_residual_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = T, method = 'gam') + 
  geom_hline(aes(yintercept = 0), lty = 2) + 
  scale_colour_manual(values = c('red', 'blue')) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5) +
  ylab('Shapley value') + 
  xlab('Environmental value') + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  ylim(c(-0.2, 0.2))

ggplot(data = shap_Bb_Cc, aes(x = ecoF_eco_mean, y = ecoF_eco_mean_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = T, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5) + 
  ylab('Shapley value') + 
  xlab('Environmental value') + 
  scale_colour_manual(values = c('gray50', 'gray75'))


S# make wide for plot of shapley values for each species against eachother
shap_Lp_Ts_wide <- shap_Lp_Ts %>% 
  pivot_wider(names_from = 'species_name', values_from = 'stars_lu_crp_p_e_ele_residual_SHAP')

ggplot(data = shap_Lp_Ts_wide, 
       aes(x = `Lampetra planeri`, y = `Telestes souffia`)) + 
  geom_jitter(col = 'gray75') + 
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) + 
  scale_colour_manual(values =  c('#2200c9', '#00b5c9', '#bd0f06', '#fa8eec', 'gray75')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', aspect.ratio = 0.5) + 
  xlab('Shapley values: Lampetra planeri') + 
  ylab('Shapley values: Telestes souffia') + 
  geom_abline() + 
  xlim(c(-0.2, 0.2)) + 
  ylim(c(-0.2, 0.2))



tm_shape(shap_rast_Lp['stars_lu_crp_p_e_ele_residual_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = F,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Lampetra planeri')

tm_shape(shap_rast_Ts['stars_lu_crp_p_e_ele_residual_SHAP']) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = F,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = T
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            title = 'Telestes souffia')


shap_Ts_Lp <- mean(c(shap_rast_Ts$stars_lu_crp_p_e_ele_residual_SHAP, shap_rast_Lp$stars_lu_crp_p_e_ele_residual_SHAP))

tm_shape(shap_Ts_Lp) +
  tm_raster(    style = "cont",
                palette = "RdBu",
                legend.reverse = F,
                title = 'Mean cropland effect',
                legend.is.portrait = F, 
                legend.show = T, 
                breaks = c(-0.2, -0.1, 0, 0.1, 0.2)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent")


# plot the environmental data
tm_shape(env_data['stars_lu_crp_p_e_ele_residual']) + 
  tm_raster(style = 'cont', 
            palette = '-viridis', 
            legend.is.portrait = F, 
            title = 'Cropland proportion cover', 
            breaks = c(quantile(env_data['stars_lu_crp_p_e_ele_residual'][], 0.025, na.rm = T), 
                       quantile(env_data['stars_lu_crp_p_e_ele_residual'][], 0.975, na.rm = T)), 
            midpoint = NA) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            legend.title.size = 1)

