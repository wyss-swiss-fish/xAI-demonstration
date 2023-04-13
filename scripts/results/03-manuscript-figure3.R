### Figure 3. Look at ways how disconnection occurs between species ----

### Pick interesting examples of pairs of species and variables to investigate ----

### Barbu or Abramis and Cyprinus carpio for ecomorphology (divergent responses)
### Cottus gobio and Anguilla anguilla for connectivity (similar responses but different magnitude)

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
RUN <- "ubelix_SDM_RF_MARCH_v1"

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

# load actual random forests
rf_pa <- paste0(sdm_dirs, "/output/final_models_rf_pa.RDS")

# read in environmental data
env_data <- rast(paste0(dd_env, '/ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif'))

#### Read in spatial data and assign simpler names to variables ----

# load spatial objects
source('scripts/results/00-load-spatial.R')

vars <- c('ecoF_discharge_max_log10', 
          'ecoF_slope_min_log10', 
          'ecoF_flow_velocity_mean', 
          'stars_t_mx_m_c', 
          'stars_t_mn_m_c',
          'local_asym_cl_log10',
          'local_dis2lake',
          'ecoF_eco_mean', 
          'stars_lu_crp_p_e_ele_residual',
          'local_tcd_0.1_log10_ele_residual', 
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
                 'connectivity', 'distance to lake', 'ecomorphology', 'cropland', 'tree cover', 'urbanisation',
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


#### Figure directories ----

contrast_dir <- 'figures-march2023-v2/constrast_species_examples/'
dir.create(contrast_dir, recursive = T)

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


# contrasting species response curves
Bb_Cc <- 'figures-march2023-v2/constrast_species_examples/Bb_Cc/'
dir.create(Bb_Cc, recursive =  T)

pdf(paste0(Bb_Cc, 'response_curves.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Bb_Cc, aes(x = ecoF_eco_mean, y = ecoF_eco_mean_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = F, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4)) + 
  geom_hline(yintercept = 0, lwd = 0.25) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.2,0.8), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  ylab('Shapley value') + 
  xlab('Environmental value') + 
  scale_colour_manual(values = c('gray50', 'gray75'))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_Bb_Cc_wide <- shap_Bb_Cc %>% 
  pivot_wider(names_from = 'species_name', values_from = 'ecoF_eco_mean_SHAP')

shap_Bb_Cc_wide$col <- ifelse(shap_Bb_Cc_wide$`Barbus barbus` >= shap_Bb_Cc_wide$`Cyprinus carpio` & 
                                shap_Bb_Cc_wide$`Barbus barbus` > 0 & shap_Bb_Cc_wide$`Cyprinus carpio` < 0,
                              1, 
                              ifelse(shap_Bb_Cc_wide$`Barbus barbus` <= shap_Bb_Cc_wide$`Cyprinus carpio` & 
                                       shap_Bb_Cc_wide$`Cyprinus carpio` > 0 & shap_Bb_Cc_wide$`Barbus barbus` < 0, 
                                     2, 3))

pdf(paste0(Bb_Cc, 'shapley_biplot.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Bb_Cc_wide, 
       aes(x = `Barbus barbus`, y = `Cyprinus carpio`, col = as.factor(col))) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  geom_vline(aes(xintercept = 0), lwd = 0.25) + 
  geom_abline(lwd = 0.25) + 
  scale_colour_manual(values = c('#bd0f06', '#2200c9', 'gray90')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  xlab('Shapley values: Barbus barbus') + 
  ylab('Shapley values: Cyprinus carpio')
dev.off()

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

pdf(paste0(Bb_Cc, 'contrast_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_div_Bb_Cc) +
  tm_raster(n = 2,
            palette = c('#bd0f06', '#2200c9'), 
            labels = c('Barbus barbus', 'Cyprinus carpio'), 
            title = "Ecomorphology positive for:") + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F, 
            legend.text.size = 1, 
            legend.title.size = 1.5)
dev.off()

# plot the environmental data
pdf(paste0(Bb_Cc, 'raw_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(env_data['ecoF_eco_mean']) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = 'Ecomorphological modification', 
            legend.show = F) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

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

# contrasting species response curves
Aa_Cg <- 'figures-march2023-v2/constrast_species_examples/Aa_Cg/'
dir.create(Aa_Cg, recursive =  T)

pdf(paste0(Aa_Cg, 'response_curves.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Aa_Cg, aes(x = local_asym_cl_log10, y = local_asym_cl_log10_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black',se = F, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4)) + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.2,0.8), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  ylab('Shapley value') + 
  xlab('Environmental value')
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_Aa_Cg_wide <- shap_Aa_Cg %>% 
  pivot_wider(names_from = 'species_name', values_from = 'local_asym_cl_log10_SHAP')

shap_Aa_Cg_wide$col <- ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` > 0 & shap_Aa_Cg_wide$`Cottus gobio` > 0, 1, 
                              ifelse(shap_Aa_Cg_wide$`Anguilla anguilla` < 0 & shap_Aa_Cg_wide$`Cottus gobio` < 0, 2, 3))

pdf(paste0(Aa_Cg, 'shapley_biplot.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Aa_Cg_wide, 
       aes(x = `Anguilla anguilla`, y = `Cottus gobio`, col = (`Anguilla anguilla` + `Cottus gobio`) / 2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  geom_vline(aes(xintercept = 0), lwd = 0.25) + 
  scale_colour_gradientn(
    colors=c('#bd0f06','gray90', '#2200c9'),
    values=scales::rescale(c(-0.1, -0.02, 0, 0.02 , 0.1)),
    limits=c(-max(abs(range(shap_Aa_Cg$local_asym_cl_log10_SHAP, na.rm = T))), 
             max(abs(range(shap_Aa_Cg$local_asym_cl_log10_SHAP, na.rm = T))))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        aspect.ratio = 0.5,  
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  xlab('Shapley values: Anguilla anguilla') + 
  ylab('Shapley values: Cottus gobio') + 
  geom_abline(lwd = 0.25)
dev.off()

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

shap_mean_Aa_Cg <- c(shap_rast_Aa['local_asym_cl_log10_SHAP'], shap_rast_Cg['local_asym_cl_log10_SHAP'] ) %>% mean

pdf(paste0(Aa_Cg, 'contrast_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_mean_Aa_Cg) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = 'Mean connectivity shapley value',
            legend.is.portrait = F, 
            legend.show = T, 
            breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
            midpoint = 0
  ) + 
  
  #tm_raster(palette = c('#2200c9', '#00b5c9', '#bd0f06', '#fa8eec'), 
  #          style = 'cont',) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# plot the environmental data
pdf(paste0(Aa_Cg, 'raw_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(env_data['local_asym_cl_log10']) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = 'Connectivity', 
            legend.show = F) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()


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

# contrasting species response curves
Lp_Ts <- 'figures-march2023-v2/constrast_species_examples/Lp_Ts/'
dir.create(Lp_Ts, recursive =  T)

pdf(paste0(Lp_Ts, 'response_curves.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Lp_Ts, aes(x = stars_lu_crp_p_e_ele_residual, 
                              y = stars_lu_crp_p_e_ele_residual_SHAP)) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = F, method = 'gam') + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.2,0.8), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  ylab('Shapley value') + 
  xlab('Environmental value') + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  ylim(c(-0.1, 0.1))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_Lp_Ts_wide <- shap_Lp_Ts %>% 
  pivot_wider(names_from = 'species_name', values_from = 'stars_lu_crp_p_e_ele_residual_SHAP')

pdf(paste0(Lp_Ts, 'shapley_biplot.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_Lp_Ts_wide, 
       aes(x = `Lampetra planeri`, y = `Telestes souffia`, col = (`Lampetra planeri` + `Telestes souffia`)/2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  geom_vline(aes(xintercept = 0), lwd = 0.25) + 
  scale_colour_gradientn(
    colors=c('#bd0f06','gray90', '#2200c9'),
    values=scales::rescale(c(-0.2, -0.1, 0, 0.1, 0.2)),
    limits=c(-max(abs(range(shap_Lp_Ts$stars_lu_crp_p_e_ele_residual_SHAP, na.rm = T))), 
             max(abs(range(shap_Lp_Ts$stars_lu_crp_p_e_ele_residual_SHAP, na.rm = T))))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        aspect.ratio = 0.5, 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  xlab('Shapley values: Lampetra planeri') + 
  ylab('Shapley values: Telestes souffia') + 
  geom_abline(lwd = 0.25) + 
  xlim(c(-0.1, 0.1)) + 
  ylim(c(-0.1, 0.1))
dev.off()

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

pdf(paste0(Lp_Ts, 'contrast_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_Ts_Lp) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = 'Mean cropland shapley value',
            legend.is.portrait = F, 
            legend.show = T, 
            breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
            midpoint = 0
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(Lp_Ts, 'raw_map.pdf'), width = 5, height = 5, bg = 'transparent')
# plot the environmental data
tm_shape(env_data['stars_lu_crp_p_e_ele_residual']) + 
  tm_raster(style = 'cont', 
            palette = c('#bd0f06','gray90', '#2200c9'), 
            legend.is.portrait = F, 
            title = 'Cropland proportion cover', 
            breaks = c(quantile(env_data['stars_lu_crp_p_e_ele_residual'][], 0.025, na.rm = T), 
                       quantile(env_data['stars_lu_crp_p_e_ele_residual'][], 0.975, na.rm = T)), 
            midpoint = NA, 
            legend.show = F) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            legend.title.size = 1, 
            frame = F)
dev.off()


#### Catchment force plots ----


# get the teilenzugsgebeit for the sense
sense <- data.frame(catchment = 'Sense', TEZGNR40 = c(100184, 101484, 105044, 105663, 102204, 102377))
emme <- data.frame(catchment = 'Emme', TEZGNR40 = c(106570, 106299, 101728, 101280, 106548, 103982, 108894))
catchments <- rbind(sense, emme)

# read in subcatchments, transform, union
subcatchments_sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% catchments$TEZGNR40) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

mapview::mapview(subcatchments_sense_union)

subcatchment_2km_names <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% catchments$TEZGNR40) %>%
  st_drop_geometry() %>% 
  select(TEILEZGNR, TEZGNR40)
subcatchment_2km_names

# join with focal catchments
subcatchment_2km_names <- left_join(subcatchment_2km_names, catchments)

# read in subcatchment data and filter
shapley_focal_list <- lapply(1:length(sp_list), function(x){ 
  
  if(file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    
    return(shap_x)
    
  }
}
)

# bind outputs for species
shapley_focal <- bind_rows(shapley_focal_list) %>% tibble

# join in subcatchment information
shapley_focal <- left_join(shapley_focal, subcatchment_2km_names)

# read in baseline predictions from raw models 
baseline_values <- sapply(1:length(sp_list), function(x){ 
  
  if(file.exists(rf_pa[x])){
    
    # read in raw shapleys
    baseline_model <- readRDS(rf_pa[x])
    baseline_value <- mean(as.numeric(predict(baseline_model, type = 'prob')[,2]), na.rm = T)
    
    return(baseline_value)
    
  }
}
)

# make baseline prediction dataframe
baseline_prediction <- data.frame(species_name = sp_list, baseline_prediction = baseline_values)

# get mean prediction per focal catchment
catchment_prediction <- shapley_focal %>% 
  filter(!is.na(catchment)) %>% 
  group_by(species_name, catchment) %>% 
  do(catchment_prediction = mean(.$suitability, na.rm = T)) %>% 
  unnest(c(catchment_prediction))

# quick example plot to see if the data are in the correct format
ggplot(data = left_join(catchment_prediction, baseline_prediction)) + 
  geom_point(aes(y = species_name, x = baseline_prediction), col = 'black') + 
  geom_point(aes(y = species_name, x = catchment_prediction), col = 'blue') + 
  facet_wrap(~catchment)

# join catchment shapley summaries to the suitability data
shapley_catchment_summary <- shapley_focal %>% 
  tibble %>% 
  filter(!is.na(catchment)) %>% 
  select(species_name, catchment, matches('_SHAP')) %>% 
  pivot_longer(col = c(-species_name, -catchment)) %>% 
  group_by(species_name, catchment , name) %>% 
  do(mean_shap = mean(.$value, na.rm = T), 
     sign_shap = sign(mean(.$value, na.rm = T))) %>% 
  unnest(c(mean_shap, sign_shap)) %>% 
  left_join(., baseline_prediction) %>% 
  left_join(., catchment_prediction)

# quick example plot to see if the data are in the correct format
ggplot(data = shapley_catchment_summary) + 
  geom_point(aes(y = species_name, x = (baseline_prediction)), col = 'black') + 
  geom_point(aes(y = species_name, x = (catchment_prediction)), col = 'blue') + 
  facet_wrap(~catchment)

# join in renaming object
shapley_catchment_summary <- left_join(shapley_catchment_summary, data.frame(vars_renamed), by = c('name' = 'vars_shap'))

# order by mean value
levels <- shapley_catchment_summary %>% 
  group_by(vars_renamed) %>% 
  summarise(mean_shap2 = mean(mean_shap, na.rm = T)) %>% 
  arrange(mean_shap2) %>% 
  pull(vars_renamed)

shapley_catchment_summary$vars_renamed <- factor(shapley_catchment_summary$vars_renamed, 
                                                 levels = levels)

ggplot(data = shapley_catchment_summary) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        axis.line.x = element_line(), 
        panel.border = element_rect(fill = 'transparent'), 
        axis.text.x = element_text(size = 7), 
        axis.ticks = element_line(), 
        strip.text.x = element_text(face = 'italic'), 
        strip.text.y = element_text(face = 'bold', size = 12), 
        aspect.ratio = 1.75) + 
  facet_grid(catchment~species_name, scales = 'free_x') + 
  geom_segment(aes(xend = mean_shap + baseline_prediction, 
                   x = baseline_prediction, 
                   y = vars_renamed, yend = vars_renamed, 
                   col = mean_shap, lwd = abs(mean_shap)),
               lineend = 'butt', linejoin = 'bevel') + 
  geom_point(aes(x = mean_shap + baseline_prediction, 
                 y = vars_renamed, 
                 col = mean_shap, 
                 size = abs(mean_shap)*2)) + 
  geom_vline(aes(xintercept = baseline_prediction), size = 1, col = 'gray75') + 
  scale_colour_gradientn(
                           colors=c('#bd0f06','gray90', '#2200c9'),
                           values=scales::rescale(c(-0.2, -0.02, 0, 0.02 , 0.2)),
                           limits=c(-max(abs(range(shapley_catchment_summary$mean_shap, na.rm = T))), 
                                    max(abs(range(shapley_catchment_summary$mean_shap, na.rm = T))))) +
  scale_size_continuous(range = c(0.1, 3)) + 
  xlab('Change in local prediction from average prediction') + 
  ylab(NULL)



  
#### Make catchment specific response curves ----

shap_focal_catch_rc <- shapley_focal %>%  
  left_join(all_env_subcatchments) %>% 
  select(TEILEZGNR, species_name, matches('local_flood'), matches('local_asym_cl_log10')) %>% 
  left_join(subcatchment_2km_names)

ggplot(data = shap_focal_catch_rc) + 
  geom_hline(yintercept = 0) + 
  geom_point(aes(x = local_flood, y = local_flood_SHAP), 
             alpha = 0.1, stroke = 0, pch = 19, col = 'gray50') + 
  geom_point(data = shap_focal_catch_rc %>% filter(!is.na(catchment)) %>% sample_n(., size = nrow(.)), 
             aes(x = local_flood, y = local_flood_SHAP, col = catchment),
             pch = 19, stroke = 0, size = 2) + 
  geom_boxplot(data = shap_focal_catch_rc %>% filter(!is.na(catchment)) %>% sample_n(., size = nrow(.)), 
               aes(x = local_flood, y = -0.1, col = catchment, group = catchment), 
               width = 0.05) + 
  stat_smooth(data = shap_focal_catch_rc, 
              aes(x = local_flood, y = local_flood_SHAP), col = 'gray5') + 
  facet_wrap(~species_name, nrow = 1) + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(face = 'italic'), 
        axis.text = element_text(size = 7)) + 
  xlab('Floodplain proportion') + 
  ylab('Floodplain proportion \n shapley value') + 
  scale_colour_manual(values = c('#9FC131FF', '#FDD20EFF')) +
  scale_x_continuous(breaks = c(0, 0.5, 1))


ggplot(data = shap_focal_catch_rc) + 
  geom_hline(yintercept = 0) + 
  geom_point(aes(x = local_asym_cl_log10, y = local_asym_cl_log10_SHAP), 
             alpha = 0.1, stroke = 0, pch = 19, col = 'gray50') + 
  geom_point(data = shap_focal_catch_rc %>% filter(!is.na(catchment)) %>% sample_n(., size = nrow(.)), 
             aes(x = local_asym_cl_log10, y = local_asym_cl_log10_SHAP, col = catchment),
             pch = 19, stroke = 0, size = 2) + 
  geom_boxplot(data = shap_focal_catch_rc %>% filter(!is.na(catchment)) %>% sample_n(., size = nrow(.)), 
             aes(x = local_asym_cl_log10, y = 0.1, col = catchment, group = catchment), 
             width = 0.1) + 
  stat_smooth(data = shap_focal_catch_rc, 
              aes(x = local_asym_cl_log10, y = local_asym_cl_log10_SHAP), col = 'gray5') + 
  facet_wrap(~species_name, nrow = 1) + 
  theme_bw() +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(face = 'italic'), 
        axis.text = element_text(size = 7)) + 
  xlab('Connectivity') + 
  ylab('Connectivity \n shapley value') + 
  scale_colour_manual(values = c('#9FC131FF', '#FDD20EFF')) +
  scale_x_continuous(breaks = c(-3, -2, -1))



