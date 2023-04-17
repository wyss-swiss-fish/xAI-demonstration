### Figure 3. Look at ways how disconnection occurs between species ----

# Here we explore two categories of responses that are interesting to explore from shapley values:
# - coupling or decoupling of species responses to threats
# - convergence or divergence of threat responses between species 

# This can lead to four scenrious
# coupled and convergent responses between species
# coupled and divergent responses between species
# decoupling responses for both species
# decoupling responses for one species and divergence. 

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

#### 1. Set directories for loading and saving objects ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_MARCH_v6/"
contrast_dir <- paste0(fig_dir, 'constrast_species_examples/')
dir.create(contrast_dir, recursive = T)

# get run to mak figures for
RUN <- "ubelix_SDM_RF_MARCH_v5"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS('figures/ubelix_SDM_RF_MARCH_v6/evaluations/subset_sp.RDS')

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN), full.names = T)
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

# load spatial objects
source('scripts/results/00-load-spatial.R')

#### 2. Set variable names ----

# set vector of focal variables for assigning new names, this must be consistent with the models fitted
vars <- c('ecoF_discharge_max_log10', 
          'ecoF_slope_min_log10', 
          'ecoF_flow_velocity_mean', 
          'stars_t_mx_m_c', 
          'stars_t_mn_m_c',
          'local_asym_cl_log10',
          'local_dis2lake',
          'ecoF_eco_mean', 
          'local_imd_log10',
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

#### 3. Read in all shapley response curves and plot for investigation across our focal species together ----

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

all_shap_env <- bind_rows(all_shap_env)


# make plot of all species and all response curves
pdf(file = paste0(contrast_dir, 'all_response_curves_linear.pdf'), width = 15, height = 15)
ggplot(data = all_shap_env %>% 
         group_by(species_name, vars_renamed)) + 
  stat_smooth(aes(x = env, y = shapley, col = species_name), method = 'lm') +
  facet_wrap(~vars_renamed, scales = 'free') + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()

pdf(file = paste0(contrast_dir, 'all_response_curves_nonlinear.pdf'), width = 15, height = 15)
ggplot(data = all_shap_env %>% 
         group_by(species_name, vars_renamed)) + 
  stat_smooth(aes(x = env, y = shapley, col = species_name)) +
  facet_wrap(~vars_renamed, scales = 'free') + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()


#### Example 1. Contrast response curves: coupled but divergent -----

# contrast species
CouDiv1 <- 'Thymallus thymallus'
CouDiv2 <- 'Perca fluviatilis'

# variable name
CouDiv_var <- 'local_dis2lake'
# shapley name
CouDiv_shap <- 'local_dis2lake_SHAP'

# environmental label
env_label = 'distance to lake'
shap_label = expression(phi1~'distance to lake')

# plot titles
title_CouDiv_mean = expression('Mean distance to lake' ~ phi1)
title_CouDiv_sd   = expression('S.D. distance to lake' ~ phi1)
title_CouDiv_raw  = 'Raw distance to lake'


# read in rasters for focal variables
shap_rast_CouDiv1 <- rast(paste0(shap_dirs[grepl(CouDiv1, shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_CouDiv2 <- rast(paste0(shap_dirs[grepl(CouDiv2, shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_CouDiv1 <- readRDS(shap_pa[grepl(CouDiv1, shap_pa)])
shap_CouDiv2 <- readRDS(shap_pa[grepl(CouDiv2, shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_CouDiv1 <- left_join(shap_CouDiv1 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_CouDiv2 <- left_join(shap_CouDiv2 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_CouDiv12 <- bind_rows(shap_CouDiv1, shap_CouDiv2) %>% 
  select(species_name, TEILEZGNR, matches(CouDiv_var))


# contrasting species response curves
CouDiv <- paste0(contrast_dir, '/Coupled_Divergent/')
dir.create(CouDiv, recursive =  T)

pdf(paste0(CouDiv, 'response_curves.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_CouDiv12%>% sample_n(., nrow(.)), 
       aes(x = .data[[CouDiv_var]], y = .data[[CouDiv_shap]])) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black', se = F, method = 'gam') + 
  geom_hline(yintercept = 0, lwd = 0.25) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.5,0.8), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  ylab(shap_label) + 
  xlab(env_label) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_CouDiv12_wide <- shap_CouDiv12 %>% 
  pivot_wider(names_from = 'species_name', values_from = CouDiv_shap)

shap_CouDiv12_wide$col <- ifelse(shap_CouDiv12_wide[CouDiv1] >= shap_CouDiv12_wide[CouDiv2] & 
                                shap_CouDiv12_wide[CouDiv1] > 0 & shap_CouDiv12_wide[CouDiv2] < 0,
                              1, 
                              ifelse(shap_CouDiv12_wide[CouDiv1] <= shap_CouDiv12_wide[CouDiv2] & 
                                       shap_CouDiv12_wide[CouDiv2] > 0 & shap_CouDiv12_wide[CouDiv1] < 0, 
                                     2, 3))

# set the limits
min_axis <- min(c(shap_CouDiv12_wide[CouDiv2], shap_CouDiv12_wide[CouDiv1]), na.rm = T)
max_axis <- max(c(shap_CouDiv12_wide[CouDiv2], shap_CouDiv12_wide[CouDiv1]), na.rm = T)

pdf(paste0(CouDiv, 'shapley_biplot.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_CouDiv12_wide, 
       aes(x = .data[[CouDiv1]], y = .data[[CouDiv2]], col = as.factor(col))) + 
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
  xlab(bquote(phi1~ .(CouDiv1))) + 
  ylab(bquote(phi1~ .(CouDiv2))) + 
  ylim(c(min_axis, max_axis)) + 
  xlim(c(min_axis, max_axis))
dev.off()

pdf(paste0(CouDiv, '/CouDiv1_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_rast_CouDiv1[CouDiv_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouDiv1,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 5), 1)
            ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(CouDiv, '/CouDiv2_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_rast_CouDiv2[CouDiv_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouDiv2,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 5), 1)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()


# read in presence predictions for both species
pres_CouDiv1 <- rast(sp_raster_pres_pa[grepl(CouDiv1, sp_raster_pres_pa)])
pres_CouDiv2 <- rast(sp_raster_pres_pa[grepl(CouDiv2, sp_raster_pres_pa)])
pres_CouDiv12 <- mosaic(pres_CouDiv1, pres_CouDiv12)

shap_div_CouDiv1 <- shap_rast_CouDiv1[CouDiv_shap] > shap_rast_CouDiv2[CouDiv_shap] & shap_rast_CouDiv1[CouDiv_shap] > 0
shap_div_CouDiv2 <- shap_rast_CouDiv2[CouDiv_shap] > shap_rast_CouDiv1[CouDiv_shap] & shap_rast_CouDiv2[CouDiv_shap] > 0

shap_div_CouDiv1[shap_div_CouDiv1[]==F] <- NA
shap_div_CouDiv1[shap_div_CouDiv1[]==T] <- 1
shap_div_CouDiv2[shap_div_CouDiv2[]==F] <- NA
shap_div_CouDiv2[shap_div_CouDiv2[]==T] <- 2

shap_div_CouDiv12 <- mosaic(shap_div_CouDiv1, shap_div_CouDiv2)
pres_CouDiv12 <- resample(pres_CouDiv12, shap_div_CouDiv12)
shap_div_CouDiv12[is.na(pres_CouDiv12)] <- NA

pdf(paste0(CouDiv, '/contrast_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_div_CouDiv12) +
  tm_raster(n = 2,
            palette = c('#bd0f06', '#2200c9'), 
            labels = c(CouDiv1, CouDiv2), 
            title =expression("positive" ~ phi1 ~ "for distance to lake")) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F, 
            legend.text.size = 1, 
            legend.title.size = 1.5,
            legend.width = 2)
dev.off()


## take average of species shapley values and plot
shap_mean_CouDiv12 <- app(c(shap_rast_CouDiv1[CouDiv_shap], shap_rast_CouDiv2[CouDiv_shap] ), mean)

pdf(paste0(CouDiv, 'mean_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_mean_CouDiv12) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = title_CouDiv_mean,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0
            #breaks = c(-0.2, -0.1, 0, 0.1, 0.2)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

## take the sd of species shapley values
## take average of species shapley values and plot
shap_sd_CouDiv12 <- app(c(shap_rast_CouDiv1[CouDiv_shap], shap_rast_CouDiv2[CouDiv_shap] ), sd)

pdf(paste0(CouDiv, 'sd_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_sd_CouDiv12) +
  tm_raster(style = "cont",
            palette = c('gray90', '#bd0f06'),
            legend.reverse = F,
            title = title_CouDiv_sd,
            legend.is.portrait = F, 
            legend.show = T, 
            breaks = c(0, 0.1, 0.2, 0.25)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# plot the environmental data
pdf(paste0(CouDiv, 'raw_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(env_data[CouDiv_var]) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = title_CouDiv_raw, 
            legend.show = T) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

#### Example 2. Contrast response curves for connnectivity across species ----

# contrast species
CouCon1 <- 'Lampetra planeri'
CouCon2 <- 'Barbus barbus'

# variable name
CouCon_var <- 'local_asym_cl_log10'
# shapley name
CouCon_shap <- 'local_asym_cl_log10_SHAP'

# environmental label
env_label = 'connectivity'
shap_label = expression(phi1~'conectivity')

# plot titles
title_CouCon_mean = expression('Mean connectivity' ~ phi1)
title_CouCon_sd   = expression('S.D. connectivity' ~ phi1)
title_CouCon_raw  = 'Raw connectivity'


# read in rasters for focal variables
shap_rast_CouCon1 <- rast(paste0(shap_dirs[grepl(CouCon1, shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_CouCon2 <- rast(paste0(shap_dirs[grepl(CouCon2, shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_CouCon1 <- readRDS(shap_pa[grepl(CouCon1, shap_pa)])
shap_CouCon2 <- readRDS(shap_pa[grepl(CouCon2, shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_CouCon1 <- left_join(shap_CouCon1 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_CouCon2 <- left_join(shap_CouCon2 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_CouCon12 <- bind_rows(shap_CouCon1, shap_CouCon2) %>% 
  select(species_name, TEILEZGNR, matches(CouCon_var))

# contrasting species response curves
CouCon <- paste0(contrast_dir, 'Coupled_Convergent/')
dir.create(CouCon, recursive =  T)

pdf(paste0(CouCon, 'response_curves.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_CouCon12 %>% sample_n(., nrow(.)), 
       aes(x = .data[[CouCon_var]], y = .data[[CouCon_shap]])) + 
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
  ylab(shap_label) + 
  xlab(env_label) + 
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_CouCon12_wide <- shap_CouCon12 %>% 
  pivot_wider(names_from = 'species_name', values_from = all_of(CouCon_shap))

shap_CouCon12_wide$col <- ifelse(shap_CouCon12_wide[[CouCon1]] > 0 & shap_CouCon12_wide[[CouCon2]] > 0, 1, 
                              ifelse(shap_CouCon12_wide[[CouCon1]] < 0 & shap_CouCon12_wide[[CouCon2]] < 0, 2, 3))


# set the limits
max_axis <- max(c(shap_CouCon12_wide[CouCon2], shap_CouCon12_wide[CouCon1]), na.rm = T)
min_axis <- min(c(shap_CouCon12_wide[CouCon2], shap_CouCon12_wide[CouCon1]), na.rm = T)

pdf(paste0(CouCon, 'shapley_biplot.pdf'), width = 5, height = 5, bg = 'transparent')
ggplot(data = shap_CouCon12_wide, 
       aes(x = .data[[CouCon1]], y = .data[[CouCon2]], col = (.data[[CouCon1]] + .data[[CouCon2]]) / 2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  geom_vline(aes(xintercept = 0), lwd = 0.25) + 
  scale_colour_gradientn(
    colors=c('#bd0f06','gray90', '#2200c9'),
    values=scales::rescale(c(-0.1, -0.02, 0, 0.02 , 0.1)),
    limits=c(-max(abs(range(shap_CouCon12_wide[c(CouCon1, CouCon2)], na.rm = T))), 
             max(abs(range(shap_CouCon12_wide[c(CouCon1, CouCon2)], na.rm = T))))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        aspect.ratio = 0.5,  
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank()) + 
  xlab(bquote(phi1~ .(CouCon1))) + 
  ylab(bquote(phi1~ .(CouCon2))) + 
  geom_abline(lwd = 0.25) + 
  xlim(c(min_axis, max_axis)) + 
  ylim(c(min_axis, max_axis))
dev.off()

pdf(paste0(CouCon, '/CouCon1_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_rast_CouCon1[CouCon_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouCon1,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 5), 1)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(CouCon, '/CouCon2_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_rast_CouCon2[CouCon_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouCon2,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 5), 1)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# read in presence predictions for both species
pres_CouCon1 <- rast(sp_raster_pres_pa[grepl(CouCon1, sp_raster_pres_pa)])
pres_CouCon2 <- rast(sp_raster_pres_pa[grepl(CouCon2, sp_raster_pres_pa)])
pres_CouCon12 <- mosaic(pres_CouCon1, pres_CouCon2)

# take the mean of the shapley values
shap_mean_CouCon12 <- app(c(shap_rast_CouCon1[CouCon_shap], shap_rast_CouCon2[CouCon_shap]), mean)

# make maps of mean shapley values across species
pdf(paste0(CouCon, 'mean_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_mean_CouCon12) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = title_CouCon_mean,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# take the sd of the shapley values
shap_sd_CouCon12 <- app(c(shap_rast_CouCon1[CouCon_shap], shap_rast_CouCon2[CouCon_shap]), sd)

# make maps of mean shapley values across species
pdf(paste0(CouCon, 'sd_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(shap_sd_CouCon12) +
  tm_raster(style = "cont",
            palette = c('gray90', '#bd0f06'),
            legend.reverse = F,
            title = title_CouCon_sd,
            legend.is.portrait = F, 
            legend.show = T, 
            breaks = c(0, 0.1, 0.2, 0.25)
  ) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_borders(col = "gray75", lwd = 0.01) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()


# plot the environmental data
pdf(paste0(CouCon, 'raw_map.pdf'), width = 5, height = 5, bg = 'transparent')
tm_shape(env_data['local_asym_cl_log10']) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = 'Connectivity', 
            legend.show = T) + 
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



