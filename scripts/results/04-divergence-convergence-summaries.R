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
dd     <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch  <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_APRIL_V1_02/"
contrast_dir <- paste0(fig_dir, 'constrast_species_examples/')
dir.create(contrast_dir, recursive = T)

# get run to mak figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN), full.names = T)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]

# create shapley folders per species
shap_po <- paste0(shap_dirs, "/shapley_rf_po.RDS")
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")

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

#### 2. Set variable names ----

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
                 'morph. mod.', 
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


## PLOTS OF OVERALL RELATIONSHIPS BETWEEN SHAPLEY VALUES AND ENVIRONMENTAL VALUES AS A FIRST LOOK
# linear response curves
pdf(file = paste0(contrast_dir, 'all_response_curves_linear.pdf'), width = 15, height = 15)
ggplot(data = all_shap_env %>% 
         group_by(species_name, vars_renamed)) + 
  stat_smooth(aes(x = env, y = shapley, col = species_name), method = 'lm') +
  facet_wrap(~vars_renamed, scales = 'free') + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()

# non-linear response curves
pdf(file = paste0(contrast_dir, 'all_response_curves_nonlinear.pdf'), width = 15, height = 15)
ggplot(data = all_shap_env %>% 
         group_by(species_name, vars_renamed)) + 
  stat_smooth(aes(x = env, y = shapley, col = species_name)) +
  facet_wrap(~vars_renamed, scales = 'free') + 
  geom_hline(aes(yintercept = 0)) + 
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()



#### 4. Overall convergence, divergence, coupling and decoupling ----

# calculate correlation between variables
all_cors <- all_shap_env %>% 
  group_by(species_name, vars_renamed) %>% 
  do(cor = as.numeric(cor(x = .$shapley, y = .$env, method = 'pearson'))) %>% 
  unnest()
  
# make wider and fill with 0s
all_cors <- pivot_longer(pivot_wider(all_cors, 
                                     names_from = 'vars_renamed', values_from = cor, values_fill = 0), 
                         cols = -species_name, 
                         names_to = 'vars_renamed', values_to = 'cor')

# get levels
relevel_shap <- all_shap_env %>% 
  group_by(vars_renamed) %>% 
  do(mean_shap = mean(abs(.$shapley))) %>% 
  unnest() %>% 
  arrange(mean_shap)

# relevel factors
all_cors$vars_renamed <- factor(all_cors$vars_renamed, levels = relevel_shap$vars_renamed)

#install.packages("seecolor")
library(seecolor)
pal <- c('#B436FF', '#E83197', '#FF6242', '#E88C31', '#FFCC35', 
         '#E2EB9E', '#65EBA5', '#4DE3FF', '#6772FF' )
print_color(pal)

pdf(paste0(contrast_dir, 'correlation_env_shap.pdf'), height = 6, width = 8)
ggplot(data = all_cors) + 
  geom_violin(data = all_cors %>% filter(cor != 0), aes(x = cor, y = vars_renamed), 
              scale = 'width', fill = 'gray75') + 
  geom_point(aes(x = cor, y = vars_renamed, fill = species_name),
             pch = 21, size=4, position = position_jitter(width = 0, height = 0.2)) + 
  geom_vline(aes(xintercept = 0)) +
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 13), 
        axis.text.x = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13)) +
  xlab(expression('Correlation between environment and'~phi)) +
  ylab(NULL) + 
  scale_fill_manual(' ', values = pal) + 
  guides(fill = guide_legend(override.aes = list(size=5)))
dev.off()


  

#### 5. Example 1. Contrast response curves: coupled but divergent -----

# contrast species
CouDiv1 <- 'Oncorhynchus mykiss'
CouDiv2 <- 'Squalius cephalus'

CouDiv1_lab <- 'O. mykiss'
CouDiv2_lab <- 'S. cephalus'

# variable name
CouDiv_var <- 'ecoF_eco_mean_ele_residual'
# shapley name
CouDiv_shap <- 'ecoF_eco_mean_ele_residual_SHAP'

# environmental label
env_label = 'morph.'
shap_label = expression(phi1~'morph.')

# plot titles
title_CouDiv_mean = expression('Mean morphological mod.' ~ phi1)
title_CouDiv_sd   = expression('S.D. morphological mod.' ~ phi1)
title_CouDiv_raw  = 'morph.'


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

# shorten names
shap_CouDiv12 <- shap_CouDiv12 %>% mutate(species_name_lab = case_when(.$species_name == CouDiv1 ~ CouDiv1_lab, 
                                                                       .$species_name == CouDiv2 ~ CouDiv2_lab))

# contrasting species response curves
CouDiv <- paste0(contrast_dir, '/Coupled_Divergent/')
dir.create(CouDiv, recursive =  T)

pdf(paste0(CouDiv, 'response_curves.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_CouDiv12%>% sample_n(., nrow(.)), 
       aes(x = .data[[CouDiv_var]], y = .data[[CouDiv_shap]])) + 
  geom_jitter(aes(pch = species_name_lab, col = species_name_lab)) + 
  stat_smooth(aes(group = species_name_lab, lty = species_name_lab), 
              col = 'black', se = F, method = 'gam', size = 0.5) + 
  geom_hline(yintercept = 0, lwd = 0.25) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.6,0.9), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank(), 
        text = element_text(size = 14)) + 
  ylab(shap_label) + 
  xlab(env_label) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_CouDiv12_wide <- shap_CouDiv12 %>%
  select(-species_name_lab) %>% 
  pivot_wider(names_from = 'species_name', values_from = CouDiv_shap)

# set the limits
max_axis <- max(c(shap_CouDiv12_wide[[CouDiv2]], shap_CouDiv12_wide[[CouDiv1]]), na.rm = T)
min_axis <- min(c(shap_CouDiv12_wide[[CouDiv2]], shap_CouDiv12_wide[[CouDiv1]]), na.rm = T)


# make biplot of shapley
pdf(paste0(CouDiv, 'shapley_biplot.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_CouDiv12_wide, 
       aes(x = .data[[CouDiv1]], y = .data[[CouDiv2]], col = (.data[[CouDiv1]] + .data[[CouDiv2]]) / 2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25, col = 'gray50') + 
  geom_vline(aes(xintercept = 0), lwd = 0.25, col = 'gray50') + 
  geom_abline(lwd = 0.25, col = 'gray50') + 
  #scale_colour_manual(values = c('#bd0f06', '#2200c9', 'gray90')) + 
  scale_colour_gradientn(
    colors=c('#bd0f06','gray90', '#2200c9'),
    values=scales::rescale(c(-0.1, -0.02, 0, 0.02 , 0.1)),
    limits=c(-max(abs(range(shap_CouDiv12_wide[c(CouDiv1, CouDiv2)], na.rm = T))), 
             max(abs(range(shap_CouDiv12_wide[c(CouDiv1, CouDiv2)], na.rm = T))))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank(), 
        text = element_text(size = 14)) + 
  xlab(bquote(phi1~ .(paste0(' ', CouDiv1_lab)))) + 
  ylab(bquote(phi1~ .(paste0(' ', CouDiv2_lab)))) + 
  ylim(c(min_axis, max_axis)) + 
  xlim(c(min_axis, max_axis))
dev.off()

pdf(paste0(CouDiv, '/CouDiv1_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_CouDiv1[CouDiv_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouDiv1,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(min(shap_rast_CouDiv1[CouDiv_shap][], na.rm = T), max(shap_rast_CouDiv1[CouDiv_shap][], na.rm = T), length.out = 4),1)) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(CouDiv, '/CouDiv2_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_CouDiv2[CouDiv_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouDiv2,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(min(shap_rast_CouDiv2[CouDiv_shap][], na.rm = T), max(shap_rast_CouDiv2[CouDiv_shap][], na.rm = T), length.out = 4),1)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()


# read in presence predictions for both species
pres_CouDiv1 <- rast(sp_raster_pres_pa[grepl(CouDiv1, sp_raster_pres_pa)])
pres_CouDiv2 <- rast(sp_raster_pres_pa[grepl(CouDiv2, sp_raster_pres_pa)])
pres_CouDiv12 <- mosaic(pres_CouDiv1, pres_CouDiv2)

shap_div_CouDiv1 <- shap_rast_CouDiv1[CouDiv_shap] > shap_rast_CouDiv2[CouDiv_shap] & shap_rast_CouDiv1[CouDiv_shap] > 0
shap_div_CouDiv2 <- shap_rast_CouDiv2[CouDiv_shap] > shap_rast_CouDiv1[CouDiv_shap] & shap_rast_CouDiv2[CouDiv_shap] > 0

shap_div_CouDiv1[shap_div_CouDiv1[]==F] <- NA
shap_div_CouDiv1[shap_div_CouDiv1[]==T] <- 1
shap_div_CouDiv2[shap_div_CouDiv2[]==F] <- NA
shap_div_CouDiv2[shap_div_CouDiv2[]==T] <- 2

shap_div_CouDiv12 <- mosaic(shap_div_CouDiv1, shap_div_CouDiv2)
pres_CouDiv12 <- resample(pres_CouDiv12, shap_div_CouDiv12)
shap_div_CouDiv12[is.na(pres_CouDiv12)] <- NA


## take average of species shapley values and plot
shap_mean_CouDiv12 <- app(c(shap_rast_CouDiv1[CouDiv_shap], shap_rast_CouDiv2[CouDiv_shap] ), mean)
pdf(paste0(CouDiv, 'mean_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_mean_CouDiv12) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = title_CouDiv_mean,
            legend.is.portrait = F, 
            legend.show = T, 
            midpoint = 0, 
            breaks = signif(seq(min(shap_mean_CouDiv12[], na.rm = T), max(shap_mean_CouDiv12[], na.rm = T), length.out = 4),1)
            #breaks = c(-0.2, -0.1, 0, 0.1, 0.2)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

## take the sd of species shapley values
shap_sd_CouDiv12 <- app(c(shap_rast_CouDiv1[CouDiv_shap], shap_rast_CouDiv2[CouDiv_shap] ), sd)
pdf(paste0(CouDiv, 'sd_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_sd_CouDiv12) +
  tm_raster(style = "cont",
            palette = c('gray90', '#bd0f06'),
            legend.reverse = F,
            title = title_CouDiv_sd,
            legend.is.portrait = F, 
            legend.show = T,
            breaks = signif(seq(min(shap_sd_CouDiv12[], na.rm = T), max(shap_sd_CouDiv12[], na.rm = T), length.out = 4),1)
            #breaks = c(0, 0.1, 0.2, 0.25)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# plot the environmental data
pdf(paste0(CouDiv, 'raw_map.pdf'), width = 3, height = 3, bg = 'transparent')
tm_shape(env_data[CouDiv_var]) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = title_CouDiv_raw, 
            legend.show = T, 
            midpoint = 0, 
            breaks = signif(seq(min(env_data[CouDiv_var][], na.rm = T), max(env_data[CouDiv_var][], na.rm = T), length.out = 4),1)) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

#### 6. Example 2. Contrast response curves for connnectivity across species ----

# contrast species
CouCon1 <- 'Lampetra planeri'
CouCon2 <- 'Barbus barbus'

CouCon1_lab <- 'L. planeri'
CouCon2_lab <- 'B. barbus'

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

shap_CouCon12 <- shap_CouCon12 %>% mutate(species_name_lab = case_when(.$species_name == CouCon1 ~ CouCon1_lab, 
                                                                       .$species_name == CouCon2 ~ CouCon2_lab))

# contrasting species response curves
CouCon <- paste0(contrast_dir, 'Coupled_Convergent/')
dir.create(CouCon, recursive =  T)

pdf(paste0(CouCon, 'response_curves.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_CouCon12 %>% sample_n(., nrow(.)), 
       aes(x = .data[[CouCon_var]], y = .data[[CouCon_shap]])) + 
  geom_jitter(aes(pch = species_name_lab, col = species_name_lab)) + 
  stat_smooth(aes(group = species_name_lab, lty = species_name_lab), 
              col = 'black',se = F, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4), 
              size = 0.5) + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.3,0.9), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank(), 
        text = element_text(size = 14)) + 
  ylab(shap_label) + 
  xlab(env_label) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_CouCon12_wide <- shap_CouCon12 %>% 
  select(-species_name_lab) %>% 
  pivot_wider(names_from = 'species_name', values_from = all_of(CouCon_shap))

shap_CouCon12_wide$col <- ifelse(shap_CouCon12_wide[[CouCon1]] > 0 & shap_CouCon12_wide[[CouCon2]] > 0, 1, 
                              ifelse(shap_CouCon12_wide[[CouCon1]] < 0 & shap_CouCon12_wide[[CouCon2]] < 0, 2, 3))


# set the limits
max_axis <- max(c(shap_CouCon12_wide[[CouCon2]], shap_CouCon12_wide[[CouCon1]]), na.rm = T)
min_axis <- min(c(shap_CouCon12_wide[[CouCon2]], shap_CouCon12_wide[[CouCon1]]), na.rm = T)

pdf(paste0(CouCon, 'shapley_biplot.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_CouCon12_wide, 
       aes(x = .data[[CouCon1]], y = .data[[CouCon2]], col = (.data[[CouCon1]] + .data[[CouCon2]]) / 2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25, col = 'gray50') + 
  geom_vline(aes(xintercept = 0), lwd = 0.25, col = 'gray50') + 
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
        plot.background = element_blank(),
        text = element_text(size = 14)) + 
  xlab(bquote(phi1~ .(paste0(' ', CouCon1_lab)))) + 
  ylab(bquote(phi1~ .(paste0(' ', CouCon2_lab)))) + 
  geom_abline(lwd = 0.25, col = 'gray50') + 
  xlim(c(min_axis, max_axis)) + 
  ylim(c(min_axis, max_axis))
dev.off()

pdf(paste0(CouCon, '/CouCon1_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_CouCon1[CouCon_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouCon1,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 4), 1)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(CouCon, '/CouCon2_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_CouCon2[CouCon_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = CouCon2,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 4), 1)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# read in presence predictions for both species
pres_CouCon1 <- rast(sp_raster_pres_pa[grepl(CouCon1, sp_raster_pres_pa)])
pres_CouCon2 <- rast(sp_raster_pres_pa[grepl(CouCon2, sp_raster_pres_pa)])
pres_CouCon12 <- mosaic(pres_CouCon1, pres_CouCon2)


# plot the environmental data
pdf(paste0(CouCon, 'raw_map.pdf'), width = 3, height = 3, bg = 'transparent')
tm_shape(env_data['local_asym_cl_log10']) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = 'Connectivity', 
            legend.show = T, 
            breaks = signif(seq(from = min(env_data['local_asym_cl_log10'][], na.rm = T), to = max(env_data['local_asym_cl_log10'][], na.rm = T), length.out = 4),3)) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()



#### 7. Example 3. Contrast response curves for decoupled effects ----

# contrast species
DeCou1 <- 'Cottus gobio'
DeCou2 <- 'Gobio gobio'

# variable name
DeCou_var <- 'local_imd_log10_ele_residual'
# shapley name
DeCou_shap <- 'local_imd_log10_ele_residual_SHAP'

# environmental label
env_label = 'urbanisation'
shap_label = expression(phi1~'urbanisation')

# plot titles
title_DeCou_mean = expression('Mean urbanisation' ~ phi1)
title_DeCou_sd   = expression('S.D. urbanisation' ~ phi1)
title_DeCou_raw  = 'Raw urbanisation'


# read in rasters for focal variables
shap_rast_DeCou1 <- rast(paste0(shap_dirs[grepl(DeCou1, shap_dirs)], "/shap_raster_pa.TIF"))
shap_rast_DeCou2 <- rast(paste0(shap_dirs[grepl(DeCou2, shap_dirs)], "/shap_raster_pa.TIF"))

# read in raw shapleys
shap_DeCou1 <- readRDS(shap_pa[grepl(DeCou1, shap_pa)])
shap_DeCou2 <- readRDS(shap_pa[grepl(DeCou2, shap_pa)])

# join together shapley values to environmental data for each TEILEZGNR
shap_DeCou1 <- left_join(shap_DeCou1 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
shap_DeCou2 <- left_join(shap_DeCou2 %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get variable of interest
shap_DeCou12 <- bind_rows(shap_DeCou1, shap_DeCou2) %>% 
  select(species_name, TEILEZGNR, matches(DeCou_var))

# contrasting species response curves
DeCou <- paste0(contrast_dir, 'Decoupled/')
dir.create(DeCou, recursive =  T)

pdf(paste0(DeCou, 'response_curves.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_DeCou12 %>% sample_n(., nrow(.)), 
       aes(x = .data[[DeCou_var]], y = .data[[DeCou_shap]])) + 
  geom_jitter(aes(pch = species_name, col = species_name)) + 
  stat_smooth(aes(group = species_name, lty = species_name), 
              col = 'black',se = F, method = 'lm', formula = y ~ x + I(x^2) + I(x^3) + I(x^4), 
              size = 0.5) + 
  geom_hline(aes(yintercept = 0), lwd = 0.25) + 
  scale_colour_manual(values = c('gray50', 'gray75')) + 
  theme_bw() +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        rect = element_rect(fill = "transparent"), 
        legend.position = c(0.6,0.9), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank(), 
        text = element_text(size = 14)) + 
  ylab(shap_label) + 
  xlab(env_label) + 
  guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

# make wide for plot of shapley values for each species against eachother
shap_DeCou12_wide <- shap_DeCou12 %>% 
  pivot_wider(names_from = 'species_name', values_from = all_of(DeCou_shap))

shap_DeCou12_wide$col <- ifelse(shap_DeCou12_wide[[DeCou1]] > 0 & shap_DeCou12_wide[[DeCou2]] > 0, 1, 
                                 ifelse(shap_DeCou12_wide[[DeCou1]] < 0 & shap_DeCou12_wide[[DeCou2]] < 0, 2, 3))


# set the limits
max_axis <- max(c(shap_DeCou12_wide[[DeCou2]], shap_DeCou12_wide[[DeCou1]]), na.rm = T)
min_axis <- min(c(shap_DeCou12_wide[[DeCou2]], shap_DeCou12_wide[[DeCou1]]), na.rm = T)

pdf(paste0(DeCou, 'shapley_biplot.pdf'), width = 3, height = 3, bg = 'transparent')
ggplot(data = shap_DeCou12_wide, 
       aes(x = .data[[DeCou1]], y = .data[[DeCou2]], col = (.data[[DeCou1]] + .data[[DeCou2]]) / 2)) + 
  geom_jitter() + 
  geom_hline(aes(yintercept = 0), lwd = 0.25, col = 'gray50') + 
  geom_vline(aes(xintercept = 0), lwd = 0.25, col = 'gray50') + 
  scale_colour_gradientn(
    colors=c('#bd0f06','gray90', '#2200c9'),
    values=scales::rescale(c(-0.1, -0.02, 0, 0.02 , 0.1)),
    limits=c(-max(abs(range(shap_DeCou12_wide[c(DeCou1, DeCou2)], na.rm = T))), 
             max(abs(range(shap_DeCou12_wide[c(DeCou1, DeCou2)], na.rm = T))))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        aspect.ratio = 0.5,  
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        plot.background = element_blank(), 
        text = element_text(size = 14)) + 
  xlab(bquote(phi1~ .(paste0(' ', DeCou1)))) + 
  ylab(bquote(phi1~ .(paste0(' ', DeCou2)))) + 
  geom_abline(lwd = 0.25, col = 'gray50') + 
  xlim(c(min_axis, max_axis)) + 
  ylim(c(min_axis, max_axis))
dev.off()

pdf(paste0(DeCou, '/DeCou1_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_DeCou1[DeCou_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = DeCou1,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 4), 1)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

pdf(paste0(DeCou, '/DeCou2_map.pdf'), width = 4, height =4, bg = 'transparent')
tm_shape(shap_rast_DeCou2[DeCou_shap]) +
  tm_raster(style = "cont",
            palette = c('#bd0f06','gray90', '#2200c9'),
            legend.reverse = F,
            title = DeCou2,
            legend.is.portrait = F, 
            legend.show = F, 
            midpoint = 0,
            breaks = signif(seq(from = signif(min_axis,2), to = signif(max_axis,2), length.out = 4), 1)
  ) + 
  tm_layout(bg.color = "transparent", 
            frame = F)
dev.off()

# read in presence predictions for both species
pres_DeCou1 <- rast(sp_raster_pres_pa[grepl(DeCou1, sp_raster_pres_pa)])
pres_DeCou2 <- rast(sp_raster_pres_pa[grepl(DeCou2, sp_raster_pres_pa)])
pres_DeCou12 <- mosaic(pres_DeCou1, pres_DeCou2)


# plot the environmental data
pdf(paste0(DeCou, 'raw_map.pdf'), width = 3, height = 3, bg = 'transparent')
tm_shape(env_data[DeCou_var]) + 
  tm_raster(style = 'cont', 
            palette = rev(c('#bd0f06','gray90', '#2200c9')), 
            legend.is.portrait = F, 
            title = title_DeCou_raw, 
            legend.show = T, 
            midpoint = 0, 
            breaks = signif(seq(min(env_data[DeCou_var][], na.rm = T), max(env_data[DeCou_var][], na.rm = T), length.out = 4),1)) + 
  tm_layout(bg.color = "transparent", 
            frame = F, 
            legend.hist.width = 0.1)
dev.off()





#### 8. Catchment force plots ----

# get the teilenzugsgebeit for the sense
sense <- data.frame(catchment = 'Sense', 
                    TEZGNR40 = c(100184, 101484, 105044, 105663, 102204, 102377))
emme <- data.frame(catchment = 'Emme', 
                   TEZGNR40 = c(106570, 106299, 101728, 101280, 106548, 103982))
catchments <- rbind(sense, emme)

# testing taking the 2km2 subcatchments containig the main stems
sense <- data.frame(catchment = 'Sense', 
                    TEILEZGNR = c(54206, 49238, 36039, 53576, 92800, 
                                  73564, 79104, 55745, 61040, 21140, 
                                  76158, 40046, 95372, 78251, 21268, 
                                  96970, 1492, 64847, 85482, 21499, 
                                  5407, 72778, 98155, 18709, 9298))
emme <- data.frame(catchment = 'Emme', 
                    TEILEZGNR = c(48935, 519, 76613, 38042, 23267,64482, 
                                  17698, 39457, 59083, 8682, 42744, 
                                  14206, 88728, 80592, 44084, 92441,
                                  62780, 47744, 19469, 80840, 66777, 
                                  54624, 54833, 26718, 81132, 31679, 
                                  4435, 94151))
catchments <- rbind(sense, emme)


# read in subcatchments, transform, union
subcatchments_sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% catchments$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

# get non-unioned objects to explore plots
subcatchments_sense_nounion <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% catchments$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  st_zm()
mapview::mapview(subcatchments_sense_nounion, zcol = 'TEILEZGNR')

# get the areas of the focal regions
st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% catchments$TEILEZGNR) %>% 
  left_join(., catchments) %>% 
  group_by(catchment) %>% 
  do(area = sum(st_area(.))) %>%  
  unnest()

# join with focal catchments
# subcatchment_2km_names <- left_join(subcatchment_2km_names, catchments)

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
shapley_focal <- left_join(shapley_focal, catchments)


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


shapley_catchment_summary %>% 
  filter(species_name == 'Perca fluviatilis', 
         catchment == 'Sense') %>% View


# quick example plot to see if the data are in the correct format
ggplot(data = shapley_catchment_summary) + 
  geom_point(aes(y = species_name, x = (baseline_prediction)), col = 'black') + 
  geom_point(aes(y = species_name, x = (catchment_prediction)), col = 'blue') + 
  facet_wrap(~catchment)

# join in renaming object
shapley_catchment_summary <- left_join(shapley_catchment_summary, 
                                       data.frame(vars_renamed), 
                                       by = c('name' = 'vars_shap'))

# order by mean value
levels <- shapley_catchment_summary %>% 
  group_by(vars_renamed) %>% 
  summarise(mean_shap2 = mean(mean_shap, na.rm = T)) %>% 
  arrange(mean_shap2) %>% 
  pull(vars_renamed)

shapley_catchment_summary$vars_renamed <- factor(shapley_catchment_summary$vars_renamed, 
                                                 levels = levels)


## Read in raw data used to fit the models
pa_data <- read.csv(paste0(dd, '/sdm-pipeline/species-records-final/fish-presenceAbsence_2010.csv'))

# take only presences
pa_data_2 <- pa_data %>% 
  filter(occ == 1, 
         species_name %in% sp_list) %>% 
  st_as_sf(., coords = c('X', 'Y'), crs = 'wgs84') %>% 
  # transform to ETRS89
  st_transform(., crs = 'epsg:3035')

# get the raw catchments
all_catchments <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% catchments$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  select(TEILEZGNR)

# join in catchment data
pa_data_catchments <- st_join(pa_data_2, all_catchments) %>% 
  left_join(., catchments)

# create column to colour by
pa_data_catchment_simple <- pa_data_catchments %>% 
  filter(!is.na(catchment)) %>% 
  select(species_name, catchment) %>% 
  st_drop_geometry() %>% 
  unique() %>% 
  mutate(present_in_sense = ifelse(.$catchment == 'Sense', 1, 0), 
         present_in_emme  = ifelse(.$catchment == 'Emme', 1, 0))

# join with summary species
shapley_catchment_summary <- left_join(shapley_catchment_summary, pa_data_catchment_simple) %>% 
  mutate(present_in_sense = ifelse(is.na(.$present_in_sense), 0, .$present_in_sense), 
         present_in_emme = ifelse(is.na(.$present_in_emme), 0, .$present_in_emme))
  

# make output directories
force_dir <- paste0(fig_dir, 'local_force_plots')
dir.create(paste0(force_dir), recursive = T)

shapley_catchment_summary$species_name <- plyr::revalue(shapley_catchment_summary$species_name, 
                                                        c('Alburnoides bipunctatus' = 'A. bipunctatus', 
                                                          'Barbus barbus' = 'B. barbus', 
                                                          'Cottus gobio' = 'C. gobio', 
                                                          'Gobio gobio' = 'G. gobio', 
                                                          'Lampetra planeri' = 'L. planeri', 
                                                          'Oncorhynchus mykiss' = 'O. mykiss', 
                                                          'Perca fluviatilis' = 'P. fluviatilis', 
                                                          'Squalius cephalus' = 'S. cephalus', 
                                                          'Thymallus thymallus' = 'T. thymallus'))


# remove wetland due to generally weak effects
shapley_catchment_summary <- shapley_catchment_summary %>% filter(vars_renamed != 'wetland')

# check directions of effects
pos <- '#528F70'
neg <- '#FF8000'
direction <- rev(c(pos, neg, pos, neg, neg, neg, pos, pos, neg, pos))

# Set up colours for direction plot
shapley_catchment_summary <- shapley_catchment_summary %>% 
  mutate(bar_cols = case_when(
    sign_shap == 1  & present_in_emme == 0 & catchment == 'Emme' ~ '#9287cc', 
    sign_shap == -1 & present_in_emme == 0 & catchment == 'Emme' ~ '#e89e99', 
    sign_shap == 1  & present_in_emme == 1 & catchment == 'Emme' ~ '#2200c9', 
    sign_shap == -1 & present_in_emme == 1 & catchment == 'Emme' ~ '#bd0f06', 
    sign_shap == 1  & present_in_sense == 0 & catchment == 'Sense' ~ '#9287cc', 
    sign_shap == -1 & present_in_sense == 0 & catchment == 'Sense' ~ '#e89e99', 
    sign_shap == 1  & present_in_sense == 1 & catchment == 'Sense' ~ '#2200c9', 
    sign_shap == -1 & present_in_sense == 1 & catchment == 'Sense' ~ '#bd0f06'))


# make force plots contrasting mean responses across emme and sense
pdf(paste0(force_dir, '/emme_sense_force.pdf'), width = 10, height = 5)
ggplot(data = shapley_catchment_summary) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'none', 
        axis.line.x = element_line(), 
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), 
        axis.ticks = element_line(), 
        strip.text.x = element_text(face = 'italic'), 
        strip.text.y = element_text(face = 'bold', size = 12),
        #panel.border = element_rect(fill = 'transparent'), 
        aspect.ratio = 1.75, 
        text = element_text(colour = 'black'), 
        panel.border = element_blank(), 
        panel.spacing.y = unit(2, "lines"), 
        panel.spacing.x = unit(1, "lines"), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(colour = direction)) + 
  
  facet_grid(catchment~species_name, scales = 'free_x') + 
  
  #draw segments
  geom_segment(data = shapley_catchment_summary %>% 
                 filter(catchment == 'Emme'), 
               aes(xend = mean_shap + baseline_prediction, 
                   x = baseline_prediction, 
                   y = vars_renamed, yend = vars_renamed, 
                   col = bar_cols
                   ),
               lwd = 4, lineend = 'butt', linejoin = 'mitre') + 
  geom_segment(data = shapley_catchment_summary %>% 
                 filter(catchment == 'Sense'), 
               aes(xend = mean_shap + baseline_prediction, 
                   x = baseline_prediction, 
                   y = vars_renamed, yend = vars_renamed, 
                   col = bar_cols
                   ),
               lwd = 4, lineend = 'butt', linejoin = 'mitre') + 
  
  # draw vertical lines
  geom_vline(data = shapley_catchment_summary %>% 
               filter(catchment == 'Emme'), 
             aes(xintercept = baseline_prediction), size = 0.2) + 
  geom_vline(data = shapley_catchment_summary %>% 
               filter(catchment == 'Sense'), 
             aes(xintercept = baseline_prediction), size = 0.2) + 
  
 #  # draw points
 #  geom_point(data = shapley_catchment_summary %>% 
 #               filter(catchment == 'Emme'), 
 #             aes(x = mean_shap + baseline_prediction, 
 #                 y = vars_renamed, 
 #                 col = bar_cols, 
 #                 size = abs(mean_shap)*2)) + 
 #  geom_point(data = shapley_catchment_summary %>% 
 #               filter(catchment == 'Sense'), 
 #             aes(x = mean_shap + baseline_prediction, 
 #                 y = vars_renamed, 
 #                 col = bar_cols, 
 #                 size = abs(mean_shap)*2)) + 
  
  
  # scales
  scale_colour_manual(values = na.omit(unique(shapley_catchment_summary$bar_cols))[c(3,1,4,2)]) + 
  xlab('Change in local prediction from average prediction') + 
  ylab(NULL) +
  coord_cartesian(clip = "off")
dev.off()

#e3b1af

#### 9. Make catchment specific response curves ----

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

#### 10. Plot environmental data for each catchment ----

subcatchment_file <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/swiss-2km-subcatchments/EZG_Gewaesser.gdb'

# read in subcatchments, transform, union
emme_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% emme$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_as_sf() %>% 
  mutate(name = 'Emme')

sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEILEZGNR %in% sense$TEILEZGNR) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm() %>% 
  st_as_sf() %>% 
  mutate(name = 'Sense')

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

vars_renamed = c('discharge', 
                 'slope', 
                 'flow velocity', 
                 'temperature max', 
                 'temperature min', 
                 'connectivity', 
                 'distance to lake', 
                 'morph. mod.', 
                 'urbanisation',
                 'wetland', 
                 'floodplains')

# reset environmental names
vars_renamed_V2 <- cbind(vars, vars_renamed)

emme_env <- extract(env_data, vect(emme_union)) %>% 
  select(all_of(vars)) %>% 
  pivot_longer(colnames(.)) %>% 
  mutate(catchment = 'Emme')

sense_env <- extract(env_data, vect(sense_union)) %>% 
  select(all_of(vars)) %>% 
  pivot_longer(colnames(.)) %>% 
  mutate(catchment = 'Sense')

catchment_envs <- rbind(emme_env, sense_env)

# organise variable names
catchment_envs <- left_join(catchment_envs, data.frame(vars_renamed_V2), by = c('name' = 'vars'))

# take sample of all 
env_df <- data.frame(env_data) %>% na.omit 

sample_env_data <-  env_df %>% 
  filter(elevation < 900) %>% 
  sample_n(., 10000) %>% 
  pivot_longer(colnames(.)) %>% 
  mutate(catchment = 'all')
sample_env_data <- left_join(sample_env_data, data.frame(vars_renamed_V2), by = c('name' = 'vars'))
sample_env_data <- na.omit(sample_env_data)

pdf(paste0(force_dir, '/emme_sense_env.pdf'), width = 9)
ggplot(data = catchment_envs) + 
  geom_boxplot(aes(x = catchment, y = value, col = catchment)) + 
  geom_boxplot(data = sample_env_data,
               aes(x = catchment, y = value, col = catchment)) + 
  facet_wrap(~vars_renamed, scale = 'free') + 
  theme_bw() + 
  theme(panel.grid = element_blank())
dev.off()

