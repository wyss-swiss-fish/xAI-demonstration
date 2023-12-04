### Figure 2. Explore all Shapley values for a given species

### Alburnoides bipunctatus maps of spatial distribution of key variables

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

# get run to mak figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- c("Alburnoides bipunctatus")

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
                        full.names = T
)
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


#### 4. Estimate variable importance and variable effect distribution ----

# read in shapley values for presence absence models
shap_rast_pa_i <- rast(paste0(shap_dirs, "/shap_raster_pa.TIF"))

# clean variable names
names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']

# remove non_na catchments
shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]

# calculate average variable importance across the entire landscape from shapley values
all_shap_df <- data.frame(shap_rast_pa_i)
var_imp <- sapply(all_shap_df, function(x) mean(abs(x), na.rm = T))
var_imp_raw  <- sapply(all_shap_df, function(x) mean(x, na.rm = T))
var_imp <- data.frame(var_imp, var_imp_raw, var_names = stringr::str_replace(stringr::str_replace(names(var_imp), '\\.', ' '), '\\.', ' '))
var_imp$var_names <- factor(var_imp$var_names, levels = var_imp$var_names[order(var_imp$var_imp)])

print(tmaptools::get_brewer_pal("YlOrRd", n = 7))

dir.create(paste0(fig_dir, '/spatial_shap_example/'), recursive = T)

pdf(paste0(fig_dir, '/spatial_shap_example/variable_importance_', sp_list[1], '.pdf'), width = 4, height = 4)
ggplot(data = var_imp) +
  geom_bar(aes(x = var_imp, y = var_names, fill = var_imp), stat = 'identity') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.position = 'none', aspect.ratio = 1.5, 
        axis.text.y = element_text(colour = 'black', size = 10, vjust = 0.2)) + 
  scale_fill_gradient(low = '#FFF7BA', high = '#DF171C') + 
  ylab(NULL) + 
  xlab('mean(|Shapley values|)') + 
  ggtitle('(b) variable importance')
dev.off()


# pivot shapley values longer for creating full distribution of values
all_shap_long <- all_shap_df %>% pivot_longer(cols = names(.))
all_shap_long$name <- gsub('\\.', ' ', all_shap_long$name)
all_shap_long <- left_join(all_shap_long, var_imp, by = c('name' = 'var_names'))
all_shap_long$name <- factor(all_shap_long$name, levels = levels(var_imp$var_names))


pdf(paste0(fig_dir, '/spatial_shap_example/shap_distribution_', sp_list[1], '.pdf'), width = 4, height = 4)
ggplot(data = all_shap_long) + 
  geom_jitter(data = all_shap_long %>% group_by(name) %>% sample_n(., 5000), 
              aes(x = value, y = name, col = var_imp), pch = 19, alpha = 0.5, stroke = NA, size = 2) + 
  geom_violin(data = all_shap_long %>% group_by(name), 
              aes(x = value, y = name), fill = 'transparent', scale = 'width') + 
  scale_fill_gradient(low = '#FFF7BA', high = '#DF171C') + 
  scale_colour_gradient(low = '#FFF7BA', high = '#DF171C') + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.position = 'none', 
        aspect.ratio = 1.5, 
        axis.text.y = element_text(colour = 'black', size = 10, vjust = 0.2)) + 
  ylab(NULL) + 
  xlab('Shapley values') + 
  ggtitle('(c) variable effect distribution') + 
  geom_vline(aes(xintercept=0))
dev.off()


#### 5. Map spatial distribution of shapley values for Alburnoides bipunctatus ---- 

# read in shapley values for presence-absence models
shap_rast_pa_i <- rast(paste0(shap_dirs, "/shap_raster_pa.TIF"))

# clean variable names
names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']

# remove non_na
shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]

# get mean values to order by
shap_rast_pa_i <- shap_rast_pa_i[[order(var_imp$var_imp, decreasing = T)]]

# get raster of species presences and absences
suit_i <- rast(sp_raster_suit_pa)
presence <- rast(sp_raster_pres_pa[1])
suit_i <- resample(suit_i, presence)

# convert non-predicted locations to NAs
pres_suit <- suit_i
pres_suit[is.na(presence)] <- NA
names(pres_suit) <- 'presence'

# convert non-predicted locations to NAs
abs_suit <- suit_i
abs_suit[!is.na(presence)] <- NA
names(abs_suit) <- 'absence'

# combine presences and absences
pa_map <- c(pres_suit, abs_suit)
pa_map <- trim(pa_map)

# make plot of locations predicted as present
occ_plot <- 
  tm_shape(ch) + 
  tm_borders(col = 'black') + 
  tm_shape(suit_i) + 
  tm_raster(style = 'cont', 
            legend.is.portrait = F, 
            legend.show = F, 
            palette = 'YlOrRd') + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, 
           col = 'black', 
           lwd = 0.1) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = T,
            legend.outside.size = 0.1,
            legend.outside.position = c('bottom'), 
            panel.labels=c('(a) relative occurrence'), 
            panel.label.bg.color = 'white', 
            frame = FALSE, 
            bg.color = "transparent") + 
  tm_shape(ch) + 
  tm_borders(col = 'black') + 
  tm_facets(ncol = 1) +
  tm_scale_bar(width  = 0.25, 
               text.size = 0.5, 
               position = c('center', 'bottom'))
  
  
# plot shapley values as a faceted tmap  
shap_rast_plot <- shap_rast_pa_i
names(shap_rast_plot) <- paste0('(', letters[3:(2+nlyr(shap_rast_plot))], ') ', names(shap_rast_plot))
shap_plot <- tm_shape(shap_rast_plot) +
  tm_raster(
    style = "cont",
    palette = c('#bd0f06','gray90', '#2200c9'), 
    midpoint = 0,
    legend.reverse = F,
    title = '',
    legend.is.portrait = F, 
    legend.show = F
  ) +
  tm_facets(ncol = 4, 
            free.scales.raster = T) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, 
           col = 'black', 
           lwd = 0.1) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = F, 
            legend.outside.size = 0.1,
            legend.outside.position = c('bottom'), 
            panel.labels=names(shap_rast_plot), 
            panel.label.bg.color = 'white', 
            frame = FALSE,
            bg.color = "transparent", 
            panel.label.size = 1.2)


pdf(paste0(fig_dir, '/spatial_shap_example/presence_', sp_list[1], '.pdf'), width = 6*1.5, height = 2*1.5, 
    bg = 'transparent')
print(occ_plot) 
dev.off()

pdf(paste0(fig_dir, '/spatial_shap_example/presence_legend', sp_list[1], '.pdf'), width = 2, height = 2, 
    bg = 'transparent')
print(tm_shape(pa_map[[1]]) + 
        tm_raster(style = 'cont', 
                  title = '',
                  legend.is.portrait = F, 
                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                  labels = as.character(c(0, 0.25, 0.5, 0.75, 1))) +
        tm_layout(legend.only = T)) 
dev.off()


# make pdf of shapley distributions
pdf(paste0(fig_dir, '/spatial_shap_example/shap_', sp_list[1], '.pdf'), width = 8, height = 6)
print(shap_plot) 
dev.off()

# make png to save space
png(paste0(fig_dir, '/spatial_shap_example/shap_', sp_list[1], '.png'), width = 3500, height = 1800, res = 300)
print(shap_plot) 
dev.off()

# pdf of legend
pdf(paste0(fig_dir, '/spatial_shap_example/shap_legend_', sp_list[1], '.pdf'), width = 2, height = 2)
print(tm_shape(shap_rast_plot[[1]]) +
  tm_raster(
    style = "cont",
    palette = c('#bd0f06','gray90', '#2200c9'), 
    midpoint = 0,
    legend.reverse = FALSE,
    title = 'Shapley values',
    legend.is.portrait = F) + 
  tm_layout(legend.only = T))
dev.off()


#### 6. Make response curves per variable ----

# read in raw shapleys
shap_i <- readRDS(shap_pa)
# join together shapley values to environmental data for each TEILEZGNR
shap_i <- left_join(shap_i %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')

# get the shapley values
grouped_rc_shap <- shap_i %>% 
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
  data.frame %>% 
  arrange(vars_renamed)

# get the minimum and maximum
min_y <- min(grouped_rc_shap$shapley, na.rm = T)
max_y <- max(grouped_rc_shap$shapley, na.rm = T)

grouped_rc_shap$vars_renamed <- factor(grouped_rc_shap$vars_renamed, 
                                       levels = var_imp$var_names[order(var_imp$var_imp, decreasing = T)])


# make shapley plot for all variables for a single species
dir.create(paste0(fig_dir, 'spatial_shap_example/response-curves/'), recursive = T)

for(i in 1:length(unique(grouped_rc_shap$vars_renamed))){
  var_i <- unique(grouped_rc_shap$vars_renamed)[i]
  
  data_i <- grouped_rc_shap %>% 
    filter(vars_renamed == var_i) %>% 
    sample_n(5000)
  
pdf(paste0(fig_dir, 'spatial_shap_example/response-curves/', var_i, '.pdf'), height = 1.5, width = 3, bg = 'transparent')
  print(ggplot(data = data_i, 
       aes(x=env, y = shapley, col = shapley)) + 
    geom_point() +
  stat_smooth(col = 'black', se = F) +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.5, 
        axis.title = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_blank(), 
        legend.position = 'none', 
        panel.background = element_blank(), 
        plot.background = element_blank()) +
  scale_colour_gradient2(low  = '#bd0f06',
                         mid  = 'gray90',
                         high = '#2200c9',
                         midpoint = 0))
  dev.off()
  
  }


#### 7. Correlation between shapley values ----
library(GGally)
library(corrplot)
# create unique wide matrix
unique_shap_wide <- unique(all_shap_df)

# change names
names(unique_shap_wide) <- gsub('\\.', '  ', names(unique_shap_wide))

# create correlation matrix
cor_shap <- cor(unique_shap_wide, method = 'spearman')
upr_tri_shap <- cor_shap[upper.tri(cor_shap)]
median(abs(upr_tri_shap),na.rm = T)
summary(abs(upr_tri_shap),na.rm = T)
sd(abs(upr_tri_shap),na.rm = T)

# example correlation for text beween discharge and connectivity
cor_discharge_connectivity <- signif(cor_shap['discharge', 'connectivity'], 1)
cor_median_abs_all <- signif(median(abs(upr_tri_shap),na.rm = T),2)
cor_mean_sd_all <- paste0(signif(mean(abs(upr_tri_shap),na.rm = T),2),' ± ',signif(sd(abs(upr_tri_shap),na.rm = T),2))


# plot correlation matrix
pdf(paste0(fig_dir, 'spatial_shap_example/', 'shap_correlations.pdf'), width = 7, height = 7)
ggcorr(unique_shap_wide, 
       method = c("everything", "spearman"), 
       limits = c(-1, 1), 
       label = T, 
       label_round = 2, 
       angle = 0, 
       hjust = 0.75, 
       low = "#bd0f06",
       mid = "gray90",
       high = "#2200c9")
dev.off()


#### 8. Estimate area of threatened natural range ----

# custom function to estimate shadow distribution properties and output as a shapefile
source('scripts/functions/generate_SD_outputs.R')

# define natural niche factors
natural_niche_factors = c('ecoF_discharge_max_log10_SHAP', 
                          'stars_t_mn_m_c_SHAP', 
                          'stars_t_mx_m_c_SHAP', 
                          'ecoF_flow_velocity_mean_SHAP', 
                          'local_dis2lake_SHAP', 
                          'ecoF_slope_min_log10_SHAP')
# define habitat 'threat' factors
habitat_factors = c('local_wet_SHAP', 
                    'local_flood_SHAP', 
                    'local_imd_log10_ele_residual_SHAP', 
                    'ecoF_eco_mean_ele_residual_SHAP')
# define connectivity 'threat' factors
conn_factors = 'local_asym_cl_log10_SHAP'

### Apply function that estimates the properties of the shadow distributions
sp_shap <- generate_SD_outputs(sdm_input_data = sdm_dirs[1], 
                               raster_data = sp_raster_suit_pa, 
                               shap = shap_pa, 
                               natural_niche_factors = natural_niche_factors,
                               habitat_factors = habitat_factors,
                               conn_factors = conn_factors,
                               species = sp_list[1],
                               output_folder = paste0('figures/ubelix_SDM_RF_APRIL_V1_02/shadow_dist_summaries/'))

AB_shad_summary <- read.csv('figures/ubelix_SDM_RF_APRIL_V1_02/shadow_dist_summaries/Alburnoides bipunctatus.csv')


# create scatter plot of shapley values against suitability classified by threat
pdf(paste0(fig_dir, '/spatial_shap_example/natural_threat_suitability_biplot.pdf'), 
    width = 3, height = 3, bg = 'transparent')
ggplot() + 
  geom_point(data = sp_shap %>% arrange(niche_categories), 
             aes(y = suitability, 
                 x = natural_niche_value, 
                 col = niche_categories, 
                 size = niche_categories,
                 alpha = is.na(presence)), 
             stroke = 0) + 
  xlab(expression('abiotic niche' ~phi~ 'values')) + 
  ylab('predicted habitat suitability') + 
  scale_colour_manual('', values = c('#E8E8E8','#06C4C4','#D1D111', '#D16111', '#F70000')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none', 
        aspect.ratio = 1, 
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line()) + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  scale_size_manual(values = c(1,2,2,2,2)) + 
  scale_alpha_manual(values = c(1,0.4))
dev.off()


#### 9. Plot areas that are inside the ecological niche but are threatened ----

## 1 - outside of ecological niche
## 2 - inside of ecological niche and few threats
## 3 - natural areas with connectivity negative
## 4 - natural areas with flood negative
## 5 - natural areas with both connectivity and floodplains negative
threat_niche <- tm_shape(sp_shap) + 
  tm_fill(col = 'niche_categories', 
          style = 'fixed',
          palette = c('#E8E8E8','#06C4C4','#D1D111', '#D16111', '#F70000'), 
          title = '') +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F) + 
  tm_layout(frame = F,
            bg.color = "transparent")

# map the threat maps
pdf(paste0(fig_dir, '/spatial_shap_example/natural_threat_range.pdf'))
threat_niche
dev.off()


#### 10. Look at natural factor constrains inside and outside of species predicted areas ----

###
# find the catchments where natural factors act positively and identify the factor
# with the highest positive effect
which_pos <- st_drop_geometry(sp_shap)[,which(names(sp_shap) %in% natural_niche_factors)]
which_pos[which_pos<0] <- NA
sp_shap$which_pos <- unlist(apply(which_pos, 1, function(x) ifelse(length(which.max(x))==0, NA, names(which_pos)[which.max(x)])))

# subset to catchments inside the natural niche of the species and rename
sp_shap_which_pos <- sp_shap %>% 
  filter(natural_niche == T) %>% 
  left_join(., data.frame(vars_renamed), by = c('which_pos' = 'vars_shap'))

###
# find the catchments where natural factors act negatively, and identify the facto
# with the lowest negative effect
which_neg <- st_drop_geometry(sp_shap)[,which(names(sp_shap) %in% natural_niche_factors)]
which_neg[which_neg>0] <- NA
sp_shap$which_neg <- unlist(apply(which_neg, 1, function(x) ifelse(length(which.max(x))==0, NA, names(which_neg)[which.min(x)])))

# subset to catchments inside the natural niche of the species and rename
sp_shap_which_neg <- sp_shap %>% 
  filter(natural_niche == F) %>% 
  left_join(., data.frame(vars_renamed), by = c('which_neg' = 'vars_shap'))


### make plots of distributed niche effects
# plot positive niche effects
natural_niche_pos <- tm_shape(sp_shap_which_pos) +
  tm_fill(col = 'vars_renamed', 
          style = 'fixed',
            palette = c('#FAE100', '#59CAFA', '#845CFA', '#FA522A'), 
            title = 'Main positive niche effect') +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.5) + 
  tm_layout(frame = F,
            bg.color = "transparent")

# plot negative niche effects
natural_niche_neg <- tm_shape(sp_shap_which_neg) +
  tm_fill(col = 'vars_renamed', 
          style = 'fixed',
          palette = c('#FAE100', '#59CAFA', '#845CFA', '#FA522A'), 
            title = 'Main negative niche effect') +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.5) + 
  tm_layout(frame = F,
            bg.color = "transparent")

# combine into joint plot for manuscript figure
pdf(paste0(fig_dir, '/spatial_shap_example/natural_range_pos_neg.pdf'), height = 4*0.6, width = 14*0.6)
tmap_arrange(natural_niche_pos, natural_niche_neg, ncol = 2)
dev.off()



#### 11. Plot quantitative shadow distribution for A. bipunctatus ---- 

# define factor categories as global objects
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

source('scripts/functions/generate_quant_SD.R')
method = 3
ab_shad <- generate_quant_SD(sp_shap)

pdf(paste0(fig_dir, '/spatial_shap_example/quantitative_shadow.pdf'), height = 4*0.6, width = 14*0.6)
tm_shape(ab_shad %>% filter(!is.na(SD_OratioE))) +
  tm_fill(col = 'SD_OratioE', 
          style = 'cont',
          palette = 'Spectral',
          title = 'Quantitative shadow distribution', 
          colorNA = 'transparent', 
          #breaks = c(0, 0.5, 1), 
          legend.is.portrait = F) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.5) + 
  tm_layout(frame = F,
            bg.color = "transparent")
dev.off()


## save as an object
saveRDS(ab_shad %>% select(TEILEZGNR, nn_sum_mask, nn_sum, 
                           expected_distribution, observed_distribution, 
                           SD_OratioE) %>% 
          st_drop_geometry(), 
        file = 'alburnoides_SD.rds')

#### Results reporting ---- 

save(
  # correlations between shapley values and environment
  cor_discharge_connectivity,
  cor_median_abs_all, 
  cor_mean_sd_all, 
  # summary of shadow distributions for text
  AB_shad_summary,
  file = 'data/03-results-text.RData')

