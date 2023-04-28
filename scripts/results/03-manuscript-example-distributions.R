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
fig_dir <- "figures/ubelix_SDM_RF_MARCH_v6/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_MARCH_v6"

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


#### 4. Estimate variable importance and effect distribution ----

# read in presence absence data
shap_rast_pa_i <- rast(paste0(shap_dirs, "/shap_raster_pa.TIF"))

# clean names
names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']

# remove non_na
shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]

# varible importance across entire landscape
all_shap_df <- data.frame(shap_rast_pa_i)
var_imp <- sapply(all_shap_df, function(x) mean(abs(x), na.rm = T))
var_imp_raw  <- sapply(all_shap_df, function(x) mean(x, na.rm = T))
var_imp <- data.frame(var_imp, var_imp_raw, var_names = stringr::str_replace(stringr::str_replace(names(var_imp), '\\.', ' '), '\\.', ' '))
var_imp$var_names <- factor(var_imp$var_names, levels = var_imp$var_names[order(var_imp$var_imp)])

print(tmaptools::get_brewer_pal("YlOrRd", n = 7))

pdf(paste0(fig_dir, '/spatial_shap_example/variable_importance', sp_list[1], '.pdf'), width = 4, height = 4)
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

pdf(paste0(fig_dir, '/spatial_shap_example/shap_distribution', sp_list[1], '.pdf'), width = 4, height = 4)
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

# read in presence absence data
shap_rast_pa_i <- rast(paste0(shap_dirs, "/shap_raster_pa.TIF"))

# clean names
names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']

# remove non_na
shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]

# get mean values to order by
shap_rast_pa_i <- shap_rast_pa_i[[order(var_imp$var_imp, decreasing = T)]]

# get raster of species presences and absences
suit_i <- rast(sp_raster_suit_pa)
presence <- rast(sp_raster_pres_pa[1])
suit_i <- resample(suit_i, presence)

pres_suit <- suit_i
pres_suit[is.na(presence)] <- NA
names(pres_suit) <- 'presence'

abs_suit <- suit_i
abs_suit[!is.na(presence)] <- NA
names(abs_suit) <- 'absence'

pa_map <- c(pres_suit, abs_suit)
pa_map <- trim(pa_map)

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
           col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
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
  
  

# plot overall 
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
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = F, 
            legend.outside.size = 0.1,
            legend.outside.position = c('bottom'), 
            panel.labels=names(shap_rast_plot), 
            panel.label.bg.color = 'white', 
            frame = FALSE,
            bg.color = "transparent", 
            panel.label.size = 1.2)


dir.create(paste0(fig_dir, '/spatial_shap_example/'), recursive = T)
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


# make pdf
pdf(paste0(fig_dir, '/spatial_shap_example/shap_', sp_list[1], '.pdf'), width = 8, height = 6)
print(shap_plot) 
dev.off()

# make png to save space
png(paste0(fig_dir, '/spatial_shap_example/shap_', sp_list[1], '.png'), width = 3500, height = 1800, res = 300)
print(shap_plot) 
dev.off()


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
        #axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(), 
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

# plot correlation matrix
pdf(paste0(fig_dir, 'spatial_shap_example/', 'shap_correlations.pdf'), width = 6, height = 6)
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

### ESTIMATE SHAPLEY AREAS BASED ON CATCHMENTS 
## get the catchment level shapley values
# read in the RDS
sp_shap <- readRDS(shap_pa)

# join shapley values to the teilenzugsgebeit
sp_shap <- left_join(
  subcatchments_final %>%
    select(TEILEZGNR),
  sp_shap
) %>%
  filter(!is.na(suitability))

sp_shap$natural_niche <- (sp_shap$ecoF_discharge_max_log10_SHAP>0) + (sp_shap$stars_t_mn_m_c_SHAP>0) + (sp_shap$ecoF_flow_velocity_mean_SHAP>0) == 3
sp_shap$natural_niche_value <- (sp_shap$ecoF_discharge_max_log10_SHAP + sp_shap$stars_t_mn_m_c_SHAP + sp_shap$ecoF_flow_velocity_mean_SHAP) / 3

round(table(sp_shap$natural_niche) / nrow(sp_shap), 2)
round(table(sp_shap$ecoF_discharge_max_log10_SHAP>0) / nrow(sp_shap), 2)
round(table(sp_shap$stars_t_mn_m_c_SHAP>0) / nrow(sp_shap), 2)
round(table(sp_shap$ecoF_flow_velocity_mean_SHAP>0) / nrow(sp_shap), 2)

## find the sites that fit our 5 categories of outside, inside, inside + negative con, inside + negative habitat, inside + negative both

sp_shap$neg_con <- sp_shap$local_asym_cl_log10_SHAP < 0
sp_shap$neg_habitat <- ((sp_shap$local_flood_SHAP + 
  sp_shap$local_imd_log10_ele_residual_SHAP + 
  sp_shap$ecoF_eco_mean_ele_residual_SHAP) /3) < 0

sp_shap$niche_categories <- as.factor(ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == F & sp_shap$neg_habitat == F, '2. inside ecological niche',
                                      ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == T & sp_shap$neg_habitat == F, '3. poor connectivity', 
                                             ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == F & sp_shap$neg_habitat == T, '4. poor habitat', 
                                                    ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == T & sp_shap$neg_habitat == T, '5. poor connectivity and habitat', 
                                                           ifelse(sp_shap$natural_niche == F, '1. outside ecological niche', NA))))))

# create plot of shapley values against suitability classified by threat
pdf(paste0(fig_dir, '/spatial_shap_example/natural_threat_suitability_biplot.pdf'), 
    width = 3, height = 3, bg = 'transparent')
ggplot() + 
  geom_point(data = sp_shap %>% arrange(niche_categories), 
             aes(y = suitability, 
                 x = natural_niche_value, 
                 col = niche_categories, 
                 size = niche_categories)) + 
  xlab(expression('abiotic niche' ~phi~ 'values')) + 
  ylab('predicted relative occurrence') + 
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
  scale_size_manual(values = c(1,2,2,2,2))
dev.off()
  

## what % of catchments inside the niche are negative connectivity (local_asym_cl_log10_SHAP)
neg_con <- sp_shap %>% filter(natural_niche == T) %>% .$local_asym_cl_log10_SHAP %>%  .[] < 0
round((table(neg_con) / length(neg_con))*100)

## what % of catchments inside the niche are negative ecomorphology (ecoF_eco_mean_ele_residual_SHAP)
neg_ecoF <- sp_shap %>% filter(natural_niche == T) %>% .$ecoF_eco_mean_ele_residual_SHAP  %>%  .[] < 0
round((table(neg_ecoF) / length(neg_ecoF))*100)

## what % of catchments inside the niche are negative ecomorphology (local_imd_log10_ele_residual_SHAP)
neg_urban <- sp_shap %>% filter(natural_niche == T) %>% .$local_imd_log10_ele_residual_SHAP  %>%  .[] < 0
round((table(neg_urban) / length(neg_urban))*100)

## what % of catchments inside the niche are negative ecomorphology (local_flood_SHAP)
neg_flood <- sp_shap %>% filter(natural_niche == T) %>% .$local_flood_SHAP  %>%  .[] < 0
round((table(neg_flood) / length(neg_flood))*100)

# estiate % of catchments jointly affected by habitat threats
round((table((neg_ecoF + neg_urban + neg_flood) >= 1) / length(neg_ecoF))*100)
round((table((neg_ecoF + neg_urban + neg_flood) >= 2) / length(neg_ecoF))*100)
round((table((neg_ecoF + neg_urban + neg_flood) >= 3) / length(neg_ecoF))*100)

round((table((((neg_ecoF + neg_urban + neg_flood) >= 1) + neg_con) == 2) / length(neg_ecoF))*100)



#### find areas where suitability estimates are higher than predicted occurrence threshold
# get mean of raster layers
rast_suit_pa <- mean(rast(sp_raster_suit_pa), na.rm = T)
model_data   <- readRDS(paste0(sdm_dirs, '/data/sdm_input_data.rds'))@pa_data$full_data
pred_occ <- terra::extract(rast_suit_pa, model_data[c('X', 'Y')])
threshold <- ecospat::ecospat.max.tss(pred_occ$mean, model_data$occ)$max.threshold

# define if the catchment is occupied based on mean model predictions
sp_shap$present <- sp_shap$suitability >= threshold
table(sp_shap$present)

# get presence or absence inside niche
niche_pa <- sp_shap %>% filter(natural_niche == T) %>% .$present
round((table(niche_pa) / length(niche_pa))*100)

#### 9. Plot areas that are inside the ecological niche but are threatened ----

### PLOT SHAPLEY AREAS BASED ON RASTERS (MORE CONVININIENT PLOTTING)

# find areas where shapley values for natural niche factors (discharge, flow velocity and temperature) 
# are positive for Alburnoides bipunctatus
natural <- sum(c(shap_rast_pa_i[[grep('discharge', names(shap_rast_pa_i))]] > 0, 
                 shap_rast_pa_i[[grep('flow velocity', names(shap_rast_pa_i))]] > 0, 
                 shap_rast_pa_i[[grep('temperature', names(shap_rast_pa_i))]]) > 0)
all_natural <- natural == 3
plot(natural)
plot(all_natural)

# get an index of 'naturalness' which is the sum of the shapley contributions of all natural factors
natural_index <- sum(c(shap_rast_pa_i[[grep('discharge', names(shap_rast_pa_i))]], 
                       shap_rast_pa_i[[grep('flow velocity', names(shap_rast_pa_i))]], 
                       shap_rast_pa_i[[grep('temperature', names(shap_rast_pa_i))]]))
natural_index[natural_index < 0] <- 0

# which areas are shapley vales negative for connectivity or floodplain area
threatened <- sum(c(shap_rast_pa_i[[grep('connectivity', names(shap_rast_pa_i))]] < 0, 
                    app(shap_rast_pa_i[[grep('floodplains|ecomorphology|urban', names(shap_rast_pa_i))]], mean) < 0))
 
# find areas where shapley values are negative for 'threats' but shapley values are positive for natural niche factors
neg_con <- shap_rast_pa_i[[grep('connectivity', names(shap_rast_pa_i))]] < 0 & natural == 3
neg_habitat<- app(shap_rast_pa_i[[grep('floodplains|ecomorphology|urban', names(shap_rast_pa_i))]], mean) < 0 & natural == 3
plot(threatened != 0 & natural == 3)
plot(threatened == 0 & natural == 3)
plot(neg_con)
plot(neg_habitat)

# assign categories of range constraints based on values of the above rasters
natural_threat_rast <- shap_rast_pa_i[[1]]
values(natural_threat_rast) <- ifelse(natural[] == 3 & neg_con[] == F & neg_habitat[] == F, '2. inside ecological niche',
                                      ifelse(natural[] == 3 & neg_con[] == T & neg_habitat[] == F, '3. poor connectivity', 
                                             ifelse(natural[] == 3 & neg_con[] == F & neg_habitat[] == T, '4. poor habitat', 
                                                    ifelse(natural[] == 3 & neg_con[] == T & neg_habitat[] == T, '5. poor connectivity and habitat', 
                                                           ifelse(natural[] != 3, '1. outside ecological niche', NA)))))
names(natural_threat_rast) <- 'threat map'

## 1 - outside of ecological niche
## 2 - inside of ecological niche and few threats
## 3 - natural areas with connectivity negative
## 4 - natural areas with flood negative
## 5 - natural areas with both connectivity and floodplains negative
threat_niche <- tm_shape(natural_threat_rast, raster.downsample = F) +
  tm_raster(style = 'fixed',
            breaks = 1:4,
            palette = c('#E8E8E8','#06C4C4','#D1D111', '#D16111', '#F70000'), 
            title = '') +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F) + 
  tm_layout(frame = F,
            bg.color = "transparent")


# map the threat maps
pdf(paste0(fig_dir, '/spatial_shap_example/natural_threat_range.pdf'))
threat_niche
dev.off()

# get area of classes
natural_threat_sf <- terra::as.polygons(natural_threat_rast) %>% st_as_sf() 
natural_threat_sf$area <- st_area(natural_threat_sf) 

# get area of natural range potential
natural_area <- terra::as.polygons(all_natural) %>% st_as_sf() %>% st_area %>% .[2]

# percentage of natural area that is under different threat classes
natural_threat_sf$area_natural <- round(natural_threat_sf$area / natural_area, 2)

#### 10. Look at natural factor constrains inside and outside of species predicted areas ----

# mainly discharge = discharge is highest positive
# mainly flow velocity = flow velocity is highest negative
# mainly temperature = temperature is highest positive

# constrained discharge = discharge is highest negative
# constrained flow velocity = flow velocity is highest negative
# constrained temperature = temperature is highest negative

# find areas where natural niche effects are positive
natural_pos <- sum(shap_rast_pa_i[[grep('discharge', names(shap_rast_pa_i))]] > 0, 
                  shap_rast_pa_i[[grep('flow velocity', names(shap_rast_pa_i))]] > 0, 
                  shap_rast_pa_i[[grep('temperature', names(shap_rast_pa_i))]] > 0)

# find areas where natural niche effects are negative
natural_neg<- sum(shap_rast_pa_i[[grep('discharge', names(shap_rast_pa_i))]] < 0, 
                   shap_rast_pa_i[[grep('flow velocity', names(shap_rast_pa_i))]] < 0, 
                   shap_rast_pa_i[[grep('temperature', names(shap_rast_pa_i))]] < 0)

# estimate which niche effect is the MOST positive
which_pos <- shap_rast_pa_i[[grep('discharge|flow velocity|temperature', names(shap_rast_pa_i))]]
which_pos[natural_pos==0] <- NA
which_pos <- terra::app(which_pos, which.max)
values(which_pos) <- ifelse(which_pos[] == 1, 'discharge',
                            ifelse(which_pos[] == 2, 'temperature', 
                                   ifelse(which_pos[]==3, 'flow velocity', NA)))
which_pos[natural_pos!=3] <- NA

# estimate which niche effect is MOST negative
which_neg <- shap_rast_pa_i[[grep('discharge|flow velocity|temperature', names(shap_rast_pa_i))]]
which_neg[natural_neg==0] <- NA
which_neg <- terra::app(which_neg, which.min)
values(which_neg) <- ifelse(which_neg[] == 1, 'discharge',
                            ifelse(which_neg[] == 2, 'temperature', 
                                   ifelse(which_neg[]==3, 'flow velocity', NA)))
which_neg[natural_neg==0] <- NA

# plot positive niche effects
natural_niche_pos <- tm_shape(which_pos) +
  tm_raster(style = 'fixed',
            palette = c('#6a64ed', '#64edd2', '#e6e61c'), 
            title = 'Main positive niche effect') +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'black', lwd = 0.5) + 
  tm_shape(lakes) +
  tm_polygons(border.col = "black", col = "white", legend.show = F, lwd = 0.5) + 
  tm_layout(frame = F,
            bg.color = "transparent")

# plot negative niche effects
natural_niche_neg <- tm_shape(which_neg) +
  tm_raster(style = 'fixed',
            palette = c('#4d0266', '#0297c4', '#aeb800'), 
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
