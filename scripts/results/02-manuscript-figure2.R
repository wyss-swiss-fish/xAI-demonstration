### Figure 2. Explore all Shapley values for a given species

### Alburnoides bipunctatus maps of spatial distribution of key variables ----


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
sp_list <- c("Alburnoides bipunctatus")
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



#### Maps of shapley values ordered by mean mod of shapleys ----

# load spatial objects
source('scripts/results/00-load-spatial.R')

# set vector of focal variables for assigning new names
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


#### save rasters of shapley values ----

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


#### Map spatial distribution of shapley values for Alburnoides bipunctatus ---- 

# read in presence absence data
shap_rast_pa_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))

# clean names
names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']

# remove non_na
shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]

# get mean values to order by
abs_shap <- sapply(abs(shap_rast_pa_i), function(x) mean(values(x), na.rm = T))
abs_shap_filter <- which(sort(abs_shap, decreasing = T) > 0.01)
shap_rast_pa_i <- shap_rast_pa_i[[order(abs_shap, decreasing = T)[abs_shap_filter]]]

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
  tm_shape(pa_map) + 
  tm_raster(style = 'cont', 
            title = 'Occurrence probability', 
            legend.is.portrait = F, 
            legend.show = F) + 
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, 
           col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = T,
            legend.outside.size = 0.1,
            legend.outside.position = c('bottom'), 
            panel.labels=c('presence (p>0.5)', 'absence (p<0.5)'), 
            panel.label.bg.color = 'white', 
            frame = FALSE, 
            bg.color = "transparent") + 
  tm_facets(ncol = 1)
  
  

# plot overall 
shap_plot <- tm_shape(shap_rast_pa_i) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = TRUE,
    title = 'Shapley values',
    legend.is.portrait = F, 
    legend.show = F
  ) +
  tm_shape(river_intersect_lakes) + 
  tm_lines(legend.show = F, col = 'gray75') + 
  tm_shape(lakes) +
  tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = T, 
            legend.outside.size = 0.1,
            legend.outside.position = c('bottom'), 
            panel.labels=names(shap_rast_pa_i), 
            panel.label.bg.color = 'white', 
            frame = FALSE,
            bg.color = "transparent", 
            panel.label.size = 1.2) + 
  tm_facets(ncol = 3)


dir.create(paste0(fig_dir, '/spatial_shap_example/'), recursive = T)
pdf(paste0(fig_dir, '/spatial_shap_example/presence_', sp_list[1], '.pdf'), width = 6, height = 4, 
    bg = 'transparent')
print(occ_plot) 
dev.off()

pdf(paste0(fig_dir, '/spatial_shap_example/presence_legend', sp_list[1], '.pdf'), width = 2, height = 2, 
    bg = 'transparent')
print(tm_shape(pa_map[[1]]) + 
        tm_raster(style = 'cont', 
                  title = 'Occurrence probability', 
                  legend.is.portrait = F, 
                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                  labels = as.character(c(0, 0.25, 0.5, 0.75, 1))) +
        tm_layout(legend.only = T)) 
dev.off()


pdf(paste0(fig_dir, '/spatial_shap_example/shap_', sp_list[1], '.pdf'), width = 6*3, height = 4)
print(shap_plot) 
dev.off()

pdf(paste0(fig_dir, '/spatial_shap_example/shap_legend_', sp_list[1], '.pdf'), width = 2, height = 2)
print(tm_shape(shap_rast_pa_i[[1]]) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    legend.reverse = FALSE,
    title = 'Shapley values',
    legend.is.portrait = F) + 
  tm_layout(legend.only = T))
dev.off()
# cowplot::plot_grid














