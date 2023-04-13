### Aim of script is to generate maps of SHAPLEY values for focal species ----
#### 1. Load in packages data ----

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

#### 2. Set directories for loading and saving objects ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/figures-march2023-v5"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_MARCH_v5"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
# sp_list <- list.files(paste0('D:/sdm-pipeline/sdm-run/', RUN))
# sp_list <- unique(records_table$species_name)
sp_list <- sort(c('Thymallus thymallus', 'Squalius cephalus', 
             'Barbus barbus', 'Gobio gobio', 
             'Alburnoides bipunctatus', 'Perca fluviatilis', 
             'Lampetra planeri'))
 

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
sp_raster_pres_po <- paste0(sdm_dirs, "/output/raster_maps/presence_MCC_stack_po_rf.TIF")
sp_raster_pres_pa <- paste0(sdm_dirs, "/output/raster_maps/presence_stack_pa_rf.TIF")

# read in environmental data
env_data <- rast(paste0(dd_env, '/ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif'))



#### 3. Load all spatial objects ----

# read in elevation raster
base_rast <- rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif"))
base_rast <- base_rast[[1]]

# all env rast
all_env_rast <- trim(rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif")))

# set crs
target_crs <- "epsg:3035"

# read in lakes
lakes <- paste0(dd_env, "ch-lake-reference.shp")
lakes <- st_union(st_read(lakes))

# read in rivers
rivers <- paste0(dd_env, "CH_HYDRO.shp")
rivers <- st_read(rivers)
rivers <- rivers %>%
  filter(CUM_LEN > 500, STRAHLE > 1) %>%
  mutate(CUM_LEN_LOG = log(.$CUM_LEN))

# read in catchment object
subcatchment_file <- paste0(dd_ch, "swiss-2km-subcatchments/EZG_Gewaesser.gdb")
catchments_rhine <- c("Aare", "Reuss", "Limmat", "Rhein")

# read in subcatchments, transform, union
subcatchments_rhine_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(FLUSSGB %in% catchments_rhine) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

# mask over rhein
base_rast <- terra::trim(terra::mask(base_rast, vect(subcatchments_rhine_union)))
lakes <- st_intersection(lakes, subcatchments_rhine_union)
rivers <- st_intersection(rivers, subcatchments_rhine_union)
rivers_int <- st_union(rivers %>% dplyr::select(geometry))
river_intersect_lakes <- st_difference(rivers_int, lakes)

# make non-unioned object
subcatchments_rhine <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(FLUSSGB %in% catchments_rhine) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # remove z and m properties that can cause errors later
  st_zm()

# create croping area to inside of switzerland
cropping_sf <-
  read_sf(
    dsn = paste0(dd_ch, "hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95"),
    layer = "basis04"
  ) %>%
  st_union() %>%
  st_transform(., crs = target_crs)

# get intersection between rhine subcatchments and switzerland borders
subcatchments_rhine_v2 <- st_intersection(subcatchments_rhine, cropping_sf)
subcatchments_rhine    <- subcatchments_rhine_v2

# remove lakes from objects
subcatchments_final <- st_difference(subcatchments_rhine, lakes)
river_intersect_lakes <- st_intersection(river_intersect_lakes, cropping_sf)

## extract environmental data
# extract polygon values over environmental data
all_env_subcatchments <- terra::extract(all_env_rast, 
                                        terra::vect(subcatchments_final), 
                                        fun = function(x) mean(x, na.rm = T))

all_env_subcatchments$TEILEZGNR <- subcatchments_final$TEILEZGNR
all_env_subcatchments$ID <- NULL

#### 4. Set up variable names ----

# set vector of focal variables for assigning new names
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
vars_renamed = c('discharge', 'slope', 'flow velocity', 
                 'temperature max', 'temperature min', 
                 'connectivity', 'distance to lake', 
                 'ecomorphology', 'urbanisation',
                 'wetland', 'floodplains')

vars_renamed <- cbind(vars_shap, vars_renamed)


#### 5. Convert shapley polygons to rasters for easier plotting ----

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


#### 6. Plot all variables at once ----

# create output directory
shap_mean_dir <- paste0(fig_dir,'/',  unique(sp_shap_i$species_name), '/')
dir.create(shap_mean_dir, recursive = T)

for (i in 1:length(sp_list)){
  
  do_po <- file.exists(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
  do_pa <- file.exists(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
  
  if(do_po){
  # read in shapley raster
  shap_rast_po_i <- rast(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
  
  # clean names
  names(shap_rast_po_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_po_i)), 'vars_renamed']
  
  # remove non_na
  shap_rast_po_i <- shap_rast_po_i[[which(global(shap_rast_po_i, fun="notNA")!=0)]]
  
  # plot overall 
  shap_all_plot_po <- tmap_grob(tm_shape(shap_rast_po_i) +
          tm_raster(
            style = "cont",
            palette = "RdBu",
            title = "",
            legend.reverse = TRUE,
            legend.is.portrait = T
          ) +
          tm_shape(lakes) +
          tm_polygons(border.col = "gray50", col = "gray95", legend.show = F, lwd = 0.01) + 
          tm_layout(legend.outside = T, 
                    legend.position = c('left', 'top'), 
                    panel.labels=names(shap_rast_po_i)) + 
    tm_facets(ncol = 2))
  
  }else{shap_all_plot_po <- NULL}
  
  if(do_pa){
    
  # read in presence absence data
  shap_rast_pa_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
  
  # clean names
  names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']
  
  # remove non_na
  shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]
  
  # plot overall 
  shap_all_plot_pa <- tmap_grob(tm_shape(shap_rast_pa_i) +
    tm_raster(
      style = "cont",
      palette = "RdBu",
      legend.reverse = TRUE,
      title = '',
      legend.is.portrait = T
    ) +
    tm_shape(lakes) +
    tm_polygons(border.col = "gray50", col = "gray95", legend.show = F, lwd = 0.01) + 
    tm_layout(legend.outside = T, 
              legend.position = c('left', 'top'), 
              panel.labels=names(shap_rast_pa_i)) + 
    tm_facets(ncol = 2))
  
  
  }else{shap_all_plot_pa <- NULL}
  
  shap_fig_dir <- paste0(fig_dir, "/", sp_list[i], "/")
  dir.create(shap_fig_dir, recursive = T)
  pdf(paste0(shap_fig_dir, "shap_all", ".pdf"), width = 14, height = 10)
  print(cowplot::plot_grid(shap_all_plot_po, shap_all_plot_pa))
  dev.off()
  
}


#### 7. Create single pdf of all shapely across each species ----

fig_dir_allsp <- paste0(fig_dir, '-allsp/')
dir.create(fig_dir_allsp, recursive = T)

pdf(paste0(fig_dir_allsp, "shap_maps_all_po", ".pdf"), width = 14, height = 10)
for (i in 1:length(sp_list)){
  
    do_po <- file.exists(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
    
    if(do_po){
      # read in shapley raster
      shap_rast_po_i <- rast(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
      
      # clean names
      names(shap_rast_po_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_po_i)), 'vars_renamed']
      
      # remove non_na
      shap_rast_po_i <- shap_rast_po_i[[which(global(shap_rast_po_i, fun="notNA")!=0)]]
      
      # plot overall 
      shap_all_plot_po <- tm_shape(shap_rast_po_i) +
                                      tm_raster(
                                        style = "cont",
                                        palette = "RdBu",
                                        title = sp_list[i],
                                        legend.reverse = TRUE,
                                        legend.is.portrait = T
                                      ) +
                                      tm_shape(lakes) +
                                      tm_polygons(border.col = "gray50", col = "gray95", legend.show = F, lwd = 0.01) + 
                                      tm_layout(legend.outside = T, 
                                                legend.position = c('left', 'top'), 
                                                panel.labels=names(shap_rast_po_i)) + 
                                      tm_facets(ncol = 3)
      
    }else{shap_all_plot_po <- NULL}
    
    print(shap_all_plot_po)
    
}
dev.off() 
   
pdf(paste0(fig_dir_allsp, "shap_maps_all_pa", ".pdf"), width = 14, height = 10)
for(i in 1:length(sp_list)){
  
  do_pa <- file.exists(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
  
    if(do_pa){
      
      # read in presence absence data
      shap_rast_pa_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
      
      # clean names
      names(shap_rast_pa_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_pa_i)), 'vars_renamed']
      
      # remove non_na
      shap_rast_pa_i <- shap_rast_pa_i[[which(global(shap_rast_pa_i, fun="notNA")!=0)]]
      
      # order by name
      shap_rast_pa_i <- shap_rast_pa_i[[sort(names(shap_rast_pa_i))]]
      
      # plot overall 
      shap_all_plot_pa <- tm_shape(shap_rast_pa_i) +
                                      tm_raster(
                                        style = "cont",
                                        palette = "RdBu",
                                        legend.reverse = TRUE,
                                        title = sp_list[i],
                                        legend.is.portrait = T
                                      ) +
                                      tm_shape(lakes) +
                                      tm_polygons(border.col = "gray50", col = "gray95", legend.show = F, lwd = 0.01) + 
                                      tm_layout(legend.outside = T, 
                                                legend.position = c('left', 'top'), 
                                                panel.labels=names(shap_rast_pa_i)) + 
                                      tm_facets(ncol = 3)
      
      
    }else{shap_all_plot_pa <- NULL}
    
    print(shap_all_plot_pa)

}
dev.off()


#### 8. Create shapley based response curves for interpretations ----

pdf(paste0(fig_dir_allsp, "shap_rc_all_po", ".pdf"), width = 14, height = 14)
for (i in 1:length(sp_list)) {
  
  if (all(sapply(c(paste0(shap_dirs[i], "/shap_raster_po.TIF"), rc_po[i], sp_raster_pres_po[i]), file.exists))) {
    # read in shapley raster
    shap_rast_i <- rast(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
    # read in raw shapleys
    shap_i <- readRDS(shap_po[i])
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
      data.frame
    
    # get the minimum and maximum
    min_y <- min(grouped_rc_shap$shapley, na.rm = T)
    max_y <- max(grouped_rc_shap$shapley, na.rm = T)
    
    print(head(grouped_rc_shap))
    
    # make shapley plot for all variables for a single species
    rc_shap_po_plot <- ggplot(data = grouped_rc_shap %>% group_by(name) %>% sample_n(1000), aes(x=env, y = shapley)) + 
      facet_wrap(~vars_renamed, scales = 'free') +
      geom_point(col = 'gray50') +
      stat_smooth(col = 'black', se = F) +
      geom_hline(aes(yintercept = 0)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            aspect.ratio = 0.5) +
      scale_y_continuous(limits = c(min_y, max_y)) +
      ggtitle(sp_list[i])
    
    print(rc_shap_po_plot)
  
  }
  print(i)
}
dev.off()



# make plots of response curves for presence-absence data
pdf(paste0(fig_dir_allsp, "shap_rc_all_pa", ".pdf"), width = 14, height = 14)
for (i in 1:length(sp_list)) {
  
  if (all(sapply(c(paste0(shap_dirs[i], "/shap_raster_pa.TIF"), rc_pa[i], sp_raster_pres_pa[i]), file.exists))) {
    # read in shapley raster
    shap_rast_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
    # read in raw shapleys
    shap_i <- readRDS(shap_pa[i])
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
      data.frame
    
    # get the minimum and maximum
    min_y <- min(grouped_rc_shap$shapley, na.rm = T)
    max_y <- max(grouped_rc_shap$shapley, na.rm = T)
    
    print(head(grouped_rc_shap))
    
    # make shapley plot for all variables for a single species
    rc_shap_po_plot <- ggplot(data = grouped_rc_shap %>% group_by(name) %>% sample_n(1000), aes(x=env, y = shapley)) + 
      facet_wrap(~vars_renamed, scales = 'free') +
      geom_point(col = 'gray50') +
      stat_smooth(col = 'black', se = F) +
      geom_hline(aes(yintercept = 0)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            aspect.ratio = 0.5) +
      scale_y_continuous(limits = c(min_y, max_y)) +
      ggtitle(sp_list[i])
    
    print(rc_shap_po_plot)
    
  }
  print(i)
}
dev.off()



#### 9. Make plot comparing all species shapley values together ----

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

# bind all together
all_shap_env_pa <- bind_rows(all_shap_env) %>% 
  mutate(vars_renamed = replace(vars_renamed, vars_renamed %in% c("temperature min", "temperature max"), "temperature"))

# estimate variable importance (mean(abs(shapley))) for species and variable
all_shap_env_pa_simple <- all_shap_env_pa %>% 
  group_by(vars_renamed, species_name) %>% 
  do(mean_var_imp = mean(abs(.$shapley), na.rm=T), 
     range_var_imp = sd(.$shapley), na.rm = T) %>% 
  unnest(c(mean_var_imp, range_var_imp))

focal_variables <- c('connectivity', 
                     'cropland', 
                     'floodplains', 
                     'temperature',
                     'urbanisation', 
                     'insecticide', 
                     'livestock')

focal_species <- c('Thymallus thymallus',
                   'Telestes souffia', 
                   'Barbus barbus', 
                   'Lampetra planeri', 
                   'Alburnoides bipunctatus', 
                   'Squalius cephalus', 
                   'Scardinius erythropthalmus', 
                   'Anguilla anguilla', 
                   'Esox lucius')

pdf(paste0(fig_dir_allsp, "shap_select_varimp", ".pdf"), width = 7, height = 7)
ggplot(data = all_shap_env_pa_simple %>% 
         filter(vars_renamed %in% focal_variables, 
                species_name %in% focal_species) %>% 
         mutate(species_name = factor(species_name, levels = rev(sort(unique(species_name)))))) +
  geom_raster(aes(x = vars_renamed, y = species_name, fill = sqrt(mean_var_imp))) +
  scale_fill_viridis_c(name = 'variable importance') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'bottom') + 
  xlab(NULL) + 
  ylab(NULL) 
dev.off()

pdf(paste0(fig_dir_allsp, "shap_select_rc", ".pdf"), width = 13, height = 13)
ggplot(data = all_shap_env_pa %>% 
         filter(vars_renamed %in% focal_variables, 
                species_name %in% focal_species) %>% 
         left_join(., all_shap_env_pa_simple) %>% 
         group_by(vars_renamed, species_name) %>% 
         sample_n(5000)) +
  geom_point(aes(x = env, y = shapley, col = sqrt(mean_var_imp))) +
  #geom_ma(aes(x = env, y = shapley), col = 'black') + 
  facet_grid(species_name ~ vars_renamed, scales = 'free_x') +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        legend.position = 'bottom') +
  scale_y_continuous(limits = c(min_y, max_y)) +
  scale_colour_viridis_c(name = 'variable importance') + 
  xlab('environmental variable') + 
  ylab('effect on habitat suitability \n (shapley value)')
dev.off()


#### 10. Extract certain regions and see focal species responses ----

# get the teilenzugsgebeit for the emme
focal_area <- c(106570, 106299, 101728, 101280, 106548, 103982, 108894)

# read in subcatchments, transform, union
subcatchments_emme_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

subcatchment_2km_names <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>% 
  pull(TEILEZGNR)
subcatchment_2km_names


shapley_emme <- lapply(1:length(sp_list), function(x){ 
  
  if(file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_x <- left_join(shap_x %>% 
                          filter(TEILEZGNR %in% subcatchment_2km_names) %>% 
                          unique(), 
                        all_env_subcatchments %>% 
                          unique(), 
                        by = 'TEILEZGNR')
    
    return(shap_x)
    
    }
  }
)

shapley_emme_focal <- shapley_emme[which(sp_list %in% focal_species)]

shapley_emme_focal <- bind_rows(shapley_emme_focal)

shapley_emme_data <- shapley_emme_focal %>% 
  tibble %>% 
  select(species_name, matches('_SHAP')) %>% 
  pivot_longer(col = -species_name) %>% 
  group_by(species_name, name) %>% 
  do(mean_shap = mean(.$value, na.rm = T)) %>% 
  unnest()

shapley_emme_data$name <- vars_renamed[,2][match(shapley_emme_data$name, data.frame(vars_renamed)[,1])]

shapley_emme_data <- shapley_emme_data %>%   
  mutate(name = replace(name,
                        name %in% c("temperature min", "temperature max"), 
                        "temperature")) %>% 
  filter(name %in% focal_variables)
  

pdf(paste0(fig_dir_allsp, "shap_local_emme", ".pdf"), width = 10, height = 5)
ggplot(data = shapley_emme_data) + 
  geom_bar(aes(y = name, 
               x = mean_shap, 
               fill = name), stat = 'identity', 
           position = 'dodge') + 
  geom_vline(aes(xintercept = 0)) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        aspect.ratio = 1) + 
  facet_wrap(~species_name, nrow = 2) + 
  xlab('Local variable importance (shapley)') + 
  ylab(NULL)
dev.off()


#### 11. Extract limmat river and see focal species responses ----

# get the teilenzugsgebeit for the limmat
focal_area <- c(100134)

# read in subcatchments, transform, union
subcatchments_limmat_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

mapview::mapview(subcatchments_limmat_union)

subcatchment_2km_names <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>% 
  pull(TEILEZGNR)
subcatchment_2km_names


shapley_limmat <- lapply(1:length(sp_list), function(x){ 
  
  if(file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_x <- left_join(shap_x %>% 
                          filter(TEILEZGNR %in% subcatchment_2km_names) %>% 
                          unique(), 
                        all_env_subcatchments %>% 
                          unique(), 
                        by = 'TEILEZGNR')
    
    return(shap_x)
    
  }
}
)

shapley_limmat_focal <- shapley_limmat[which(sp_list %in% focal_species)]

shapley_limmat_focal <- bind_rows(shapley_limmat_focal)

shapley_limmat_data <- shapley_limmat_focal %>% 
  tibble %>% 
  select(species_name, matches('_SHAP')) %>% 
  pivot_longer(col = -species_name) %>% 
  group_by(species_name, name) %>% 
  do(mean_shap = mean(.$value, na.rm = T)) %>% 
  unnest()

shapley_limmat_data$name <- vars_renamed[,2][match(shapley_limmat_data$name, data.frame(vars_renamed)[,1])]

shapley_limmat_data <- shapley_limmat_data %>%   
  mutate(name = replace(name,
                        name %in% c("temperature min", "temperature max"), 
                        "temperature")) %>% 
  filter(name %in% focal_variables)


pdf(paste0(fig_dir_allsp, "shap_local_limmat", ".pdf"), width = 10, height = 5)
ggplot(data = shapley_limmat_data) + 
  geom_bar(aes(y = name, 
               x = mean_shap, 
               fill = name), stat = 'identity', 
           position = 'dodge') + 
  
  geom_vline(aes(xintercept = 0)) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        aspect.ratio = 1) + 
  facet_wrap(~species_name, nrow = 2) + 
  xlab('Local variable importance (shapley)') + 
  ylab(NULL)
dev.off()


#### 12. Extract reuss river and see focal species responses ----

# get the teilenzugsgebeit for the reuss
focal_area <- c(100675)

# read in subcatchments, transform, union
subcatchments_reuss_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

mapview::mapview(subcatchments_reuss_union)

subcatchment_2km_names <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>% 
  pull(TEILEZGNR)
subcatchment_2km_names


shapley_reuss <- lapply(1:length(sp_list), function(x){ 
  
  if(file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_x <- left_join(shap_x %>% 
                          filter(TEILEZGNR %in% subcatchment_2km_names) %>% 
                          unique(), 
                        all_env_subcatchments %>% 
                          unique(), 
                        by = 'TEILEZGNR')
    
    return(shap_x)
    
  }
}
)

shapley_reuss_focal <- shapley_reuss[which(sp_list %in% focal_species)]

shapley_reuss_focal <- bind_rows(shapley_reuss_focal)

shapley_reuss_data <- shapley_reuss_focal %>% 
  tibble %>% 
  select(species_name, matches('_SHAP')) %>% 
  pivot_longer(col = -species_name) %>% 
  group_by(species_name, name) %>% 
  do(mean_shap = mean(.$value, na.rm = T)) %>% 
  unnest()

shapley_reuss_data$name <- vars_renamed[,2][match(shapley_reuss_data$name, data.frame(vars_renamed)[,1])]

shapley_reuss_data <- shapley_reuss_data %>%   
  mutate(name = replace(name,
                        name %in% c("temperature min", "temperature max"), 
                        "temperature")) %>% 
  filter(name %in% focal_variables)


pdf(paste0(fig_dir_allsp, "shap_local_reuss", ".pdf"), width = 10, height = 5)
ggplot(data = shapley_reuss_data) + 
  geom_bar(aes(y = name, 
               x = mean_shap, 
               fill = name), stat = 'identity', 
           position = 'dodge') + 
  
  geom_vline(aes(xintercept = 0)) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        aspect.ratio = 1) + 
  facet_wrap(~species_name, nrow = 2) + 
  xlab('Local variable importance (shapley)') + 
  ylab(NULL)
dev.off()


#### 12. Extract sense river and see focal species responses ----

# get the teilenzugsgebeit for the sense
focal_area <- c(100184)

# read in subcatchments, transform, union
subcatchments_sense_union <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>%
  # convert to target crs
  st_transform(., crs = target_crs) %>%
  # union together
  st_union() %>%
  # remove z and m properties that can cause errors later
  st_zm()

mapview::mapview(subcatchments_sense_union)

subcatchment_2km_names <- st_read(subcatchment_file, layer = "Teileinzugsgebiet") %>%
  filter(TEZGNR40 %in% focal_area) %>% 
  pull(TEILEZGNR)
subcatchment_2km_names


shapley_sense <- lapply(1:length(sp_list), function(x){ 
  
  if(file.exists(shap_pa[x])){
    
    # read in raw shapleys
    shap_x <- readRDS(shap_pa[x])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_x <- left_join(shap_x %>% 
                          filter(TEILEZGNR %in% subcatchment_2km_names) %>% 
                          unique(), 
                        all_env_subcatchments %>% 
                          unique(), 
                        by = 'TEILEZGNR')
    
    return(shap_x)
    
  }
}
)

shapley_sense_focal <- shapley_sense[which(sp_list %in% focal_species)]

shapley_sense_focal <- bind_rows(shapley_sense_focal)

shapley_sense_data <- shapley_sense_focal %>% 
  tibble %>% 
  select(species_name, matches('_SHAP')) %>% 
  pivot_longer(col = -species_name) %>% 
  group_by(species_name, name) %>% 
  do(mean_shap = mean(.$value, na.rm = T)) %>% 
  unnest()

shapley_sense_data$name <- vars_renamed[,2][match(shapley_sense_data$name, data.frame(vars_renamed)[,1])]

shapley_sense_data <- shapley_sense_data %>%   
  mutate(name = replace(name,
                        name %in% c("temperature min", "temperature max"), 
                        "temperature")) %>% 
  filter(name %in% focal_variables)


pdf(paste0(fig_dir_allsp, "shap_local_sense", ".pdf"), width = 10, height = 5)
ggplot(data = shapley_sense_data) + 
  geom_bar(aes(y = name, 
               x = mean_shap, 
               fill = name), stat = 'identity', 
           position = 'dodge') + 
  
  geom_vline(aes(xintercept = 0)) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        aspect.ratio = 1) + 
  facet_wrap(~species_name, nrow = 2) + 
  xlab('Local variable importance (shapley)') + 
  ylab(NULL)
dev.off()



#### 10 . Plot shapley values and response curves in a 4 panel plot (env data, response curve, shapleys x2) ----

for (i in 1:length(sp_list)) {
  
  ### RUN PRESENCE ONLY SHAPLEY ANALYSIS PLOTS
  # check if all the correct files are available
  if (all(sapply(c(paste0(shap_dirs[i], "/shap_raster_po.TIF"), rc_po[i], sp_raster_pres_po[i]), file.exists))) {
    # read in shapley raster
    shap_rast_i <- rast(paste0(shap_dirs[i], "/shap_raster_po.TIF"))

    # read in response curves
    rc_po_i <- readRDS(rc_po[i])
    
    # read in range raster
    pres_po <- mean(rast(sp_raster_pres_po[i]))
    
    pres_po_vect <- st_as_sf(as.polygons(pres_po)) %>% st_union()

    # mask presence and absences
    shap_rast_po_pres <- mask(shap_rast_i, resample(pres_po, shap_rast_i))
    shap_rast_po_abse <- mask(shap_rast_i, resample(pres_po, shap_rast_i), inverse = T)
    
    # read in raw shapleys
    shap_i <- readRDS(shap_po[i])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_i <- left_join(shap_i %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
    
    # create output directory
    shap_plot_po <- list()
    for (j in names(shap_rast_po_pres)) {
      
      # get variable name
      var_name <- as.character(vars_renamed[which(j==vars_renamed[,1]),2])
      
      # plot of environmental data distribution
      env_j <- tmap_grob(tm_shape(all_env_rast[[gsub('_SHAP', '', j)]]) +
                  tm_raster(
                    style = "cont",
                    palette = "RdBu",
                    title = "",
                    legend.reverse = TRUE, 
                    midpoint = median(values(all_env_rast[[gsub('_SHAP', '', j)]]),na.rm = T)
                  ) +
                  tm_layout(title = var_name,
                            title.size = 1) +
                  tm_shape(river_intersect_lakes) +
                  tm_lines(col = "gray50") +
                  tm_shape(lakes) +
                  tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))
      
      
      # shapley when present plot
      shap_plot_pres_j <- tmap_grob(tm_shape(shap_rast_po_pres[[j]]) +
        tm_raster(
          style = "cont",
          palette = "RdBu",
          title = "",
          legend.reverse = TRUE
        ) +
        tm_layout(title = paste0(var_name, ": presence"),
                  title.size = 1) +
        tm_shape(river_intersect_lakes) +
        tm_lines(col = "gray50") +
        tm_shape(lakes) +
        tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))

      # shapley when absent plot
      shap_plot_abse_j <- tmap_grob(tm_shape(shap_rast_po_abse[[j]]) +
        tm_raster(
          style = "cont",
          palette = "RdBu",
          title = "",
          legend.reverse = TRUE
        ) +
        tm_shape(pres_po_vect) +
        tm_polygons(col = "black") +
        tm_layout(title = paste0(var_name, ": absence"),
                  title.size = 1) +
        tm_shape(river_intersect_lakes) +
        tm_lines(col = "gray50") +
        tm_shape(lakes) +
        tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))

      # response curve plot
      rc_po_plot_j <- ggplot(data = rc_po_i %>% filter(variable == gsub("_SHAP", "", j))) +
        geom_line(aes(x = x, y = y_pred), lwd = 1) +
        geom_ribbon(aes(x = x, ymin = y_lwr, ymax = y_upr), alpha = 0.5) +
        geom_hline(aes(yintercept = 0), lty = 2) +
        facet_wrap(~variable, ncol = 5, scales = "free_x") +
        theme_bw() +
        theme(
          aspect.ratio = 0.5,
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none"
        ) +
        ylab("ALE habitat suitability") +
        xlab(var_name)
      
      # other type of response curve
      shap_xy <- shap_i %>% select(matches(j), matches(gsub('_SHAP', '', j)))
      rc_shap_po_plot_j <- ggplot(shap_i, aes(y = !!sym(j), x = !!sym(gsub('_SHAP', '', j)))) + 
        geom_point(col = 'gray50') + 
        stat_smooth(col = 'black', se = F, size = 2) + 
        geom_hline(aes(yintercept = 0)) +
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              aspect.ratio = 0.5) + 
        xlab(var_name) + 
        ylab(paste0(var_name, ' shapley value'))
      

    shap_plot_po[[which(j == names(shap_rast_po_pres))]] <- cowplot::plot_grid(env_j,
                                                                               cowplot::plot_grid(rc_po_plot_j, 
                                                                                                  rc_shap_po_plot_j, 
                                                                                                  nrow = 2),
                                                                               shap_plot_pres_j, 
                                                                               shap_plot_abse_j,
                                                                               ncol = 2)
    }
    
  

    shap_fig_dir <- paste0(fig_dir, "/", sp_list[i], "/")
    dir.create(shap_fig_dir, recursive = T)
    pdf(paste0(shap_fig_dir, "shap_", "PO", ".pdf"), width = 18, height = 15)
    print(shap_plot_po)
    dev.off()
    
    
    # plot variable importance from all shapley values
    shap_i_long <- shap_i %>% 
      select(TEILEZGNR, matches(vars_shap)) %>% 
      pivot_longer(., cols = matches(vars_shap))
    shap_mean <- shap_i_long %>% 
      group_by(name) %>% 
      summarize(mean_value = mean(abs(value), na.rm = T))
    shap_mean <- merge(shap_mean, vars_renamed, by.x = 'name', by.y = 'vars_shap')
    shap_mean$vars_renamed <- factor(shap_mean$vars_renamed, levels = shap_mean$vars_renamed[order(shap_mean$mean_value)])
    shap_var_imp_plot <- ggplot(data = shap_mean) +
      geom_bar(aes(y = vars_renamed, x = mean_value), stat = 'identity') + 
      theme_bw() + 
      xlab('mean(|shap|)') + 
      ylab(NULL)
    
    pdf(paste0(shap_fig_dir, "shap_var_imp_", "PO", ".pdf"))
    print(shap_var_imp_plot)
    dev.off()
    
  
  }
  
  
  ### RUN PRESENCE ABSENCE SHAPLEY ANALYSIS PLOTS
  # check if all the correct files are available
  if (all(sapply(c(paste0(shap_dirs[i], "/shap_raster_pa.TIF"), rc_pa[i], sp_raster_pres_pa[i]), file.exists))) {
    
    # read in shapley raster
    shap_rast_i <- rast(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
    
    # read in response curves
    rc_pa_i <- readRDS(rc_pa[i])
    
    # read in range raster
    pres_pa <- mean(rast(sp_raster_pres_pa[i]))
    
    pres_pa_vect <- st_as_sf(as.polygons(pres_pa)) %>% st_union()

    shap_rast_pa_pres <- mask(shap_rast_i, resample(pres_pa, shap_rast_i))
    shap_rast_pa_abse <- mask(shap_rast_i, resample(pres_pa, shap_rast_i), inverse = T)
    
    # read in raw shapleys
    shap_i <- readRDS(shap_pa[i])
    
    # join together shapley values to environmental data for each TEILEZGNR
    shap_i <- left_join(shap_i %>% unique(), all_env_subcatchments %>% unique(), by = 'TEILEZGNR')
    
    # create output directory
    shap_plot_pa <- list()
    for (j in names(shap_rast_pa_pres)) {
      
      # get variable name
      var_name <- as.character(vars_renamed[which(j==vars_renamed[,1]),2])
      
      # plot of environmental data distribution
      env_j <- tmap_grob(tm_shape(trim(all_env_rast[[gsub('_SHAP', '', j)]])) +
                           tm_raster(
                             style = "cont",
                             palette = "RdBu",
                             title = "",
                             legend.reverse = TRUE, 
                             midpoint = median(values(all_env_rast[[gsub('_SHAP', '', j)]]),na.rm = T)
                           ) +
                           tm_layout(title = var_name,
                                     title.size = 1) +
                           tm_shape(river_intersect_lakes) +
                           tm_lines(col = "gray50") +
                           tm_shape(lakes) +
                           tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))
      
      
      # shapley when present plot
      shap_plot_pres_j <- tmap_grob(tm_shape(shap_rast_pa_pres[[j]]) +
                                      tm_raster(
                                        style = "cont",
                                        palette = "RdBu",
                                        title = "",
                                        legend.reverse = TRUE
                                      ) +
                                      tm_layout(title = paste0(var_name, ": presence"),
                                                title.size = 1) +
                                      tm_shape(river_intersect_lakes) +
                                      tm_lines(col = "gray50") +
                                      tm_shape(lakes) +
                                      tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))
      
      # shapley when absent plot
      shap_plot_abse_j <- tmap_grob(tm_shape(shap_rast_pa_abse[[j]]) +
                                      tm_raster(
                                        style = "cont",
                                        palette = "RdBu",
                                        title = "",
                                        legend.reverse = TRUE
                                      ) +
                                      tm_shape(pres_pa_vect) +
                                      tm_polygons(col = "black") +
                                      tm_layout(title = paste0(var_name, ": absence"),
                                                title.size = 1) +
                                      tm_shape(river_intersect_lakes) +
                                      tm_lines(col = "gray50") +
                                      tm_shape(lakes) +
                                      tm_polygons(border.col = "gray50", col = "gray95", legend.show = F))
      
      # response curve plot
      rc_pa_plot_j <- ggplot(data = rc_pa_i %>% filter(variable == gsub("_SHAP", "", j))) +
        geom_line(aes(x = x, y = y_pred), lwd = 1) +
        geom_ribbon(aes(x = x, ymin = y_lwr, ymax = y_upr), alpha = 0.5) +
        geom_hline(aes(yintercept = 0), lty = 2) +
        facet_wrap(~variable, ncol = 5, scales = "free_x") +
        theme_bw() +
        theme(
          aspect.ratio = 0.5,
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none"
        ) +
        ylab("ALE habitat suitability") +
        xlab(var_name)
      
      
      # other type of response curve
      shap_xy <- shap_i %>% select(matches(j), matches(gsub('_SHAP', '', j)))
      rc_shap_pa_plot_j <- ggplot(shap_i, aes(y = !!sym(j), x = !!sym(gsub('_SHAP', '', j)))) + 
        geom_point(col = 'gray50') + 
        stat_smooth(col = 'black', se = F, size = 2) + 
        geom_hline(aes(yintercept = 0)) +
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              aspect.ratio = 0.5) + 
        xlab(var_name) + 
        ylab(paste0(var_name, ' shapley value'))
      
      
      shap_plot_pa[[which(j == names(shap_rast_pa_pres))]] <- cowplot::plot_grid(env_j,
                                                                                 cowplot::plot_grid(rc_pa_plot_j, 
                                                                                                    rc_shap_pa_plot_j, 
                                                                                                    nrow = 2),
                                                                                 shap_plot_pres_j, 
                                                                                 shap_plot_abse_j,
                                                                                 ncol = 2)
    }
    
    shap_fig_dir <- paste0(fig_dir, "/", sp_list[i], "/")
    dir.create(shap_fig_dir, recursive = T)
    pdf(paste0(shap_fig_dir, "shap_", "PA", ".pdf"), width = 18, height = 15)
    print(shap_plot_pa)
    dev.off()
    
    
    
    # plot variable importance from all shapley values
    shap_i_long <- shap_i %>% 
      select(TEILEZGNR, matches(vars_shap)) %>% 
      pivot_longer(., cols = matches(vars_shap))
    shap_mean <- shap_i_long %>% 
      group_by(name) %>% 
      summarize(mean_value = mean(abs(value), na.rm = T))
    shap_mean <- merge(shap_mean, vars_renamed, by.x = 'name', by.y = 'vars_shap')
    shap_mean$vars_renamed <- factor(shap_mean$vars_renamed, levels = shap_mean$vars_renamed[order(shap_mean$mean_value)])
    shap_var_imp_plot <- ggplot(data = shap_mean) +
      geom_bar(aes(y = vars_renamed, x = mean_value), stat = 'identity') + 
      theme_bw() + 
      xlab('mean(|shap|)') + 
      ylab(NULL)
    
    pdf(paste0(shap_fig_dir, "shap_var_imp_", "PA", ".pdf"))
    print(shap_var_imp_plot)
    dev.off()
    
    
  }
  
}




#### END OF CURRENT SCRIPT? 




# create output directory
shap_mean_dir <- paste0(fig_dir,'/',  unique(sp_shap$species_name), '/')
dir.create(shap_mean_dir, recursive = T)


## plot pdf of mean values per species
pdf(paste0(shap_mean_dir, 'mean_shap_', unique(sp_shap$model_type), '.pdf'), width = 10, height = 10)

# plot overall 
print(tm_shape(shap_rast) +
  tm_raster(
    style = "cont",
    palette = "RdBu",
    title = "",
    legend.reverse = TRUE,
    legend.is.portrait = F
  ) +
  tm_shape(lakes) +
  tm_polygons(border.col = "gray50", col = "gray95", legend.show = F, lwd = 0.01) + 
  tm_layout(legend.outside = T, 
            legend.position = c('left', 'bottom')))

# plot each
lapply(vars_shap, function(x) {
  if(x %in% names(sp_shap_vect)){
    shap_rast_i <- shap_rast[[x]]
    print(tm_shape(shap_rast_i) +
            tm_raster(
              style = "cont",
              palette = "RdBu",
              title = "",
              legend.reverse = TRUE
            ) +
            tm_layout(main.title = x) + 
            tm_shape(river_intersect_lakes) +
            tm_lines(col = "gray50") +
            tm_shape(lakes) +
            tm_polygons(border.col = "gray50", col = "gray95", legend.show = F) + 
            tm_shape(sp_range) + 
            tm_borders(col = "black", lwd = 0.5))
  }
})

dev.off()




#### For each species and variable, map the shapley values and the response curves -----

# get species directories
sdm_dirs <- list.files(paste0('D:/sdm-pipeline/sdm-run/', RUN), full.names = T)

# response curve paths
rc_po <- paste0(sdm_dirs, '/output/response_curves/response_curves_rf_po.rds')
rc_pa <- paste0(sdm_dirs, '/output/response_curves/response_curves_rf_pa.rds')

# subset to focal species
rc_po <- rc_po[grepl(paste0(sp_list, collapse = '|'), rc_po)]
rc_pa <- rc_pa[grepl(paste0(sp_list, collapse = '|'), rc_pa)]

# run response curve
if(file.exists(rc_po[i])){
  
  # read in response curve data
  rc_rf_po <- readRDS(rc_po[i])
  
  rc_list <- list()
  shap_list <- list()
  
  for(j in unique(rc_rf_po$variable)){
  rc_list[[j]] <- ggplot(data = rc_rf_po %>% filter(variable == j)) +
    geom_line(aes(x = x, y = y_pred), lwd = 1) +
    geom_ribbon(aes(x = x, ymin = y_lwr, ymax = y_upr), alpha = 0.5) +
    geom_hline(aes(yintercept = 0), lty = 2) +
    facet_wrap(~variable, ncol = 5, scales = 'free_x') +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, colour = "gray40", size = 7),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    ylab("ALE habitat suitability") +
    xlab("covariate value")
  
  
  
  }
}












