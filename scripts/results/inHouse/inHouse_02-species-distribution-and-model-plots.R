#### Aim of script is to generate series of figures summarise SMDs and shapley values for each species ----

# This script generates a series of figures for exploration including:
# - maps of species shapley values for presence-only and presence-absence SDMs (although we further consider only PA models)
# - response curves of shapley values for presence-only and presence-absence SDMs (although we further consider only PA models)

#### 1. Load in packages data ----

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, tmap, gridExtra, tidyquant)

#### 2. Set directories for loading and saving objects ----

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# figure directory
fig_dir <- "figures/ubelix_SDM_RF_JULY_inHouse_V1/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_JULY_inHouse_V1"

# get species of interest
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010_inHouse.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

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

# suitability map
sp_raster_suit_pa <- paste0(sdm_dirs, "/output/raster_maps/suitability_stack_pa_rf.TIF")

# read in environmental data
env_data <- rast(paste0(dd_env, '/ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif'))

# load in all spatial objects for making maps and plots
source('scripts/results/00-load-spatial.R')

#### 3. Set up variable names ----

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


#### 4. Convert shapley polygons to rasters for easier plotting ----

for (i in 1:length(sp_list)) {
  # check if PO files exist
  if (file.exists(shap_po[i]) & !file.exists(paste0(shap_dirs[i], "/shap_raster_po.TIF"))) {
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
  if (file.exists(shap_pa[i]) & !file.exists(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))) {
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

#### 5. Make single pdf spatial maps of shapely values across each species ----

dir.create(paste0(fig_dir, '/all_sp_plots/'), recursive = T)

pdf(paste0(fig_dir, "all_sp_plots/shap_maps_all_po", ".pdf"), width = 14, height = 10)
for (i in 1:length(sp_list)){
  
  
  do_po <- file.exists(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
  
  if(do_po){
    
    # read in shapley values data for a given species
    shap_rast_po_i <- rast(paste0(shap_dirs[i], "/shap_raster_po.TIF"))
    
    # clean names
    names(shap_rast_po_i) <- vars_renamed[which(vars_renamed[,'vars_shap']%in%names(shap_rast_po_i)), 'vars_renamed']
    
    # remove non_na
    shap_rast_po_i <- shap_rast_po_i[[which(global(shap_rast_po_i, fun="notNA")!=0)]]
    
    # order by name
    shap_rast_po_i <- shap_rast_po_i[[sort(names(shap_rast_po_i))]]
    
    # plot overall 
    shap_all_plot_po <- tm_shape(shap_rast_po_i) +
      tm_raster(
        style = "cont",
        palette = c('#bd0f06','gray90', '#2200c9'), 
        midpoint = 0,
        legend.reverse = F,
        legend.is.portrait = F, 
        legend.show = T,
        title = 'Shapley'
      ) +
      tm_layout(title = sp_list[i]) + 
      tm_facets(ncol = 4) +  
      tm_shape(river_intersect_lakes) + 
      tm_lines(legend.show = F, col = 'gray75') + 
      tm_shape(lakes) +
      tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
      tm_layout(panel.labels=names(shap_rast_po_i), 
                panel.label.bg.color = 'white', 
                frame = FALSE,
                bg.color = "transparent", 
                panel.label.size = 1.2, 
                legend.outside.position = 'bottom')
    
    
  }else{shap_all_plot_po <- NULL}
  
  print(shap_all_plot_po)
    
}
dev.off() 
   
pdf(paste0(fig_dir, "all_sp_plots/shap_maps_all_pa", ".pdf"), width = 24*0.75, height = 10*0.75)
for(i in 1:length(sp_list)){
  
  do_pa <- file.exists(paste0(shap_dirs[i], "/shap_raster_pa.TIF"))
  
    if(do_pa){
      
      # read in shapley values data for a given species
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
          palette = c('#bd0f06','gray90', '#2200c9'), 
          midpoint = 0,
          legend.reverse = F,
          legend.is.portrait = F, 
          legend.show = T,
          title = 'Shapley'
        ) +
        tm_layout(title = sp_list[i]) + 
        tm_facets(ncol = 4) +  
        tm_shape(river_intersect_lakes) + 
        tm_lines(legend.show = F, col = 'gray75') + 
        tm_shape(lakes) +
        tm_polygons(border.col = "gray75", col = "white", legend.show = F, lwd = 0.01) + 
        tm_layout(panel.labels=names(shap_rast_pa_i), 
                  panel.label.bg.color = 'white', 
                  frame = FALSE,
                  bg.color = "transparent", 
                  panel.label.size = 1.2, 
                  legend.outside.position = 'bottom')
      
      
    }else{shap_all_plot_pa <- NULL}
    
    print(shap_all_plot_pa)

}
dev.off()


#### 6. Make shapley based response curves for interpretations ----

pdf(paste0(fig_dir, "all_sp_plots/shap_rc_all_po", ".pdf"), width = 14, height = 14)
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
pdf(paste0(fig_dir, "all_sp_plots/shap_rc_all_pa", ".pdf"), width = 14, height = 14)
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




#### 7. Make plot comparing all species shapley values together ----

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
all_shap_env_pa <- bind_rows(all_shap_env) #%>% 
  #mutate(vars_renamed = replace(vars_renamed, vars_renamed %in% c("temperature min", "temperature max"), "temperature"))

# estimate variable importance (mean(abs(shapley))) for species and variable
all_shap_env_pa_simple <- all_shap_env_pa %>% 
  group_by(vars_renamed, species_name) %>% 
  do(mean_var_imp = mean(abs(.$shapley), na.rm=T), 
     range_var_imp = sd(.$shapley), na.rm = T) %>% 
  unnest(c(mean_var_imp, range_var_imp))

focal_variables <- vars_renamed[,'vars_renamed']

focal_species <- sp_list

# plot variable importance matrix with species X variable
pdf(paste0(fig_dir, "all_sp_plots/all_sp_varimp", ".pdf"), width = 7, height = 7)
ggplot(data = all_shap_env_pa_simple %>% 
         filter(vars_renamed %in% focal_variables, 
                species_name %in% focal_species) %>% 
         mutate(species_name = factor(species_name, levels = rev(sort(unique(species_name)))))) +
  geom_raster(aes(x = vars_renamed, y = species_name, fill = mean_var_imp)) +
  scale_fill_viridis_c(name = 'variable importance') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'bottom', 
        panel.grid = element_blank()) + 
  xlab(NULL) + 
  ylab(NULL) 
dev.off()

# plot the response curves as a large matrix across all species
pdf(paste0(fig_dir, "all_sp_plots/allsp_rc", ".pdf"), width = 13, height = 13)
ggplot(data = all_shap_env_pa %>% 
         filter(vars_renamed %in% focal_variables, 
                species_name %in% focal_species) %>% 
         left_join(., all_shap_env_pa_simple) %>% 
         group_by(vars_renamed, species_name) %>% 
         sample_n(5000)) +
  geom_point(aes(x = env, y = shapley, col = mean_var_imp)) +
  stat_smooth(aes(x = env, y = shapley), col = 'black', se = F) +
  facet_grid(species_name ~ vars_renamed, scales = 'free_x') +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        legend.position = 'bottom') +
  scale_colour_viridis_c(name = 'variable importance') + 
  xlab('environmental variable') + 
  ylab('Shapley value')
dev.off()


#### 8. Make plots of species habitat suitability ---- 

# species raster distributions
sp_rast <- rast(lapply(sp_raster_suit_pa, function(x) rast(x)))
names(sp_rast) <- sp_list

# define species thresholds based on TSS
threshold_poly <-lapply(sp_list, function(x){
  model_data   <- readRDS(paste0(sdm_dirs[grep(x,sdm_dirs)], '/data/sdm_input_data.rds'))@pa_data$full_data
  pred_occ <- terra::extract(sp_rast[x], model_data[c('X', 'Y')])
  threshold <- ecospat::ecospat.max.tss(pred_occ[[x]], model_data$occ)$max.threshold
  new_rast <- sp_rast[x]
  new_rast[new_rast<threshold] <- NA
  new_rast[new_rast>threshold] <- 1
  return(new_rast)})
threshold_poly <- rast(threshold_poly)

# define if the catchment is occupied based on mean model predictions
sp_shap$present <- sp_shap$suitability >= threshold
table(sp_shap$present)

pdf(paste0(fig_dir, 'all_sp_plots/allsp_suitability.pdf'), width = 12, height = 8)
tm_shape(sp_rast) + 
  tm_raster(style = 'cont', 
            title = 'suitability') + 
  tm_facets(ncol = 3) + 
  tm_shape(threshold_poly) + 
  tm_raster(title = 'predicted present', 
            palette = 'red')
dev.off()

