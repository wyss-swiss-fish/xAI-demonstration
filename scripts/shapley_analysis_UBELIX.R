
i <- as.numeric(commandArgs(trailingOnly = TRUE))
print(i)


#### 1. Load in packages data ----

## script to perform shapely analysis across multiple species
pacman::p_load(tidyverse, sf, terra, randomForest, fastshap, parallel)

#### 2. Set directories ----

# data import locations
# dd <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'
# dd_env <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/'
# dd_ch <- 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/'
dd     <- 'data/'
dd_env <- 'data/sdm-pipeline/env-data/' 
dd_ch  <- 'data/ch-spatial-products/' # other swiss spatial products

# set sdm run location
run_name <- paste0('ubelix_test_PARA_RF_v3', '/')
# run_dir <- paste0('D:/sdm-pipeline/sdm-run/', run_name)
run_dir  <- paste0('sdm-run/', run_name)

# define directory for all outputs
shap_dir <- paste0('shapley-run/', run_name)
dir.create(shap_dir, recursive = T)

# get the records list to determine what kind of model was fitted in the same way as it is run on UBELIX
record_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
record_table <- record_table %>% filter(species_name %in% unique(record_table$species_name)[i])

# get species name
sp_name <- unique(record_table$species_name)

# set directory for saving shapley outputs for species
save_dir <- paste0(shap_dir, '/', sp_name, '/')
print(save_dir)

dir.create(save_dir, recursive = T)


#### 3. Load all spatial objects ----

# read in elevation raster
base_rast <- rast(paste0(dd_env, "ch-rasters/final-raster-set/all_env_data_RHEIN_RESIDUALS.tif"))
base_rast <- base_rast[[1]]

# set crs
target_crs <- 'epsg:3035'

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
subcatchment_file <- paste0(dd_ch, 'swiss-2km-subcatchments/EZG_Gewaesser.gdb')
catchments_rhine <- c('Aare', 'Reuss', 'Limmat', 'Rhein')

# read in subcatchments, transform, union
subcatchments_rhine_union <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet') %>% 
  filter(FLUSSGB %in% catchments_rhine) %>% 
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  # union together
  st_union() %>% 
  # remove z and m properties that can cause errors later
  st_zm()

# mask over rhein
base_rast  <- terra::trim(terra::mask(base_rast, vect(subcatchments_rhine_union)))
lakes      <- st_intersection(lakes, subcatchments_rhine_union)
rivers     <- st_intersection(rivers, subcatchments_rhine_union)
rivers_int <- st_union(rivers %>% dplyr::select(geometry))
river_intersect_lakes <- st_difference(rivers_int, lakes)

# make non-unioned object
subcatchments_rhine <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet') %>% 
  filter(FLUSSGB %in% catchments_rhine) %>% 
  # convert to target crs
  st_transform(., crs = target_crs) %>% 
  # remove z and m properties that can cause errors later
  st_zm()

# create croping area to inside of switzerland
cropping_sf <- 
  read_sf(dsn = paste0(dd_ch, 'hydrographische-gliederungderschweiz/Hydrografische+Gliederung/Hydrografische Gliederung_LV95'),
          layer = 'basis04') %>% 
  st_union() %>% 
  st_transform(., crs = target_crs)

# get intersection between rhine subcatchments and switzerland borders
subcatchments_rhine_v2 <- st_intersection(subcatchments_rhine, cropping_sf)
subcatchments_rhine <- subcatchments_rhine_v2

# remove lakes from objects
subcatchments_final   <- st_difference(subcatchments_rhine, lakes)
river_intersect_lakes <- st_intersection(river_intersect_lakes, cropping_sf)



print('completed loading of all spatial data')

#### 4. List locations of sdm files to use as inputs----

# INPUTS
## model = final random forest models fitted in the sdm-pipeline project scripts
## model_final_data = data used to fit the models in the sdm-pipeline project scripts 
## env_data = environmental data used to fit the models in the sdm-pipeline (constrained to species known distributions)

# get species directories
sp_dir <- paste0(run_dir, sp_name)

# get the path to environmental data for all species
sp_env_path <- paste0(sp_dir, '/data/env_data.TIF')

# get path to data used to fit models
sp_data_po <- paste0(sp_dir, '/data/sdm_input_data_rf_poFinal.rds')
sp_data_pa <- paste0(sp_dir, '/data/sdm_input_data.rds')

# get path to species distribution models rasters
sp_raster_suit_po <- paste0(sp_dir, '/output/raster_maps/suitability_stack_po_rf.TIF')
sp_raster_pres_po <- paste0(sp_dir, '/output/raster_maps/presence_BIN_stack_po_rf.TIF')
sp_raster_suit_pa <- paste0(sp_dir, '/output/raster_maps/suitability_stack_pa_rf.TIF')
sp_raster_pres_pa <- paste0(sp_dir, '/output/raster_maps/presence_stack_pa_rf.TIF')

# get paths to secies distribution model object
sp_RF_po <- paste0(sp_dir, '/output/final_models_rf_po.RDS')
sp_RF_pa <- paste0(sp_dir, '/output/final_models_rf_pa.RDS')

# add logic here to define what data is available for this modelling run
fitting <- c("presence_only", "presence_absence")[c(record_table$evaluate_po, record_table$evaluate_pa)]

print('complete loading of species directory locations')


#### MAKE NEW ENVIRONMENTAL DATA OBJECT FOR SPECIES ----

# extract the environmental data values for the modelled variables
env_poly <- terra::extract(
  x = rast(sp_env_path),
  y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script
  fun = function(x) {
    mean(x, na.rm = T)
  }
)

# bind back in environmental data from extractions
subcatchments_final_env <- cbind(subcatchments_final, env_poly)

# drop geometrics and keep only catchments with data in the final objects
subcatchments_final_env_df <- st_drop_geometry(subcatchments_final_env) %>% na.omit()


#### RUN SHAPLEY ANALYSIS FOR PRESENCE-ONLY MODELS IF AVAILABLE ----


# first aggregate suitabilities to catchments to speed up estimation of shapley values
if ("presence_only" %in% fitting) {
 
 if (!file.exists(sp_RF_po)) {
   next()
 }

 #### SUITABILITY ACROSS SUBCATCHMENTS
 # get mean raster
 rast_suit_po <- mean(rast(sp_raster_suit_po), na.rm = T)

 # extract the raster values per subcatchment
 po_suit_poly <- terra::extract(
   x = rast_suit_po,
   y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script
   fun = function(x) {
     mean(x, na.rm = T)
   }
 )

 # bind back in environmental data from extractions
 subcatchments_suit_po <- cbind(subcatchments_final, po_suit_poly)
 subcatchments_suit_po_df <- st_drop_geometry(subcatchments_suit_po) %>% na.omit()

 
 #### PRESENCE ACROSS SUBCATCHMENTS
 rast_pres_po <- mean(rast(sp_raster_pres_po), na.rm = T)

 po_pres_poly <- terra::extract(
   x = rast_pres_po,
   y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script
   fun = function(x) {
     mean(x, na.rm = T)
   }
 )

 # bind back in environmental data from extractions
 subcatchments_pres_po <- cbind(subcatchments_final, po_pres_poly)
 subcatchments_pres_po_df <- st_drop_geometry(subcatchments_pres_po) %>% na.omit()

 subcatchments_sdm_po <- left_join(
   subcatchments_suit_po_df %>% select(TEILEZGNR, mean) %>% rename(., c("suitability" = "mean")),
   subcatchments_pres_po_df %>% select(TEILEZGNR, mean) %>% rename(., c("presence" = "mean"))
 )


 #### SHAPLEY ANALYSIS
 sp_RF_po_i <- readRDS(sp_RF_po)
 sp_data_po_i <- readRDS(sp_data_po)$rf

 # run shapleys over all model iterations
 sp_shapley_po <- mclapply(1:length(sp_RF_po_i), FUN = function(x) {
   
   # get the variables
   vars <- colnames(attr(sp_RF_po_i[[x]]$terms, "factors"))
   
   # define the prediction function to use
   pfun <- function(object, newdata) {
     as.numeric(as.character(predict(object, newdata = newdata, type = 'prob')[,2]))
   }
   
   # get the shapley values
   shapley <- fastshap::explain(
     
     # model object
     object = sp_RF_po_i[[x]], 
     # names of features to explain
     feature_names = vars, 
     # X data used to fit the model
     X = sp_data_po_i[[x]][vars], 
     # new data to predict on
     newdata = subcatchments_final_env_df[vars],
     # predictive function
     pred_wrapper = pfun, 
     # number of replicates
     nsim = 1000
     
   )
   
   # aggregate and clean up
   shapley <- shapley %>% data.frame # %>% mutate_if(., is.numeric, .funs = round, 2)
   
   # rename shapley data
   names(shapley) <- paste0(names(shapley), '_SHAP')
   
   # bind back with the new_data
   shapley <- cbind(subcatchments_final_env_df, shapley)
   
   return(shapley)
   
 }, mc.cores = 5)

 sp_shapley_po <- lapply(sp_shapley_po, function(x) {
   left_join(left_join(x, subcatchments_final_env), subcatchments_sdm_po)
 })

 sp_shapley_po <- bind_rows(sp_shapley_po, .id = "model_run")
 
 sp_shapley_po <- sp_shapley_po %>% dplyr::select(TEILEZGNR, suitability, presence, matches('_SHAP'))
 
 sp_shapley_po$model_type = 'PO'
 sp_shapley_po$species_name = sp_name
 sp_shapley_po$nrep = 1000

 saveRDS(sp_shapley_po, paste0(save_dir, "shapley_rf_po.RDS"))
 
}




#### RUN SHAPLEY ANALYSIS FOR PRESENCE ABSENCE MODELS IF AVAILABLE ----

if ("presence_absence" %in% fitting) {
 if (!file.exists(sp_RF_pa)) {
   next()
 }

 #### SUITABILITY ACROSS SUBCATCHMENTS
 # get mean of raster layers
 rast_suit_pa <- mean(rast(sp_raster_suit_pa), na.rm = T)

 # extract the raster values per subcatchment
 pa_suit_poly <- terra::extract(
   x = rast_suit_pa,
   y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script
   fun = function(x) {
     mean(x, na.rm = T)
   }
 )

 # bind back in environmental data from extractions
 subcatchments_suit_pa <- cbind(subcatchments_final, pa_suit_poly)
 subcatchments_suit_pa_df <- st_drop_geometry(subcatchments_suit_pa) %>% na.omit()

 
 #### PRESENCE ACROSS SUBCATCHMENTS
 # get mean of raster layers
 rast_pres_pa <- mean(rast(sp_raster_pres_pa), na.rm = T)

 pa_pres_poly <- terra::extract(
   x = rast_pres_pa,
   y = vect(subcatchments_final), # catchment data produced in the spatial_data.R script
   fun = function(x) {
     mean(x, na.rm = T)
   }
 )

 # bind back in environmental data from extractions
 subcatchments_pres_pa    <- cbind(subcatchments_final, pa_pres_poly)
 subcatchments_pres_pa_df <- st_drop_geometry(subcatchments_pres_pa) %>% na.omit()

 
 # join together sdms predictions for occupancy probability and presence 
 subcatchments_sdm_pa <- left_join(
   subcatchments_suit_pa_df %>% select(TEILEZGNR, mean) %>% rename(., c("suitability" = "mean")),
   subcatchments_pres_pa_df %>% select(TEILEZGNR, mean) %>% rename(., c("presence" = "mean"))
 )

 
 
 #### PERFORM SHAPLEY FUNCTION FOR GIVEN SET OF DATA
 sp_RF_pa_i <- readRDS(sp_RF_pa)
 sp_data_pa_i <- readRDS(sp_data_pa)@pa_data$full_data

 # run shapleys over all model iterations
   # get the variables
   vars <- colnames(attr(sp_RF_pa_i$terms, "factors"))
   
   # define the prediction function to use
   pfun <- function(object, newdata) {
     as.numeric(as.character(predict(object, newdata = newdata, type = 'prob')[,2]))
   }
   
   # get the shapley values
   shapley_pa <- fastshap::explain(
     
     # model object
     object = sp_RF_pa_i, 
     # names of features to explain
     feature_names = vars, 
     # X data used to fit the model
     X = sp_data_pa_i[vars], 
     # new data to predict on
     newdata = subcatchments_final_env_df[vars],
     # predictive function
     pred_wrapper = pfun, 
     # number of replicates
     nsim = 1000
     
     )
   
   # aggregate and clean up
   shapley_pa <- shapley_pa %>% data.frame %>% mutate_if(., is.numeric, .funs = round, 2)
   
   # rename shapley data
   names(shapley_pa) <- paste0(names(shapley_pa), '_SHAP')
   
   # bind back with the new_data
   sp_shapley_pa <- cbind(subcatchments_final_env_df, shapley_pa)
 
 sp_shapley_pa <- left_join(left_join(sp_shapley_pa, subcatchments_final_env), subcatchments_sdm_pa)
 
 sp_shapley_pa <- sp_shapley_pa %>% dplyr::select(TEILEZGNR, suitability, presence, matches('_SHAP'))
 
 sp_shapley_pa$model_type = 'PA'
 sp_shapley_pa$species_name = sp_name
 sp_shapley_pa$nrep = 1000
 
 saveRDS(sp_shapley_pa, paste0(save_dir, "shapley_rf_pa.RDS"))
 
}



