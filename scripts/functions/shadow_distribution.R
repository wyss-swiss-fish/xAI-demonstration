## function to summarise the shadow distribution

# Input
# Subcatchments file (subcatchments_final)
# Spatial file from shapley values (shap_pa)
# Model data (sdm_input_data.rds)
# Raster of predicted suitabilities (sp_raster_suit_pa)
# Vector of natural niche factors
# Vector of habitat threats
# Vector of connectivity threats

# this is the location of the sdm data input directory
 #sdm_input_data = sdm_dirs[1]

# this is the location of the rasterized suitabilities
 #raster_data <- sp_raster_suit_pa

# this is the location of the shapely value predictions per sub-catchment
 #shap <- shap_pa

# define the natural niche factors
 #natural_niche_factors <- c('ecoF_discharge_max_log10_SHAP', 'stars_t_mn_m_c_SHAP', 'ecoF_flow_velocity_mean_SHAP', 'local_dis2lake_SHAP')
 #habitat_factors <- c('local_wet_SHAP', 'local_flood_SHAP', 'local_imd_log10_ele_residual_SHAP', 'ecoF_eco_mean_ele_residual_SHAP')
 #conn_factors <- 'local_asym_cl_log10_SHAP'

# define species
 #species <- sp_list[1]

# define the output folder
 #output_folder <- paste0('figures/ubelix_SDM_RF_MARCH_v6/shadow_dist_summaries/')


shadow_distribution <- function(sdm_input_data, 
                             raster_data, 
                             shap, 
                             natural_niche_factors, 
                             habitat_factors,
                             conn_factors, 
                             species,
                             output_folder){
  
  

  # if any of the file locations do not exist then stop
  files <- c(sdm_input_data, raster_data, shap)
  stopifnot(sum(sapply(files, file.exists))==3)
  
  # if any of the files do not contain the species name then stop
  stopifnot(sum(sapply(files, function(x) grepl(species, x)))==3)

  #### READ IN ALL SETS OF DATA NEEDED TO PRODUCE SUMMARY PROPERTIES
  
  # get the suitability threshold by optimising for TSS using predicted suitabilities and raw occurrences
  # data used to fit the models
  model_data <- readRDS(paste0(sdm_input_data, '/data/sdm_input_data.rds'))@pa_data$full_data
  # spatial prediction of habitat suitability
  sp_rast <- rast(raster_data)
  # extract suitability per presence-absence observation
  pred_occ <- terra::extract(sp_rast, model_data[c('X', 'Y')])
  # define threshold based on optimising prediction of presence-absence
  threshold <- ecospat::ecospat.max.tss(pred_occ[,2], model_data$occ)$max.threshold
  
  # read in and join together the subcatchment data and the spatial shapley data
  # read in shapley values
  sp_shap <- readRDS(shap)
  sp_shap <- left_join(
    # join to subcatchment information read in from the spatial data processing scripts
    subcatchments_final %>%
      select(TEILEZGNR),
    sp_shap
  ) %>%
    # remove areas where suitability is not predicted (no covariate data)
    filter(!is.na(suitability)) %>% 
    # apply threshold suitabilities to obtain predicted presence-absence
    mutate(presence = ifelse(suitability < threshold, NA, suitability), 
           threshold = threshold)
  
  
  
  #### GENERATE SUMMARY PROPERTIES TO REPORT: ECOLOGICAL NICHE
  natural_niche_factors <- natural_niche_factors[natural_niche_factors %in% names(sp_shap)]
  habitat_factors <- habitat_factors[habitat_factors %in% names(sp_shap)]
  conn_factors <- conn_factors[conn_factors %in% names(sp_shap)]
  
  all_threats <- c(habitat_factors, conn_factors)
  all_threats <- all_threats[all_threats %in% names(sp_shap)]
  
  # define whether a subcatchment is inside or outside the environmental niche
  sp_shap$natural_niche <- rowSums(st_drop_geometry(sp_shap[,natural_niche_factors])) > 0
  sp_shap$natural_niche_value <- rowSums(st_drop_geometry(sp_shap[,natural_niche_factors]))
  sp_shap$natural_niche_all_positive <- rowSums(st_drop_geometry(sp_shap[,natural_niche_factors]) > 0) == length(natural_niche_factors)
  
  # Summary value 1: percentage of all catchments inside ecological niche 
  val_1 <- c('% sub-catchments inside ecological niche', 
             round(table(sp_shap$natural_niche == T) / nrow(sp_shap), 2)[[2]]*100)
  val_1 <- t(val_1)
  
  # Summary value 2: percentage of all catchments with positive value for all niche shapley values
  val_2 <- c('% sub-catchments with positive shapley values for all natural niche variables', 
             round(table(sp_shap$natural_niche_all_positive) / nrow(sp_shap), 2)[[2]]*100)
  val_2 <- t(val_2)
  
  # Summary value 3: percentage of all catchments positive for each individual niche factor
  val_niche <- sapply(natural_niche_factors, function(x) round(table(sp_shap[[x]] > 0) / nrow(sp_shap), 2)[[2]]*100)
  val_3 <- data.frame(paste0('% all subcatchments with a positive contribution of ', natural_niche_factors), val_niche)
  
  
  #### GENERATE SUMMARY PROPERTIES TO REPORT: THREATENED AREAS INSIDE NICHE
  
  # define sp_shap for only those subcatchments inside the natural niche
  sp_shap_nn <- sp_shap %>% filter(natural_niche == T)
  
  # Summary value 4: percentage of nn catchments negative for each individual threat factor
  val_threat <- sapply(all_threats, function(x) round(table(sp_shap_nn[[x]] < 0) / nrow(sp_shap_nn), 2)[[2]]*100)
  val_4 <- data.frame(paste0('% subcatchments inside niche with a negative contribution of ', all_threats), val_threat)
  
  # Summary values any threats: percentage of nn catchments affected negatively by any threat
  val_threat_any <- data.frame(round(table(rowSums(sapply(all_threats, function(x){sp_shap_nn[[x]] < 0})))/nrow(sp_shap_nn)*100))
  val_threat_any[,1] <- paste0('% of subcatchments inside niche with ', 0:length(all_threats),' threats negative')
  
  # Summary value % sub-catchments with at least 1 threat
  val_threat_1 <- 100 - val_threat_any[1,2]
  val_threat_1 <- t(c('% subcatchments inside niche with at least one negative threat', 
                      val_threat_1))
  # Summary values all threats: percentage of nn catchments affect negatively on average by threats
  val_threat_all <- data.frame(round(table(rowSums(sapply(all_threats, function(x){sp_shap_nn[[x]]})) < 0)/nrow(sp_shap_nn)*100)[[2]])
  val_threat_all <- t(c('% subcatchments inside niche with a negative net effect of all threats summed', 
                        val_threat_all))
  
  
  # Summary value 6: habitat suitability in areas with at least 1 threat
  sp_shap_nn$threatened <- rowSums(sapply(all_threats, function(x){sp_shap_nn[[x]] < 0})) >= 1
  val_6 <- sp_shap_nn %>% 
    group_by(threatened) %>% 
    do(mean_suit = round(mean(.$suitability, na.rm = T), 2)) %>% 
    unnest(cols = c(mean_suit))
  val_6$threatened <- c('mean suitability in unthreatened sub-catchments', 
                        'mean suitability in threatened sub-catchments')
  
  # Summary value 7: how much higher/lower is the suitability in the threatened locations
  val_7 <- t(c('% reduction in suitability in threatened areas within niche', 
               round(as.numeric((val_6[2,2] - val_6[1,2]) / val_6[1,2]),2)*100))
  
  
  # Summary value 8: what % of the expected distribution is absences in threatened vs. unthreatened locations
  table_PA_threat <- table(sp_shap_nn$threatened, !is.na(sp_shap_nn$presence))
  unthreatened_absent <- round(table_PA_threat[1,1]/sum(table_PA_threat[1,1:2]), 2)*100
  val_8a <- t(c('% unthreatened catchments in niche with predicted absence', 
              unthreatened_absent))
  threatened_absent <- round(table_PA_threat[2,1]/sum(table_PA_threat[2,1:2]), 2)*100
  val_8b <- t(c('% threatened catchments in niche with predicted absence', 
                threatened_absent))
  # number of times more absences in threatened compared to unthreatened
  val_8c <- t(c('number of times more absences in threatened compared to unthreatened base rate',
                round(threatened_absent / unthreatened_absent, 2)))
  
  # return all created values together as a long dataframe
  shadow_distribution_summaries <- bind_rows(lapply(list(val_1, val_2, val_3, val_4, 
                                                         val_threat_any, val_threat_1, val_threat_all, val_6, val_7, 
                                                         val_8a, val_8b, val_8c), function(x){
    x <- data.frame(x)
    names(x) <- c('property', 'value'); rownames(x) <- NULL
    lapply(x, as.character)}))
  
  # add in extra information for processing after
  shadow_distribution_summaries$species = species
  shadow_distribution_summaries$property_number = 1:nrow(shadow_distribution_summaries)
  
  ## create output directory
  dir.create(output_folder,recursive = T)
  
  ## save output
  write_csv(shadow_distribution_summaries, 
            file = paste0(output_folder,'/', species, '.csv'))
  
  
  
  #### CREATE SPATIAL OBJECT FOR PLOTTING

  # add in whether habitat factors are negative or positive
  if(length(habitat_factors) == 0){
    sp_shap$neg_habitat <- F
    }else{
      sp_shap$neg_habitat <- rowMeans(st_drop_geometry(sp_shap[habitat_factors])) < 0
      }
  if(length(conn_factors) == 0){
    sp_shap$neg_con <- F
    }else{
    sp_shap$neg_con <- rowMeans(st_drop_geometry(sp_shap[conn_factors])) < 0
    }
  
  # create categorisation of niche and threat combinations in a given subcatchment
  sp_shap$niche_categories <- as.factor(ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == F & sp_shap$neg_habitat == F, '2. inside ecological niche',
                                        ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == T & sp_shap$neg_habitat == F, '3. poor connectivity', 
                                               ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == F & sp_shap$neg_habitat == T, '4. poor habitat', 
                                                      ifelse(sp_shap$natural_niche == T & sp_shap$neg_con == T & sp_shap$neg_habitat == T, '5. poor connectivity and habitat', 
                                                             ifelse(sp_shap$natural_niche == F, '1. outside ecological niche', NA))))))
  
  
  return(sp_shap)
  
}
