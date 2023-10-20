## function to summarise the output of the generate_SD_output into a quantitative shadow distribution under different assumptions

generate_quant_SD <- function(generate_SD_object){
  
  ####
  ## READ IN DATA
  # get species i 
  sp <- generate_SD_object$species_name[1]
  print(sp)
  
  # get shapley shape files distribution
  shap_sp_i <- generate_SD_object
  
  
  ####
  ## DEFINE THE SUITABILITY PREDICTIONS USING THE BASELINE PREDICTION AND THE SUMMED SHAPLEY VALUES
  
  # get the baseline value as the mean of all habitat suitability predictions
  baseline_value <- mean(shap_sp_i$suitability, na.rm = T)
  
  # sum all the shapley values to get the net prediction as a deviaiton from the baseline mean prediction
  shap_sp_i$shap_all_sum <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% all_factors]), na.rm = T)
  
  # get the baseline suitability values estimated as the sum of the shapely + baseline
  # this ensure internal consistency and that all values are estimated from the shapely values
  shap_sp_i$shap_suit_baseline <- shap_sp_i$shap_all_sum + baseline_value
  # the shap_suit_baseline object is also the predicted distribution
  
  # check the correlation is very high between methods, some error is expected from randomisation proceedure in the shapley estimations.
  suit_method_correlation <- cor(shap_sp_i$shap_suit_baseline, shap_sp_i$suitability, method = 'pearson')
  print(suit_method_correlation)
  
  
  ####
  ## DEFINE THE EXPECTED DISTRIBUTION USING THE BASELINE PREDICTION AND THE POSITIVE SHAPLEY VALUES, RELEASING THREATS
  ## here we require 1. the areas with positive values for natural factors, and within 
  ## this set, the 2. areas with positive values for threats (unthreatened), 3. areas that would be suitable if
  
  ## PARTITION
  # partition out the different shapley components in different groupings
  shap_sp_i$nn_sum      <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% nn_factors]), na.rm = T)
  shap_sp_i$nn_sum_mask <- ifelse(shap_sp_i$nn_sum < 0, NA, shap_sp_i$nn_sum)
  shap_sp_i$hab_sum     <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% hab_factors]), na.rm = T)
  shap_sp_i$con_sum     <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% con_factors]), na.rm = T)
  shap_sp_i$threat_continuous <- rowSums(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)]), na.rm = T)
  
  ## EXPECTED DISTRIBUTION
  # get positive shapley values for threats (i.e., where habitat are positive contributions) 
  
  # different methods for estimating corrected values
  # method 1: convert to absolute value
  # method 2: convert to the mean value where the shapley values are postive
  # method 3: convert to the maximum value where the shapley values are positive
  # method 4: convert to 0 (remove threat but not any positive benefit of removal)
  
  threat_shaps <- st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])
  threat_shaps[threat_shaps<0] <- NA
  
  if(method == 1){
    # method 1
    corrected_threat_shaps <- rowSums(abs(st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])), na.rm = T)
  }
  
  if(method == 2){
    # method 2 
    corrected_threat_shaps <- rowSums(apply(threat_shaps, 2, function(x){
      x[is.na(x)] <- mean(x, na.rm = T)
      return(x)}))
  }
  
  if(method == 3){
    # method 3
    corrected_threat_shaps <- rowSums(apply(threat_shaps, 2, function(x){
      x[is.na(x)] <- quantile(x, 0.95, na.rm = T)
      return(x)}))
  }
  
  if(method == 4){
    # method 4
    corrected_threat_shaps <- rowSums(apply(threat_shaps, 2, function(x){
      x[is.na(x)] <- 0
      return(x)}))
  }
  
  # add back in positive effects of threats to get the expected distribution
  shap_sp_i$expected_distribution <- shap_sp_i$nn_sum + corrected_threat_shaps + baseline_value
  
  # mask the expected distribution 
  shap_sp_i$expected_distribution <- ifelse(shap_sp_i$nn_sum < 0, NA, shap_sp_i$expected_distribution)
  
  # convert to 1 if > 
  shap_sp_i$expected_distribution[shap_sp_i$expected_distribution > 1] <- 1
  
  ## OBSERVED DISTRIBUTION
  # take the habitat suitability predictions within the natural niche area
  shap_sp_i$observed_distribution <- ifelse(is.na(shap_sp_i$nn_sum_mask), NA, shap_sp_i$shap_suit_baseline)
  
  
  ## CONTRIBUTIONS TO SHADOW DISTRIBUTIONS
  # get matrix of negative shapley values
  any_neg_threat_shaps <- st_drop_geometry(shap_sp_i[,names(shap_sp_i) %in% c(hab_factors, con_factors)])
  any_neg_threat_shaps[any_neg_threat_shaps >= 0] <- NA
  
  # get the areas where any threat is negative
  shap_sp_i$sum_presence_negative_threat <- as.numeric(apply(any_neg_threat_shaps < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_con <- as.numeric(apply(any_neg_threat_shaps[names(any_neg_threat_shaps) %in% c(con_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  shap_sp_i$sum_presence_negative_habitat <- as.numeric(apply(any_neg_threat_shaps[names(any_neg_threat_shaps) %in% c(hab_factors)] < 0, 1, function(x) sum(x, na.rm = T)))
  
  # define the shadow distribution quantitatively as the observed / expected
  shap_sp_i$SD_OratioE <- shap_sp_i$observed_distribution / shap_sp_i$expected_distribution
  
  # define shadow distribution as the observed as a proportion of E
  shap_sp_i$SD_OpercentOfE <- (shap_sp_i$observed_distribution - shap_sp_i$expected_distribution) / shap_sp_i$expected_distribution
  
  # define the shadow distribution quantitatively as the observed - expected
  shap_sp_i$SD_OminusE <- shap_sp_i$observed_distribution - shap_sp_i$expected_distribution
  
  # define the shadow distribution as the areas in the expected distribution that have negative net threat effects
  shap_sp_i$SD_expected_threat_continuous <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$threat_continuous)
  
  # define the shadow distribution as the areas in the expected where any threat is negative, and the magnitude of this threat
  shap_sp_i$SD_expected_sum_presence_negative_threat <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_threat)
  shap_sp_i$SD_expected_sum_presence_negative_con <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_con)
  shap_sp_i$SD_expected_sum_presence_negative_habitat <- ifelse(is.na(shap_sp_i$expected_distribution), NA, shap_sp_i$sum_presence_negative_habitat)
  
  # test maps of species different distribution types
  # tm_shape(shap_sp_i) + tm_fill(col = 'shap_suit_baseline')
  # tm_shape(shap_sp_i) + tm_fill(col = 'observed_distribution')
  # tm_shape(shap_sp_i) + tm_fill(col = 'expected_distribution')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_OminusE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_ThreatInE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'shadow_distribution_anyNegThreatInE')
  # tm_shape(shap_sp_i) + tm_fill(col = 'expected_but_absent')
  # tm_shape(shap_sp_i) + tm_fill(col = 'observed_and_present')
  
  # NOTE:
  # the observed_distribution and expected_distribution have the same masking, 
  # such that the observed is simply the expected but with the effects of threat factors retained in the
  # habitat suitability score. In this way, we can directly compare these layers as they are both built from the 
  # partitioned suitability scores using the shapley values.
  # Note that we avoided metrics that use the threshold defined by the original suitability model. This value is inherently 
  # biased towards representing all the processes that affect a species distribution together, whereas our 
  # comparison of observed_distribution and expected_distribution are built from partitioning the contributions
  # based on different categories so is less biased towards the currently observed distribution, but we lack a clear way 
  # to define a 'presence' or 'absence' based on thresholds as is often reported in the literature. 
  
  # return the objects of interest
  return(shap_sp_i)
}
  
  