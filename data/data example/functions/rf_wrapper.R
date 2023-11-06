rf_wrapper <- function(x){
  
  # ensure occ is a factor
  x$occ <- factor(x$occ)
  
  #### Select covariates for down-sampled random forests:
  rf_vars <- variable_selection_rf(
    input_data = x,
    cov_names_select = names(env_data),
    cov_names_priority = names(env_data),
    method = var_selection_method
  )
  
  # process variable selection output depending on 'var_selection_method'
  if (var_selection_method == "boruta") {
    rf_vars <- names(rf_vars$finalDecision[rf_vars$finalDecision == "Confirmed"])
  }
  
  # make properties of the training data for saving with evaluations
  n_occ_training <- sum(x$occ == 1)
  n_abs_training <- sum(x$occ == 0)
  prevelence_training <- n_occ_training / n_abs_training
  prNum <- as.numeric(min(table(x$occ))) # number of records in smallest class
  spsize <- c("0" = prNum, "1" = prNum)  # sample size for both classes
  
  
  ### FIT RANDOM FORESTS EVALUATIONS
  # fit the optimized model to the cross validation subsets
  rf_optimized <- randomForest(occ ~ .,
                               data     = x[c("occ", rf_vars)],
                               ntree    = 1000,
                               sampsize = spsize,
                               replace  = TRUE
  ) # make sure samples are with replacement (default)
  
  return(rf_optimized)
  
}