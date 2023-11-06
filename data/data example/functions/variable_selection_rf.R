# function to check for correlations and select most explanatory variables for use in random forests

variable_selection_rf <- function(input_data,             # 
                                  cov_names_select,       # all variables to perform selection over
                                  cov_names_priority = NULL,     # variables to force keep
                                  method,                 # one of 'VSURF', 'boruta' for random forests, 'competition'
                                  cor_threshold = NULL    # if method is competition, define threshold'
){
  
  # error handling:
  # check method is correct
  stopifnot(method %in% c('VSURF', 'boruta', 'competition'))
  # check argument cor_threshold is specified when method is competition
  if(method == 'competition' & is.null(cor_threshold)){stop('specify cor_threshold argument when method = competition')}
  # check properties of input_data
  stopifnot(is.data.frame(input_data), 
            c('occ', cov_names_select, cov_names_priority) %in% names(input_data),
            is.character(cov_names_select), 
            is.character(cov_names_priority) | is.null(cov_names_priority))

  ## run different types of variable selection depending on the method chosen 
  if(method == 'VSURF'){
    stopifnot(require('VSURF', character.only = T))
    
    vsurf_out <- VSURF::VSURF(x = input_data[cov_names_select], y = input_data$occ, parallel = T)
    return(vsurf_out)
    
  }
  
  ## boruta
  if(method == 'boruta'){
    stopifnot(require('Boruta', character.only = T))
    
    ### set up weights for downsampled RF
    prNum <- as.numeric(min(table(input_data$occ)/length(input_data$occ))) # number of records in smallest class
    spfrac <- c("0" = prNum, "1" = prNum) # sample size for both classes
    
    boruta_out <- TentativeRoughFix(Boruta::Boruta(x = input_data[cov_names_select], 
                                 y = input_data$occ, 
                                 ntree = 1000, 
                                 replace = T, 
                                 sample.fraction = spfrac))
    
    return(boruta_out)
    
  }
  
  # variable selection method based on assessing performance of correlated covariates
  if(method == 'competition'){
    
    # check for correlation amongst variables 
    cor_vars <- cor(apply(input_data[cov_names_select], 2, function(x) as.numeric(as.character(x))))
    cor_vars[!upper.tri(cor_vars)] <- NA
    correlated <- data.frame(which(abs(cor_vars) > cor_threshold, arr.ind = T))
    
    while(nrow(correlated) != 0){
      var_pairs <- list()
      for(i in 1:nrow(correlated)){
        correlated_names <- correlated[i,]
        row_var <- rownames(cor_vars)[correlated_names[1,1]]
        col_var <- colnames(cor_vars)[correlated_names[1,2]]
        var_pairs[[i]] <- c(row_var, col_var)
      }
      
      # fit models competing variable pairs
      selected_var <- c()
      for(i in 1:length(var_pairs)){
        
        # if one variable in the pair is a 'keep' variable then select this one
        if(sum(var_pairs[[i]] %in% cov_names_priority) == 1){ # if this = 0 or equals 2 then they should fight!
          selected_var[i] <- which(var_pairs[[i]] %in% cov_names_priority)
        }else{
          
          
          # perform selection for rf
            # create rf data
            rf_data  <- data.frame(occ = input_data['occ'][,1], input_data[var_pairs[[i]][1:2]][,1:2])
            ### downsampled RF
            prNum <- as.numeric(min(table(rf_data$occ))) # number of records in smallest class
            spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes
            # fit rf_1
            rf_1 <- tryCatch(randomForest(occ ~ .,
                                          data = data.frame(occ = rf_data$occ, rf_data[var_pairs[[i]][1]]),
                                          ntree = 1000,
                                          sampsize = spsize,
                                          replace = TRUE), error = function(e) NA)
            if(sum(!is.na(rf_1)) != 0){
              
              rf_predict <- as.numeric(as.character(predict(rf_1, newdata=rf_data[var_pairs[[i]][1]], type = 'prob')[,2]))
              error_1    <- max_mcc(rf_predict, rf_data$occ)[[2]]
              
            }
            
            # fit rf_2
            rf_2 <- tryCatch(randomForest(occ ~ .,
                                          data = data.frame(occ = rf_data$occ, rf_data[var_pairs[[i]][2]]),
                                          ntree = 1000,
                                          sampsize = spsize,
                                          replace = TRUE), error = function(e) NA)
            if(sum(!is.na(rf_2)) != 0){
              rf_predict <- as.numeric(as.character(predict(rf_2, newdata=rf_data[var_pairs[[i]][2]], type = 'prob')[,2]))
              error_2    <- max_mcc(rf_predict, rf_data$occ)[[2]]
            }
            # select optimal covariate
            if(is.na(error_1)){selected_var[i] <- 2}
            if(is.na(error_2)){selected_var[i] <- 1}
            if(sum(is.na(c(error_1, error_2)))==0){selected_var[i] <- which.max(c(error_1, error_2))}
          
        } # end of else statement for forcing selection of 'keep' variables
      }  # end of loop through all variable pairs
      
      
      # selected variables using d-squared and error rate
      selected_var_2 <- unique(unlist(lapply(1:length(var_pairs), function(x) var_pairs[[x]][selected_var[x]])))
      
      # find all the variables that have no correlations
      cor_vars_all <- cor(apply(input_data[cov_names_select], 2, function(x) as.numeric(as.character(x))))
      cor_vars_all[cor_vars_all==1]<-NA
      var_no_cor <- c()
      for(i in 1:nrow(cor_vars_all)){
        var_no_cor[i] <- sum(na.omit(abs(cor_vars_all[i,])) > cor_threshold) == 0
      }
      var_no_cor <- rownames(cor_vars_all)[var_no_cor]
      
      # combine with the selected variables that have correlations and return as a covariate object
      cov_names_select <- c(selected_var_2, var_no_cor)
      cor_vars <- cor(apply(input_data[cov_names_select], 2, function(x) as.numeric(as.character(x))))
      cor_vars[!upper.tri(cor_vars)] <- NA
      correlated <- data.frame(which(abs(cor_vars) > cor_threshold, arr.ind = T))
      
    }
    
    
    

      return(cov_names_select)

  }
  
}
