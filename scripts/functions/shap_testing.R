# script to test Shapley value estimates

get_shapleys <- function(model,    # random forest model fitted to data
                         data,     # data used to fit the random forest
                         new_data, # new data on which to predict the shapley values
                         vars      # variables used in model
                         ){
  
  # if not available install fastshap
  # install.packages('fastshap')
  
  # define X
  X <- data[vars]
  
  # define the prediction function to use
  pfun <- function(object, newdata) {
    as.numeric(as.character(predict(object, newdata = newdata, type = 'prob')[,2]))
  }
  
  # get the shapley values
  shapley <- fastshap::explain(object = model, 
                                     feature_names = vars, 
                                     X = X, 
                                     newdata = new_data[vars],
                                     pred_wrapper = pfun, 
                                     nsim = 1)
  
  # aggregate and clean up
  shapley <- shapley %>% data.frame %>% mutate_if(., is.numeric, .funs = round, 2)
  
  # rename shapley data
  names(shapley) <- paste0(names(shapley), 'SHAP')
  
  # bind back with the new_data
  shapley <- cbind(new_data, shapley)
  
  return(shapley)

}



## install.packages('fastshap')
## library(fastshap)
## library(randomForest)
## 
## model = 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/output/final_rf.rds'
## model <- readRDS(model)[[1]]
## 
## data = 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/data/model_data_final.rds'
## final_data = readRDS(data)[[1]]
## rm(data)
## # this needs to be generalizable or saved
## final_vars <- colnames(attr(model$terms,"factors"))
## X <- final_data[final_vars]
## 
## pfun <- function(object, newdata) {
##   as.numeric(as.character(predict(object, newdata = newdata, type = 'prob')[,2]))
## }
## 
## # get the shapley values
## test_fastshap <- fastshap::explain(object = model, 
##                                    feature_names = final_vars, 
##                                    X = X, 
##                                    pred_wrapper = pfun, 
##                                    nsim = 250)
## 
## # Aggregate Shapley values
## shap_imp <- data.frame(
##   Variable = names(test_fastshap),
##   Importance = apply(test_fastshap, MARGIN = 2, FUN = function(x) sum(abs(x)))
## )
## 
## # Plot Shap-based variable importance
## ggplot(shap_imp, aes(reorder(Variable, Importance), Importance)) +
##   geom_col() +
##   coord_flip() +
##   xlab("") +
##   ylab("mean(|Shapley value|)")
## 
## # Plot shapley values against variable values
## grid.arrange(grobs = lapply(final_vars, function(x){
##   shap_dep_x3 <- data.frame(x3 = X[[x]], shap = test_fastshap[[x]])
##   return(ggplot(shap_dep_x3, aes(x3, shap)) +
##     geom_point(alpha = 0.3) +
##     geom_smooth(method = 'lm') +
##     geom_smooth(method = 'gam') +
##     ylab("Shapley value") + 
##     xlab(x)) 
## }))
## 
## 
## #### Try with spatial data ----
## 
## # read in environmental raster used to fit the model
## env_data = 'C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/sdm-run/august_01/species/Telestes souffia/data/env_data.TIF'
## env_data <- rast(env_data)
## 
## # load in 40km catchments
## subcatchment_file <- paste0(dd, 'swiss-2km-subcatchments/EZG_Gewaesser.gdb')
## subcatchments <- st_read(subcatchment_file, layer = 'Teileinzugsgebiet')
## subcatchments <- st_transform(subcatchments, crs = target_crs)
## subcatchments <- st_cast(subcatchments, 'POLYGON') %>% st_zm()
## 
## # union the subcatchments into 40km catchments
## subcatchments_40kms <- subcatchments %>%
##   mutate(TEZGNR40 = as.character(TEZGNR40)) %>% 
##   group_by(TEZGNR40) %>% 
##   dplyr::summarize(Shape = st_union(Shape))
## 
## #mapview::mapview(subcatchments_40kms)
## 
## # extract the environmental data values for the modelled variables
## env_poly <- terra::extract(x = env_data, y = vect(subcatchments_40kms), fun = function(x){mean(x, na.rm = T)})
## 
## # bind back in environmental data from extractions
## subcatchments_40kms_env <- cbind(subcatchments_40kms, env_poly)
## 
## subcatchments_40kms_env_df <- st_drop_geometry(subcatchments_40kms_env) %>% na.omit()
## 
## #### Test shapley values with new data ----
## 
## # get the shapley values
## test_fastshap <- fastshap::explain(object = model, 
##                                    feature_names = final_vars, 
##                                    X = X, 
##                                    newdata = subcatchments_40kms_env_df[final_vars],
##                                    pred_wrapper = pfun, 
##                                    nsim = 500)
## 
## test_fastshap <- test_fastshap %>% data.frame %>% mutate_if(., is.numeric, .funs = round, 2)
## 
## # rename shapley data
## names(test_fastshap) <- paste0(names(test_fastshap), 'SHAP')
## 
## # bind with environmental values
## shap_bind <- cbind(subcatchments_40kms_env_df, test_fastshap)
## 
## # bind in predictions
## shap_bind$preds <- pfun(model, shap_bind[final_vars])
## 
## # link with spatial sf object
## shap_spatial <- left_join(shap_bind, subcatchments_40kms_env %>% na.omit()) %>% st_as_sf()
## 
## final_SHAP <- names(shap_spatial)[grepl('SHAP', names(shap_spatial))]
## 
## shap_tmap <- list()
## for(SHAP in final_SHAP){
##   shap_tmap[[SHAP%in%final_SHAP]] <- tm_shape(shp = shap_spatial) +
##           tm_polygons(col = SHAP, style = 'cont', palette = 'RdBu')
## }
## 
## tmap_arrange(shap_tmap)
## 
## library(mapview)
## mapview(shap_spatial[c(names(test_fastshap), 'preds')] %>% filter(preds < 0.5), zcol = 'ecoF_discharge_maxSHAP')
## 
## 
## ### make long and look at upvoting and downvoting
## shap_positive <- shap_spatial %>% 
##   st_drop_geometry() %>% 
##   dplyr::select(TEZGNR40, preds, names(test_fastshap)) %>% 
##   pivot_longer(., cols = names(test_fastshap)) %>% 
##   group_by(TEZGNR40, preds) %>% 
##   do(var_pos = .$name[which.max(.$value)],
##      value_max = .$value[which.max(.$value)], 
##      sign_max  = sign(.$value[which.max(.$value)])) %>% 
##   filter(sign_max == 1, abs(value_max) > 0.05) %>% 
##   unnest() %>% 
##   left_join(shap_spatial %>% dplyr::select(TEZGNR40), .)
## 
## shap_negative <- shap_spatial %>% 
##   st_drop_geometry() %>% 
##   dplyr::select(TEZGNR40, preds, names(test_fastshap)) %>% 
##   pivot_longer(., cols = names(test_fastshap)) %>% 
##   group_by(TEZGNR40, preds) %>% 
##   do(var_neg = .$name[which.min(.$value)],
##      value_min = .$value[which.min(.$value)],
##      sign_min  = sign(.$value[which.min(.$value)])) %>% 
##   filter(sign_min == -1, abs(value_min) > 0.05) %>% 
##   unnest() %>% 
##   left_join(shap_spatial %>% dplyr::select(TEZGNR40), .)
## 
## mapview::mapview(shap_positive, zcol = 'var_pos')
## mapview::mapview(shap_negative, zcol = 'var_neg')

