
# read in required packages
pacman::p_load(tidyverse, gridExtra, randomForest, DALEX, sf, terra)

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"

# get species
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# get directories for random forest objects
sdm_dirs <- list.files(paste0("D:/sdm-pipeline/sdm-run/", RUN_SDM), full.names = T)
sdm_dirs <- sdm_dirs[grepl(paste0(sp_list, collapse = "|"), sdm_dirs)]
all_rfs  <- paste0(sdm_dirs, '/output/final_models_rf_pa.RDS')

# get directories for focal shapley objects
shap_dirs <- list.files(paste0("shapley-run/", RUN),
                        full.names = T
)
shap_dirs <- shap_dirs[grepl(paste0(sp_list, collapse = "|"), shap_dirs)]
shap_pa <- paste0(shap_dirs, "/shapley_rf_pa.RDS")

# get data inputs across all species
all_data <- paste0(sdm_dirs, '/data/sdm_input_data.rds')




#### Read in random forest objects ----

## LOAD OCCURRENCE AND ENVIRONMENTAL DATA

# occurrence dataset
pa <- st_read("scripts/workflow example/data example/pa.shp", quiet  = T)

# environmental dataset
env_data <- rast("scripts/workflow example/data example/env_data.tif")

# complete combined dataset
full_data <- readRDS(file = "scripts/workflow example/data example/sp_example.rds")

# load pre-made spatial objects for plotting
load(file = "scripts/workflow example/data example/all_spatial_Waldock2023.RData")

# load functions
funs <- lapply(list.files("scripts/workflow example/data example/functions", full.names = T),
               function(x) source(x, echo = F))


# Generate species distribution model 
var_selection_method = "boruta"
pa_rf_final <- rf_wrapper(full_data)
pa_rf_final

# rerun with standardized data
full_data_std <- full_data
full_data_std[names(env_data)] <- apply(full_data[names(env_data)], 2, scales::rescale)
pa_rf_final_std <- rf_wrapper(full_data_std)
pa_rf_final_std



## apply SHAP
# load in the fast shap package used to calculate shapley values
pacman::p_load(fastshap)

# define the prediction function to use in fastshap
pfun <- function(object, newdata) {
  as.numeric(as.character(predict(object, newdata = newdata, 
                                  type = "prob")[,2]))
}

# get the variables in a way that is model specific
vars <- colnames(attr(pa_rf_final$terms, "factors"))

# check terms are the same
identical(colnames(attr(pa_rf_final$terms, "factors")), colnames(attr(pa_rf_final_std$terms, "factors")))


shapley_pa     <- fastshap::explain(object = pa_rf_final,     feature_names = vars, X = full_data[vars],     pred_wrapper = pfun,  nsim = 500)
shapley_pa_std <- fastshap::explain(object = pa_rf_final_std, feature_names = vars, X = full_data_std[vars], pred_wrapper = pfun,  nsim = 500)

shap_comp <- cbind(shapley_pa %>% 
            data.frame() %>% 
            pivot_longer(1:ncol(.)) %>% 
            rename(variable_raw = value), 
            shapley_pa_std %>% 
            data.frame() %>% 
            pivot_longer(1:ncol(.)) %>% 
            rename(variable_std = value) %>% 
              select(-name))

ggplot(data = shap_comp) + 
  geom_point(aes(x = variable_raw, y = variable_std)) + 
  ylab('standardized variables') + 
  xlab('unstandardized variables')

cor(shap_comp$variable_raw, shap_comp$variable_std)



