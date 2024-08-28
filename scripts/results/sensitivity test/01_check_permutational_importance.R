#### R script to contrast permutational variable importance with SHAP variable importance

#### Read in final random forest models ----

# read in required packages
pacman::p_load(tidyverse, gridExtra, randomForest, DALEX)

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_APRIL_V1_02"
RUN_SDM <- "ubelix_SDM_RF_APRIL_V1"

# figure directory
fig_dir <- paste0("figures/", RUN, '/')
dir.create(fig_dir, recursive = T)
dir.create(paste0(fig_dir, 'perm_importance/'))

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

#### example for single species ----

all_importances <- lapply(1:length(sp_list), function(x){
     rf1 <- readRDS(all_rfs[[x]])
     data1 <- readRDS(all_data[[x]])
     shap1 <- readRDS(shap_pa[[x]])
     train1 <- data1@pa_data$full_data
     occ <- as.numeric(train1[,1])
     train1 <- train1[,-c(1:3)]
     
     # get permutational importance
     gini_imp <- rf1$importance
     gini_imp <- gini_imp %>% data.frame(Variable = rownames(gini_imp), .)
     
     # use of DALEX explainer
     explainer_rf <- DALEX::explain(model = rf1, 
                                    data = train1, 
                                    y = occ, 
                                    label = "Random Forest")
     
     
     perm_imp <- model_parts(explainer = explainer_rf, 
                 loss_function = loss_one_minus_auc, 
                 type = 'difference', 
                 N=NULL, 
                 B=10) %>% 
       data.frame()
     
     perm_imp <- perm_imp %>% 
       rename(Variable = variable) %>% 
       group_by(Variable) %>% 
       do(mean_dropout_loss = mean(.$dropout_loss)) %>% 
       unnest(mean_dropout_loss) 
       
     # shapley importance
     shap_imp <- data.frame(Variable = gini_imp$Variable, shap_mean = colMeans(abs(shap1[paste0(gini_imp$Variable, '_SHAP')])))
     
     # join all
     all_imp<- plyr::join_all(list(gini_imp, shap_imp, perm_imp))
     return(all_imp)
})

names(all_importances) <- sp_list
all_importances <- bind_rows(all_importances, .id = 'species_name')

png(paste0(fig_dir, '/perm_importance/gini_vs_shap_comparison.png'), res = 300, height = 2000, width = 2000)
ggplot(data = all_importances) + 
  geom_point(aes(x = MeanDecreaseGini, y = shap_mean)) + 
  stat_smooth(aes(x = MeanDecreaseGini, y = shap_mean, group = species_name), se= F, method = 'lm') + 
  facet_wrap(~species_name, scales = 'free') + 
  xlab('Variable importance as Gini index') + 
  ylab('SHAP variable importance as mean absolute SHAP value') +
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()


png(paste0(fig_dir, '/perm_importance/perm_vs_shap_comparison.png'), res = 300, height = 2000, width = 2000)
ggplot(data = all_importances) + 
  geom_point(aes(x = mean_dropout_loss, y = shap_mean)) + 
  stat_smooth(aes(x = mean_dropout_loss, y = shap_mean, group = species_name), se= F, method = 'lm') + 
  facet_wrap(~species_name, scales = 'free') + 
  xlab('Permutational variable importance as loss in AUC') + 
  ylab('SHAP variable importance as mean absolute SHAP value') +
  theme_bw() + 
  theme(aspect.ratio = 1)
dev.off()


cors <- all_importances %>% 
  group_by(species_name) %>% 
  do(cor_gini_shap = cor(.$MeanDecreaseGini, .$shap_mean, method = 'spearman'), 
     cor_perm_shap = cor(.$mean_dropout_loss, .$shap_mean, method = 'spearman')) %>% 
  unnest()
mean(cors$cor_gini_shap)
sd(cors$cor_gini_shap)

mean(cors$cor_perm_shap)
sd(cors$cor_perm_shap)

cors <- all_importances %>% 
  filter(species_name != 'Oncorhynchus mykiss') %>% 
  group_by(species_name) %>% 
  do(cor_gini_shap = cor(.$MeanDecreaseGini, .$shap_mean, method = 'pearson'), 
     cor_perm_shap = cor(.$mean_dropout_loss, .$shap_mean, method = 'pearson')) %>% 
  unnest()

mean(cors$cor_gini_shap)
mean(cors$cor_perm_shap)
