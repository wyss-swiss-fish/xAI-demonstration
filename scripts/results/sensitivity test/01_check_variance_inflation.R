
#### R script to check coefficient terms using simpler GLMs

#### Read in final random forest models ----

# read in required packages
pacman::p_load(tidyverse, gridExtra)

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
dir.create(paste0(fig_dir, 'glm_check/'))

# get species
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- readRDS(paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# get directories for random forest objects
sdm_dirs <- list.files(paste0("D:/sdm-pipeline/sdm-run/", RUN_SDM), full.names = T)
sdm_dirs <- sdm_dirs[grepl(paste0(sp_list, collapse = "|"), sdm_dirs)]
all_rfs  <- paste0(sdm_dirs, '/output/final_models_rf_pa.RDS')

# get data inputs across all species
all_data <- paste0(sdm_dirs, '/data/sdm_input_data.rds')


#### variable names ----

# vars
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
          'local_flood', 
          'I(stars_t_mx_m_c^2)')

# rename variables
vars_renamed = c('discharge', 
                 'slope', 
                 'flow velocity', 
                 'temperature max', 
                 'temperature min', 
                 'connectivity', 
                 'distance to lake', 
                 'morph. mod.', 
                 'urbanisation',
                 'wetland', 
                 'floodplains', 
                 'temperature max ^2')

vars_renamed <- cbind(vars, vars_renamed)




#### run loop to fit simple GLM for all species -----

pacman::p_load(tidymodels)


glms <- lapply(1:length(all_data), function(x){
  
  # get data and rf
  data1 <- readRDS(all_data[[x]])
  train1 <- data1@pa_data$full_data
  rf1 <- readRDS(all_rfs[[x]])
  vars_rf <- colnames(attr(rf1$terms, "factors"))
  vars_subset <- c('ecoF_discharge_max_log10', 'local_asym_cl_log10', 'stars_t_mx_m_c', 'local_flood', 'local_dis2lake', 'ecoF_flow_velocity_mean')
  vars_subset <- vars_subset[which(vars_subset %in% vars_rf)]
  train1[vars_subset] <- apply(train1[vars_subset], 2, scales::rescale)

  if('stars_t_mx_m_c' %in% vars_subset){
  form = formula(paste("occ ~ ",paste0(vars_subset,collapse = "+"), " + I(stars_t_mx_m_c^2)"))
  }
  
  mod <- glm(form, data = train1, family = 'binomial')
  print(glance(mod))
  return(tidy(mod))
  
})

names(glms) <- sp_list

glms <- bind_rows(glms, .id = 'species_name')

glms <- left_join(glms, data.frame(vars_renamed), by = c('term' = 'vars'))

png(paste0(fig_dir, 'glm_check/coefficients_glm.png'), res = 300, width = 2500, height = 1500)
ggplot(data = glms %>% filter(term != '(Intercept)')) + 
  geom_point(aes(x = estimate, y = species_name, col = p.value<0.05), size = 3) + 
  geom_linerange(aes(y = species_name, xmin = estimate-std.error, xmax = estimate+std.error, col = p.value<0.05)) + 
  geom_vline(aes(xintercept = 0), lty = 2) + 
  facet_wrap(~vars_renamed, scales = 'free_x', nrow = 2) + 
  ylab(NULL) + 
  theme_bw()
dev.off()
  
#### get correlation matrix ----

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# make correlations from terra object used to fit model
test_cor <- terra::layerCor(env_data, fun = 'pearson', n = 5000)

# load melting function
pacman::p_load('reshape2')
melt_cor <- na.omit(melt(test_cor[[1]]))

melt_cor_2 <- left_join(left_join(melt_cor, data.frame(vars_renamed), by = c('Var1' = 'vars')), 
                    data.frame(vars_renamed), by = c('Var2' = 'vars'))
melt_cor_2 <- melt_cor_2 %>% arrange(Var1, Var2)
melt_cor_2 %>% 
  filter(!vars_renamed.y == c('temperature min') & !vars_renamed.x == c('temperature max')) %>% 
  filter(!vars_renamed.y == c('temperature max') & !vars_renamed.x == c('temperature min')) %>% 
  filter(value != 1) %>% 
  pull(value) %>% 
  abs %>% 
  quantile(.,c( 0.5, 0.95, 1))

# plot
png(paste0(fig_dir, 'check_env_cor.png'), res = 300, width = 3000, height = 3000)
ggplot(data = melt_cor_2, aes(vars_renamed.y, vars_renamed.x, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours=c('indianred4','indianred3', 'gray90','white','gray90','darkseagreen2', 'darkseagreen4'),
                       limits = c(-1, 1), 
                       breaks=c(-1, -0.5,  0, 0.5, 1),
                       values = scales::rescale(c(-1, -0.5,-0.49, 0, 0.49, 0.5, 1)),
                       name="Pearson\nCorrelation", 
                       scale = T) +
  theme_minimal()+ 
  geom_text(data = melt_cor_2, aes(y = vars_renamed.y, x = vars_renamed.x, 
                                   label = round(value,2)), size = 3) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10))+
  coord_fixed() + 
  xlab(NULL) + 
  ylab(NULL)
dev.off()
