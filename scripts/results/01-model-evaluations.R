#### Supporting figures of model performance ----

# This script reads in the model performance objects created by our random forest SDM pipeline and produces data
# summaries as tables and plots for interpretation.

#### 1. Set directories for loading and saving objects ----

# read in required packages
pacman::p_load(tidyverse, gridExtra)

# local storage
dd <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"
dd_env <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/sdm-pipeline/env-data/"
dd_ch <- "C:/Users/cw21p621/OneDrive - Universitaet Bern/01_Wyss_Academy_for_Nature/analysis/data-dump/"

# get run to mak figures for
RUN <- "ubelix_SDM_RF_MARCH_v6"

# figure directory
fig_dir <- paste0("figures/", RUN, '/')
dir.create(fig_dir, recursive = T)
dir.create(paste0(fig_dir, 'evaluations/'))

# records table for species data summaries
records_table <- read.csv(paste0(dd, 'sdm-pipeline/species-records-final/records-overview_2010.csv'))
sp_list <- unique(records_table$species_name)

# directories holding evaluation data 
eval_dir <- paste0(sdm_dirs, "/output/sdm_evaluations_rf.csv")

#### 2. Join together all evaluation data ----

# read in the evaluation data from the sdm directories
sp_eval <- bind_rows(lapply(eval_dir, function(x) if(file.exists(x)){read.csv(x)}))

# organise the data to subset to the model evaluations that are relevant for our analysis
# here we assess the presence absence data on held out 5-fold spatial cross validations.
sp_eval_sub <- sp_eval %>% 
  # filter to only presence absence
  filter(model_data == 'PA') %>% 
  # filter to evaluation data on PA-all and PA-kfold
  filter(evaluation_data %in% c('PA-kfold')) %>% 
  # remove binary prediction method and focus on TSS and MCC binary predictions
  filter(threshold_method != 'binary_pred')

#### 3. Average model outputs over MCC thresholds and create table ----
# MCC:
# optimized for mcc based on mean values
sp_eval_mcc <- sp_eval_sub %>% 
  # filter threshold method for mcc
  filter(threshold_method == 'mcc', 
         ) %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(mean)) %>% 
  summarise_if(., .predicate = is.numeric, function(x) signif(x, digits = 2)) %>% 
  pivot_longer(., cols = TN:auc, values_to = 'mean')

# optimize for mcc based on sd values
sp_eval_mcc_sd <- sp_eval_sub %>% 
  # filter threshold method for mcc
  filter(threshold_method == 'mcc', 
  ) %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(sd)) %>% 
  summarise_if(., .predicate = is.numeric, function(x) signif(x, digits = 2)) %>% 
  pivot_longer(., cols = TN:auc, values_to = 'sd')

# join outputs together
mcc_summary <- left_join(sp_eval_mcc, sp_eval_mcc_sd) %>% 
  mutate(value = paste0(.$mean, ' (±', .$sd, ')')) %>% 
  select(-mean, -sd) %>% 
  pivot_wider(names_from = species_name, values_from = value)

# save outputs for mcc based thresholds
write.csv(mcc_summary, paste0(fig_dir, 'evaluations/evaluation_mcc.csv'))
 
#### 4. Average model outputs over TSS thresholds and create table ----
# TSS:
# optimized for tss based on mean values
sp_eval_tss <- sp_eval_sub %>% 
  # filter threshold method for mcc
  filter(threshold_method == 'tss') %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(mean)) %>% 
  summarise_if(., .predicate = is.numeric, function(x) signif(x, digits = 2)) %>% 
  pivot_longer(., cols = TN:auc, values_to = 'mean')

# optimize for tss based on sd values
sp_eval_tss_sd <- sp_eval_sub %>% 
  # filter threshold method for mcc
  filter(threshold_method == 'tss' ) %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(sd)) %>% 
  summarise_if(., .predicate = is.numeric, function(x) signif(x, digits = 2)) %>% 
  pivot_longer(., cols = TN:auc, values_to = 'sd')

# join outputs together
tss_summary <- left_join(sp_eval_tss, sp_eval_tss_sd) %>% 
  mutate(value = paste0(.$mean, ' (±', .$sd, ')')) %>% 
  select(-mean, -sd) %>% 
  pivot_wider(names_from = species_name, values_from = value)

# save outputs for mcc based thresholds
write.csv(tss_summary, paste0(fig_dir, 'evaluations/evaluation_tss.csv'))




#### 5. Plot integrative performance metrics against one another (AUC, TSS, MCC) ----

# number of records for use later
pa_records <- records_table %>% select(species_name, pa)

# calculate the average and sd for all metrics when optimized by TSS
tss_mean <- sp_eval_sub %>% 
  # filter threshold method for tss
  filter(threshold_method == 'tss') %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(mean, sd)) %>% 
  summarise_if(., .predicate = is.numeric, function(x) signif(x, digits = 2)) %>% 
  mutate(threshold = 'optimized for tss')

# calculate the average and sd for all metrics when optimized by MCC
mcc_mean <- sp_eval_sub %>% 
  # filter threshold method for mcc
  filter(threshold_method == 'mcc') %>% 
  group_by(species_name, evaluation_data) %>% 
  summarise_at(., .vars = vars(TN, FN, TP, FP, 
                               mcc, overprediction_rate, underprediction_rate, 
                               sorensen, jaccard, 
                               sensitivity, specificity, tss, auc), .fun = list(mean, sd)) %>% 
  summarise_if(., .predicate = is.numeric, list(function(x) signif(x, digits = 2))) %>% 
  mutate(threshold = 'optimized for mcc')

# bind together the relevant measures for plotting
metrics_mean <- left_join(pa_records, bind_rows(tss_mean, mcc_mean)) %>% na.omit()

# quick integrative assessment of which species might show most promising models for further evaluation
metrics_text <- metrics_mean %>% 
  group_by(threshold, species_name) %>% 
  do(mean_measures = (.$mcc_fn1 + .$tss_fn1 + .$auc_fn1) / 3,
     n_pa = .$pa) %>% 
  unnest() %>% 
  filter(n_pa > 75) %>% 
  arrange(threshold, desc(mean_measures)) 

# here we also include Thymallus to ensure we have a more cold-affinity species for comparison
subset_sp <- sort(unique(c(metrics_text$species_name, 'Thymallus thymallus')))

# save subset of species vector for future work
saveRDS(subset_sp, paste0(fig_dir, 'evaluations/subset_sp.RDS'))

# plot of auc vs. mcc
auc_mcc <- ggplot(metrics_mean) + 
  geom_linerange(aes(y = mcc_fn1, xmin = auc_fn1 - auc_fn2, xmax = auc_fn1 + auc_fn2), col = 'gray75') + 
  geom_linerange(aes(x = auc_fn1, ymin = mcc_fn1 - mcc_fn2, ymax = mcc_fn1 + mcc_fn2), col = 'gray75') + 
  geom_point(aes(x = auc_fn1, y = mcc_fn1, size = as.numeric(pa)), col = 'gray50') + 
  geom_point(data = metrics_mean %>% filter(species_name %in% subset_sp), 
             aes(x = auc_fn1, y = mcc_fn1, size = as.numeric(pa), col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(y = mcc_fn1, xmin = auc_fn1 - auc_fn2, xmax = auc_fn1 + auc_fn2, col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(x = auc_fn1, ymin = mcc_fn1 - mcc_fn2, ymax = mcc_fn1 + mcc_fn2, col = species_name)) + 
  facet_wrap(~threshold) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.position = 'right') + 
  scale_size(guide = 'none') + 
  ylim(c(0,1)) + 
  xlim(c(0.5,1)) + 
  geom_abline(aes(slope = 2,
                  intercept = -1)) +
  scale_colour_discrete('') + 
  xlab('AUC') + ylab('MCC')


# plot of auc vs. tss
auc_tss <- ggplot(metrics_mean) + 
  geom_linerange(aes(y = tss_fn1, xmin = auc_fn1 - auc_fn2, xmax = auc_fn1 + auc_fn2), col = 'gray75') + 
  geom_linerange(aes(x = auc_fn1, ymin = tss_fn1 - tss_fn2, ymax = tss_fn1 + tss_fn2), col = 'gray75') + 
  geom_point(aes(x = auc_fn1, y = tss_fn1, size = as.numeric(pa)), col = 'gray50') + 
  geom_point(data = metrics_mean %>% filter(species_name %in% subset_sp), 
             aes(x = auc_fn1, y = tss_fn1, size = as.numeric(pa), col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(y = tss_fn1, xmin = auc_fn1 - auc_fn2, xmax = auc_fn1 + auc_fn2, col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(x = auc_fn1, ymin = tss_fn1 - tss_fn2, ymax = tss_fn1 + tss_fn2, col = species_name)) + 
  facet_wrap(~threshold) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.position = 'right') + 
  scale_size(guide = 'none') + 
  ylim(c(0,1)) + 
  xlim(c(0.5,1)) + 
  geom_abline(aes(slope = 2,
                  intercept = -1)) + 
  scale_colour_discrete('') + 
  xlab('AUC') + ylab('TSS')


# plot of tss vs. mcc
mcc_tss <- ggplot(metrics_mean) + 
  geom_linerange(aes(y = tss_fn1, xmin = mcc_fn1 - mcc_fn2, xmax = mcc_fn1 + mcc_fn2), col = 'gray75') + 
  geom_linerange(aes(x = mcc_fn1, ymin = tss_fn1 - tss_fn2, ymax = tss_fn1 + tss_fn2), col = 'gray75') + 
  geom_point(aes(x = mcc_fn1, y = tss_fn1, size = as.numeric(pa)), col = 'gray50') + 
  geom_point(data = metrics_mean %>% filter(species_name %in% subset_sp), 
             aes(x = mcc_fn1, y = tss_fn1, size = as.numeric(pa), col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(y = tss_fn1, xmin = mcc_fn1 - mcc_fn2, xmax = mcc_fn1 + mcc_fn2, col = species_name)) + 
  geom_linerange(data = metrics_mean %>% filter(species_name %in% subset_sp), 
                 aes(x = mcc_fn1, ymin = tss_fn1 - tss_fn2, ymax = tss_fn1 + tss_fn2, col = species_name)) + 
  facet_wrap(~threshold) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.position = 'right') + 
  scale_size(guide = 'none') + 
  ylim(c(0,1)) + 
  xlim(c(0,1)) + 
  geom_abline(aes(slope = 1,
                  intercept = 0)) + 
  scale_colour_discrete('') + 
  xlab('MCC') + ylab('TSS')


# save plots showing the relationships between model performance measures and highlighting our focal species compared to all modelled species
pdf(file = paste0(fig_dir, 'evaluations/multi_metric_plot.pdf'), height = 10)
cowplot::plot_grid(auc_tss, auc_mcc, mcc_tss, ncol = 1)
dev.off()


#### 6. Plot all model evaluation metrics for our focal species ----

# subset our model evaluation data to focal integrative metrics and focal species
final_plot_data <- sp_eval_sub %>% 
  filter(species_name %in% subset_sp) %>% 
  pivot_longer(., c(overprediction_rate, underprediction_rate, 
                   sensitivity, specificity, 
                   sorensen, mcc, tss, auc))

# clean names of metrics
final_plot_data$name <- factor(final_plot_data$name, 
                               levels = c('overprediction_rate', 'underprediction_rate', 
                                                                'sensitivity', 'specificity', 
                                                                'sorensen', 'mcc', 'tss', 'auc'), 
                               labels = c('Overprediction rate', 'Underprediction rate', 
                                          'Sensitivity', 'Specificity', 
                                          'Sorensen similarity', 'MCC', 'TSS', 'AUC'))

# set up x intercept for plotting
final_plot_data$hline <- ifelse(final_plot_data$name == c('Overprediction rate', 'Underprediction rate'), 0, 1)

# plot all species together as boxplots and save output as pdf
pdf(paste0(fig_dir, 'evaluations/all_metrics_subsetSp.pdf'), width = 10, height = 6)
ggplot(data = final_plot_data) + 
  geom_boxplot(aes(fill = threshold_method, y = species_name, x = value), col = 'black', width = 0.5, lwd = 0.1) +
  facet_wrap(~name, scales = 'free_x', nrow = 2) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1), 
        strip.background = element_blank(), 
        axis.text.y = element_text(face = 'italic')) + 
  xlim(c(0,1)) + 
  geom_vline(aes(xintercept = hline), lty = 2, lwd = 0.2) + 
  scale_fill_discrete('threshold method') + 
  ylab(NULL) + 
  xlab(NULL)
dev.off()

#### 7. Remake model evaluation tables based on focal species ----

write.csv(mcc_summary %>% select(name, subset_sp), 
          paste0(fig_dir, 'evaluations/evaluation_mcc_subset.csv'))
write.csv(tss_summary %>% select(name, subset_sp), 
          paste0(fig_dir, 'evaluations/evaluation_tss_subset.csv'))


#### 8. Make table based on summarising input data to models ----

sdm_data <- read_csv(paste0(dd, 'sdm-pipeline/species-records-final/fish-presenceAbsence_2010.csv'))

summary_data <- sdm_data %>% 
  filter(occ == 1, 
         species_name %in% sp_list) %>% 
  mutate(test = case_when(dataset == 'project_fiumi' ~ '2. EAWAG Projetto Fiumi', 
                          dataset == 'wa_sampling' ~ '1. University of Bern 2022', 
                          dataset == 'vonlanthen' ~ '4. Consultancies', 
                             dataset == 'bafi_data' ~ '3. Kanton', 
                             dataset == 'alte_aare_monitoring' ~ '3. Kanton', 
                             dataset == 'lanat_3_unibe' ~ '1. University of Bern 2022')) %>% 
  group_by(species_name, test) %>% 
  count %>% 
  ungroup() %>% 
  pivot_wider(names_from = test, values_from = n) %>% 
  mutate(total = rowSums(across(where(is.numeric))))

write.csv(summary_data, file = paste0(fig_dir, 'summary_occurrences.csv'))
