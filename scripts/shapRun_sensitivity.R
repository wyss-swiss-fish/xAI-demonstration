# script with information defining the RUN for the shapley analysis script to use
# also contains the records table to call containing the species information


### Subset 1 - in house data only
# get run array ID here because this is the first script called
i <- as.numeric(commandArgs(trailingOnly = TRUE))
print(i)

RUN_NAME <- 'ubelix_SDM_RF_JULY_inHouse_V2'
RUN_NAME_APPEND <- '01'
RECORDS_TABLE <- 'records-overview_2010_focal.csv'
N_SIMS <- 1000
source('scripts/shapley/shapley_analysis_UBELIX.R')
rm(list = ls())

### Subset 2 - kantonal data only
# get run array ID here because this is the first script called
i <- as.numeric(commandArgs(trailingOnly = TRUE))
print(i)

RUN_NAME <- 'ubelix_SDM_RF_JULY_Kanton_subset_V1'
RUN_NAME_APPEND <- '01'
RECORDS_TABLE <- 'records-overview_2010_focal.csv'
N_SIMS <- 1000
source('scripts/shapley/shapley_analysis_UBELIX.R')
rm(list = ls())

### Subset 3 - all data but removing high-elevation disconnected sites
# get run array ID here because this is the first script called
i <- as.numeric(commandArgs(trailingOnly = TRUE))
print(i)

RUN_NAME <- 'ubelix_SDM_RF_JULY_Kanton_V1'
RUN_NAME_APPEND <- '01'
RECORDS_TABLE <- 'records-overview_2010_focal.csv'
N_SIMS <- 1000
source('scripts/shapley/shapley_analysis_UBELIX.R')
