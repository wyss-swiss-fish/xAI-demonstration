# script with information defining the RUN for the shapley analysis script to use
# also contains the records table to call containing the species information

# get run array ID here because this is the first script called
i <- as.numeric(commandArgs(trailingOnly = TRUE))
print(i)

RUN_NAME <- 'ubelix_SDM_RF_APRIL_V1_02'

RUN_NAME_APPEND <- '01'

RECORDS_TABLE <- 'records-overview_2010_focal.csv'

N_SIMS <- 10000

source('scripts/shapley/shapley_analysis_UBELIX.R')
