# read me

To run our analysis we use the scripts in the following order:
sbatch_runShapley_ubelix.sh -> shapRun_main.R -> shapley_analysis_UBELIX.R -> saved output files

- Initiate sbatch_runShapley_ubelix.sh on UBELIX server at University of Bern.

- The above script will run script shapRun_main.R (or sensitivity checks depending on file called). This file sets up some parameters for the shapley runs (sets up the name of the input folders from the SDM runs, an output name, and the number of simulations)

- The shapRun_main.R script runs the code shapley_analysis_UBELIX.R script which is the main script running the shapley analysis. At the end of this script is a file saved as an .RDS file containing the shapley values.