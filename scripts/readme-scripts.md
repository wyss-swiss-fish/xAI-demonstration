# read me

To run our analysis we use the scripts in the following order:
- initiate sbatch_runShapley_ubelix.sh on UBELIX server at University of Bern.
- this runs the script shapRun_main.R (or sensitivity checks depending on file called). This file sets up some parameters for the shapley runs (sets up the name of the input folders from the SDM runs, an output name, and the number of simulations)
- the shapRun_main.R script runs the code shapley_analysis_UBELIX.R script which is the main script running the shapley analysis. 