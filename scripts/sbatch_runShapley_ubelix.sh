#! /bin/bash
#SBATCH --mail-user=conor.waldock@unibe.ch
#SBATCH --mail-type=end,fail,begin
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --tmp=64G
#SBATCH --array=1-9
#SBATCH --output='error_shapley/slurm-%j.out'
#SBATCH --error='output_shapley/slurm-%j.out'
#SBATCH --chdir='/storage/homefs/cw21p621'

array_id=${SLURM_ARRAY_TASK_ID}

# Put your code below this line
module load R
Rscript scripts/shapley/shapRun_main.R $array_id