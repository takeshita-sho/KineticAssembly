#!/bin/bash -l
#SBATCH --job-name="optim_job"
#SBATCH --output="optim_job_%J.out"
#SBATCH --error="optim_job_%J.err"
#SBATCH --partition=parallel
#SBATCH --time=2-00:00:00
#SBATCH --nodes=2
#SBATCH --account=mjohn218_condo
#SBATCH --ntasks-per-node=2
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=END
#SBATCH --mail-user=stakesh1@jh.edu
##SBATCH --array=1-5  # Total combinations of nmer sizes and integrators

source ~/.bashrc

# Arrays of nmer sizes and integrators
#NMER_SIZES=(3 4 5 6 7)


#total_nmer=${#NMER_SIZES[@]}


# Calculate zero-based index
#idx=$((SLURM_ARRAY_TASK_ID - 1))

# Compute indices for nmer size and integrator
#nmer_index=$(( idx ))


# Retrieve nmer size and integrator
#nmer_size=${NMER_SIZES[$nmer_index]}


# Run the Julia script with the computed arguments
julia +1.10.0 ./efficiency.jl 7 1 .01 f "KenCarp4()"


