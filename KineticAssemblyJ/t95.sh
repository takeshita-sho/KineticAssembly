#!/bin/bash -l
#SBATCH --job-name="optim_job"
#SBATCH --output="t95_%A_%a.out"
#SBATCH --error="t95_%A_%a.err"
#SBATCH --partition=parallel
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --account=mjohn218_condo
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=END
#SBATCH --mail-user=stakesh1@jh.edu
#SBATCH --array=1-5  # Total combinations of nmer sizes and integrators

source ~/.bashrc

# Arrays of nmer sizes and integrators
NMER_SIZES=(3 4 5 6 7)
#INTEGRATORS=("Rosenbrock23()" "TRBDF2()" "QNDF()" "Rodas4P()" "Rodas5P()" "Kvaerno5()" "KenCarp4()")

total_nmer=${#NMER_SIZES[@]}
#total_integrators=${#INTEGRATORS[@]}

# Calculate zero-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

# Compute indices for nmer size and integrator
nmer_index=$(( idx ))
#integrator_index=$(( idx % total_integrators ))

# Retrieve nmer size and integrator
nmer_size=${NMER_SIZES[$nmer_index]}
#integrator=${INTEGRATORS[$integrator_index]}

# Run the Julia script with the computed arguments
julia +1.10.0 ./t95.jl "$nmer_size" 1 .01 f "KenCarp4()"