#!/bin/bash -l
#SBATCH --job-name="benchmark_job"
#SBATCH --output="bench_iters/bench_iters%A_%a.out"
#SBATCH --error="bench_iters/bench_%A_%a.err"
#SBATCH --partition=parallel
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --account=mjohn218_condo
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8gb
#SBATCH --mail-type=END
#SBATCH --mail-user=stakesh1@jh.edu
#SBATCH --array=1-4  # Total combinations of nmer sizes and integrators

source ~/.bashrc


ITERS_SIZES=(20 40 60 80)

total_iters=${#ITERS_SIZES[@]}


# Calculate zero-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))




iters_size=${ITERS_SIZES[$idx]}


# Run the Julia script with the computed arguments
#julia +1.10.0 ./bench_iters.jl 3 "$iters_size" .01 f "KenCarp4()"
julia +1.10.0 ./bench_iters.jl 3 "$iters_size" .01 r "QNDF()"