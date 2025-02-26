#!/bin/bash
#SBATCH --job-name=run_mt-SCITE
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  # Adjust CPU count as needed
#SBATCH --mem-per-cpu=4G  # Adjust memory as needed
#SBATCH --time=24:00:00 
#SBATCH --mail-type=END,FAIL

# Run Snakemake with SLURM job submission
# --keep-going --rerun-incomplete
snakemake --jobs 50 --latency-wait 60  \
    --cluster "sbatch --partition=compute --nodes=1 --ntasks=1 --cpus-per-task=8 --mem-per-cpu=4G --time=24:00:00"

