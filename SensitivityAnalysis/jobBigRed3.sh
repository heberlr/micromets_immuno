#!/bin/bash

#SBATCH -J SA_hybrid
#SBATCH -p general
#SBATCH -o SA_hybrid_%j.txt
#SBATCH -e SA_hybrid_%j.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00

module load python/3.6.11
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpu-bind=sockets python hybrid_call_c.py
