#!/bin/bash

#SBATCH -p gpu
#SBATCH --exclusive
#SBATCH --gres=gpu:2
#SBATCH --mem-per-cpu=3100
#SBATCH -N 1
#SBATCH -A lu2019-2-15

# job time, change for what your job requires
#SBATCH -t 21:00:00

# job name
#SBATCH -J qpot-m3

# filenames stdout and stderr - customise, include %j
#SBATCH -o sim.out
#SBATCH -e sim.err

cd $SLURM_SUBMIT_DIR

#module purge
#module load GCC/8.3.0  
module add intelcuda
module unload gcc
module load GCC/4.8.4
module load CUDA/10.1.243

python run.py
