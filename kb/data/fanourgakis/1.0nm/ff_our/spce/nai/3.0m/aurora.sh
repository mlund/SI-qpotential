#!/bin/bash
#SBATCH -p gpu
#SBATCH --exclusive
#SBATCH --gres=gpu:2
#SBATCH --mem-per-cpu=3100
#SBATCH -N 1
#SBATCH -A lu2017-2-5
#
# job time, change for what your job requires
#SBATCH -t 15:00:00
#
# job name
#SBATCH -J grekiss
#
# filenames stdout and stderr - customise, include %j
#SBATCH -o sim.out
#SBATCH -e sim.err

cd $SLURM_SUBMIT_DIR

#module purge
#module load GCC/5.4.0-2.26
#module load CUDA/8.0.44

module add intelcuda
module unload gcc
module load GCC/4.8.4


python run.py