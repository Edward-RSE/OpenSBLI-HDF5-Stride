#!/bin/bash

#SBATCH --partition=scavenger_mathsa100
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --exclusive

source $HOME/modules/opensbli.lmod

nvidia-smi

# export CUDA_MPS_PIPE_DIRECTORY=$HOME/mps
# export CUDA_MPS_LOG_DIRECTORY=$HOME/mps
# nvidia-cuda-mps-control -d

time mpirun -n 2 ./build/slice-example_mpi_cuda

ls -lth
