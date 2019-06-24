#!/bin/bash
#SBATCH --job-name=matlab_p100
#SBATCH --time=12:00:00
#SBATCH -p gpup100
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=6
### Using more tasks because default memory is ~5GB per core
### 'shared' will share the node with other users
### 'parallel' use entire node (24,28,48, depends on node type)
### Try the script with --ntasks-per-node=1 and see what happens
#
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------
 
ml cuda/8.0
ml gcc/5.5.0
ml # confirm modules used

for foldername in ~/work/yixuan/data/{1485..1500} ; do
    for filename in $foldername/*.h5; do
        ./k-Wave/binaries/kspaceFirstOrder3D-CUDA -i "$filename" -o "$filename" -p
    done
done
