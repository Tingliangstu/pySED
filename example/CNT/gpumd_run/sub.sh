#!/bin/sh
gpumd=/fs1/home/ouyang/software/gpumd/gpumd_v3.9.3/src/gpumd
nep=/fs1/home/ouyang/software/gpumd/gpumd_v3.9.3/src/nep

module add CUDA
yhrun -N 1 -p gpu --gpus-per-node=1 --cpus-per-gpu=1 $gpumd > gpumd.out
