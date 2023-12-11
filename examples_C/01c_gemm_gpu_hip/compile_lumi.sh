#!/bin/bash

module load PrgEnv-gnu
module load craype-accel-amd-gfx90a rocm

FileName=gemm_big_matrices_hip
cc -fopenmp -o $FileName $FileName.c -D__HIP_PLATFORM_AMD__ -lhipblas
