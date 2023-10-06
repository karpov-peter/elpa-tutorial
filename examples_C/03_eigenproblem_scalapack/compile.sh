#!/usr/bin/bash

set -v

module purge

module load PrgEnv-cray #PrgEnv-gnu
module load cray-libsci

FileName="eigenproblem_scalapack"

cc $FileName.c -o $FileName
