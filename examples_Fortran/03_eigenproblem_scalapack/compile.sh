#!/usr/bin/bash

#set -v

module purge

module load PrgEnv-gnu
module load cray-libsci

FileName="eigenproblem_scalapack"

#ftn $FileName.f90 -o $FileName
ftn -O0 -g $FileName.f90 -o $FileName
