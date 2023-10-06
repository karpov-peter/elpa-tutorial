# ScaPAPACK + ELPA tutorial
This is a tutorial for ScaLAPACK + ELPA for the NOMAD Paphos Summer school 2023

## 0. Preparation

Clone this repository to LUMI (you are also welcome to work with your favorite cluster/supercomputer or a laptop).
```
git clone https://github.com/karpov-peter/elpa-tutorial.git
```

[Optional] If you want you can install ScaLAPACK (e.g. as provided by Intel MKL) and ELPA to your laptop.

https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

https://gitlab.mpcdf.mpg.de/elpa/elpa

Consult ELPA user guide Sec. 3.4 or me for the installation of ELPA.

## 1. GEMM

Navigate to the directory ```examples_Fortran/01_gemm``` (cf ```examples_C/01_gemm```). Inspect the source code. Feel free to modify it as you like. Compile it with `bash compile.sh` and run with `sbatch job.sh`.

Check, how fast is OpenMP parallelization.

## 2. PDGEMM

Inspect the code. Check, how fast is the MPI parallel version with respect to the OpenMP one.

## 3. Eigenproblem with ScaLAPACK

Inspect and run the code. Check how fast it is.

## 4. Eigenproblem with ELPA

Compare ScaLAPACK performance with ELPA.

