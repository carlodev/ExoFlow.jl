#!/bin/sh
export PATH=$HOME/julia-1.8.5/bin/:$PATH

srun julia --project=../ -O3 --check-bounds=no -L run_ExoFlow.jl