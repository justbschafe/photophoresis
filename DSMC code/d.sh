#!/bin/bash
#SBATCH -J d1b05
#SBATCH -o d.lst
#SBATCH -e d.err
#SBATCH -n 60
#SBATCH -p shared
#SBATCH -t 0-12:00:00 
#SBATCH --mem-per-cpu=4000

ulimit -l unlimited 
mpirun -np 60 d.exe








