#!/bin/bash
#SBATCH -J SED
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o outfile_%J.vasp
#SBATCH -e errfile_%J.vasp

run_lammps_file="in.vels"

source /public3/soft/modules/module.sh
module load mpi/intel/2022.1

export PATH=/public3/home/scg5426/lammps/lammps-23Jun2022/src:${PATH}
mpirun -np 64 lmp_intel_cpu_intelmpi -in ${run_lammps_file}
