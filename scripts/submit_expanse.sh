#!/bin/bash
#SBATCH --job-name="domain_generation"
#SBATCH --output="domain_generation.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=125
#SBATCH --mem=0 
#SBATCH --account=nca125
#SBATCH --export=ALL
#SBATCH -t 12:00:00

#This job runs with 1 nodes, 128 cores per node for a total of 128 tasks.

module purge
module load cpu
#Load module file(s) into the shell environment
module load slurm
module load cpu/0.17.3b  gcc/10.2.0/npcyll4
module load openmpi/4.1.1

cd /home/tweigand/tcat-1p-exact/sim-steady-compressible
./generate_domain.sh
./run_inlet_velocities.sh
