#!/bin/bash -l
#

#SBATCH --partition=cosc
#SBATCH --job-name=ass2_mpi
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=0-10:00
#SBATCH --mem-per-cpu=2G

module load gnu
module load mpi/openmpi3_eth

# fixed problem size measurements
# for i in 2 4 6 8 10 12 14 16 18 20 22 24 26 28
# do
#   echo $i
#   mpiexec -n $i ./Assignment_mpi 10080 >> mpi_run_out.txt 
# done

# just run once, for testing
mpiexec -n 4 ./particles >> mpi_output.txt