# HPC Primordial Particle System
This project contains a serial, OpenMP and MPI C++ implementation of the Primordial Particle system described in the article published in Nature: "How life emerges from a simplistic particle motion law".

This is the repository for the COSC3500: High Performance Computing Project.

Article: https://www.nature.com/articles/srep37969
YouTube Video: https://www.youtube.com/watch?v=makaJpLvbow&ab_channel=ArtificialLifeLabGraz

Some of the radius and velocity scaling functionality was used with help from:
https://github.com/nagualdesign/Primordial-Particle-System

The MPI Implementation cartesian coordinate system was built with guidance from:
https://github.com/mokarrom/mpi-part-simulation


The project structure is broken into two source directories:

1) src: for serial / OpenMP implementations
2) mpi_src: for MPI implementation

This is for ease of use, and because of the fundamental changes made in the MPI impelementation.

USE:

For using the serial / OpenMP implementation:

1. Ensure conda Python3 is in the PATH, or you are in a conda Python3 environment with matplotlib and numpy installed.
2. Modify the parameters in the settings.json file
3. On windows, install mingw32 and add mingw32-make to the PATH, and then run 'mingw32-make all' in the src directory.
3. On Linux, simply run 'make all' in the src directory.
4. For OpenMP use, set the following variable in bash: export OMP_NUM_THREADS=n, where n is the number of threads
4. Then, run the particles executable by using  -n <num_proc> ./particles
5. A gif of the run and an image of the start state will be saved to the results directory.

For using the MPI implementation:

First run the following in your terminal: 

- module load gnu
- module load mpi/openmpi3_eth

Then:

1. Run make in the mpi_src directory
2. Modify the parameters in the settings.json file
3. Then run the particles executable by using mpiexec -n <num_proc> ./particles_mpi, where n is the desired number of MPI processes
4. A gif of the run and an image of the start state will be saved to the results directory.

NOTE - FOR RUNNING THE MPI IMPLEMENTATION!!
- The MPI implementation sometimes hangs indefinitely on larger problem sizes and low numbers of processes (due to memory scaling errors), it
is recommended if running a very large simulation (8000+ particles) to use at least 10 MPI processes. Otherwise, the density parameter can
be reduced to 4 if you wish to run large problem sizes with only a few MPI processes.
- If an error is reported in python display.py file regarding array index errors, run make clean in the mpi_src directory (this can happen
if changing the number of processes run if there are more binary files than processes)


References:
1. nagualdesign, “Primordial Particle System”. Github. https://github.com/nagualdesign/Primordial-Particle-System (accessed Sep. 2, 2022).
2. mokarrom, “mpi-part-simulation”. Github. https://github.com/mokarrom/mpi-part-simulation (accessed Oct. 21, 2022)
