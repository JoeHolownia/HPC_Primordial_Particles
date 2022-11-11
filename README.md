# HPC Primordial Particle System

This is the repository for the COSC3500: High Performance Computing Project.

This project is a serial C++ implementation of the Primordial Particle system described in the article published in Nature: "How life emerges from a simplistic particle motion law".

Article: https://www.nature.com/articles/srep37969
YouTube Video: https://www.youtube.com/watch?v=makaJpLvbow&ab_channel=ArtificialLifeLabGraz

Some of the radius and velocity scaling functionality was used with help from:
https://github.com/nagualdesign/Primordial-Particle-System

USE:

First run: 
module load gnu
module load mpi/openmpi3_eth

1. Ensure conda Python3 is in the PATH, or you are in a conda Python3 environment with matplotlib and numpy installed.
2. Modify the parameters in the settings.json file
3. On windows, install mingw32 and add mingw32-make to the PATH, and then run 'mingw32-make all' in the src directory.
3. On Linux, simply run 'make all' in the src directory.
4. Then, run the particles executable by using mpiexec -n <num_proc> ./particles
5. A gif of the run and an image of the start state will be saved to the results directory.