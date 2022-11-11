#include <cstring>
#include "universe.h"
#include "nlohmann/json.hpp"
#include "particle.h"
#include <mpi.h>
#include <chrono>
using json = nlohmann::json;

void SystemCommandCall(std::string command) {
    /**
     * @brief Function to call other processes with a given command.
     */
  #ifdef _WIN64
    // run this on windows OS
    int commLen = command.length();
    char *commandchar = new char[commLen + 1];                      // declaring character array
    int dumb = strcpy_s(commandchar, commLen + 1, command.c_str()); // copying the contents of the string to char array
    if (dumb != 0)
    {
      std::cout << "error in strcpy_s function\n";
      exit(-1);
    }
    // std::cout<<commandchar<<"\n";
    system(commandchar); // Creates the directory incase it does not exist
    delete[] commandchar;
  #else
    // run this on unix based os
    int commLen = command.length();
    char *commandchar = new char[commLen + 1];          // declaring character array
    commandchar = strcpy(commandchar, command.c_str()); // copying the contents of the string to char array
    system(commandchar);                                // Creates the directory incase it does not exist
    delete[] commandchar;
  #endif
}

// this is a helpful struct to take times of things this in a scope
double Totaltime, timeAvg;
int TimerCounter;
struct Timer {
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<float> duration;
  float DurationMilliSec;
  Timer()//Timer is started with the a scoop
  {
    start=std::chrono::high_resolution_clock::now();

  }
  ~Timer()//Timer is stoped when the scoop is existed
  {
    end=std::chrono::high_resolution_clock::now();
    duration =end-start;
    DurationMilliSec=duration.count()*1e3;

     Totaltime = Totaltime + (duration.count()*1e3);
      TimerCounter++;
    std::cout<<DurationMilliSec<<"\n";
  }
};

void calculate_grid_layout(int num_procs, int *grid_rows, int *grid_cols) {
  /**
  *  @brief Determine the number of grid cells in the cartesian space 
            based on the number of processes.
  */

    int i;
    for (i = sqrt(num_procs); i > 0; i--) {
        if (num_procs % i == 0) {
            *grid_rows = i;
            *grid_cols = num_procs / i;
            return;
        }
    }
}

int main(int argc, char *argv[]) {

	// get the current time, for benchmarking
	auto start_time = std::chrono::high_resolution_clock::now();

	// initalise MPI
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// global universe variables
	int global_num_particles, global_width, global_height, time_steps;
	float density, alpha, beta, gamma;

	// read in json settings
	if (rank == 0) {
		
		// read in json data
		std::ifstream json_settings_file("settings.json");
		json settings = json::parse(json_settings_file);
		json_settings_file.close();
		global_num_particles = settings["num_particles"].get<int>();
		global_width = settings["width"].get<int>();
		global_height = settings["height"].get<int>();
		density = settings["density"].get<float>();
		alpha = settings["alpha"].get<float>();
		beta = settings["beta"].get<float>();
		gamma = settings["gamma"].get<float>();
		time_steps = settings["time_steps"].get<int>();
	}

	// broadcast global universe vars to all processes
	MPI_Bcast(&global_num_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&density, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&beta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gamma, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&time_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("Yo I'm proc: %d and here is my Global Data - Num Particles: %d, Height: %d, Width: %d\n", rank, global_num_particles, global_height, global_width);

    // init MPI grid communicator
    MPI_Comm grid_comm;

	// init grid variables
    int grid_rows, grid_cols, reorder = 1, neighbours[NUM_NEIGHBOURS], dims[2], my_coords[2], periods[2] = {1, 1};

	// create processes in 2-D cartesian grid, which wraps.
	calculate_grid_layout(num_procs, &grid_rows, &grid_cols);
    dims[0] = grid_rows;
    dims[1] = grid_cols;
    MPI_Dims_create(num_procs, 2, dims); 
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm);
    MPI_Cart_get(grid_comm, 2, dims, periods, my_coords);

  	// calculate the local box dimension.
	box_coord_type local_box;
    int local_box_width = global_width / grid_cols;
    int local_box_height = global_height / grid_rows;
    local_box.x0 = my_coords[1] * local_box_width;
    local_box.y0 = my_coords[0] * local_box_height;
    local_box.x1 = (my_coords[1] == grid_cols - 1) ? global_width : (my_coords[1] + 1) * local_box_width;
    local_box.y1 = (my_coords[0] == grid_rows - 1) ? global_height : (my_coords[0] + 1) * local_box_height;

    // get the rank of neighbor processors.
    MPI_Cart_shift(grid_comm, 0, 1, &neighbours[DOWN], &neighbours[UP]);
    MPI_Cart_shift(grid_comm, 1, 1, &neighbours[LEFT], &neighbours[RIGHT]);

	// assign particles to each grid cell
	int local_num_particles;
	if (rank == num_procs - 1) {
		// assign the left over particles to the final process (meaning the num particles in the last proc >= all others)
    	local_num_particles = global_num_particles - (rank * (global_num_particles / num_procs));
  	} else {
      	local_num_particles = global_num_particles / num_procs;
  	}

	// instantiate universe (contains properties to be used by each miniverse)
  	Universe* universe = new Universe(global_num_particles, global_width, global_height, density, alpha, beta, gamma);

  	// instantiate miniverse, run by individual process
  	Miniverse miniverse(num_procs, rank, local_box, local_box_width, local_box_height, 
						grid_comm, neighbours, local_num_particles, universe);

	// time recorded for initialisation
	auto finish_init_time = std::chrono::high_resolution_clock::now();

    // run simulation for all time steps
    for (int i = 0; i < time_steps; i++) {

      // need to run miniverse internal step, then sharing step, then position update step
      miniverse.Step();
      miniverse.WriteLocalParticlesToOutFile();
      printf("STEP %d finished... \n", i);
    }

	// parallel portion (simulation) finish time
	auto parallel_finish_time = std::chrono::high_resolution_clock::now();

	// clean particles
	miniverse.Clean();
  	delete universe;

	// record final finish time
	auto finish_time = std::chrono::high_resolution_clock::now();

	if (rank == 0) {

		// print time details to stdout
		auto time_spent_in_init = std::chrono::duration_cast<std::chrono::microseconds>(finish_init_time - start_time);
		auto time_spent_in_parallel = std::chrono::duration_cast<std::chrono::microseconds>(parallel_finish_time - finish_init_time);
		auto total_time_spent = std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time);
		auto average_time_per_step = time_spent_in_parallel.count() / (float) time_steps;
		auto serial_time_spent = total_time_spent.count() - time_spent_in_parallel.count();
		
		std::cout << "Time spent in initialization:                     " << std::setw(12) << time_spent_in_init.count() << " us\n";
      	std::cout << "Time spent in simulation (parallel):              " << std::setw(12) << time_spent_in_parallel.count() << " us\n";
      	std::cout << "Steps:                                            " << std::setw(12) << time_steps << "\n";
      	std::cout << "Average Time Per Step:                            " << std::setw(12) << average_time_per_step << " us\n";
      	std::cout << "Total time:                                       " << std::setw(12) << total_time_spent.count() << " us\n";
      	std::cout << "Serial time:                                      " << std::setw(12) << serial_time_spent << " us\n";

    // call Python to plot data 
		// std::string command = "python3 display.py";
		// SystemCommandCall(command);
    }

	printf("I HAVE MADE IT TO THE END OF MY JOURNEY!!!\n");

	// end MPI
  MPI_Finalize();
  return 0;
}
