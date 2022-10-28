#include <cstring>
#include "universe.h"
#include "ioparser.h"
#include "nlohmann/json.hpp"
#include "particle.h"
#include <mpi.h>
using json = nlohmann::json;

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

void CalculateGridLayout(int nProcs, int *gridRows, int *gridCols); // PLEASE MOVE THIS TO A SEPERATE HEADER FILE, MAYBE MPI HELPERS??

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

void calculate_grid_layout(int num_procs, int *grid_rows, int *grid_cols)
{
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

	// initalise MPI
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// create MPI particle type // TODO: may need to compact this to just coords? but for convenience encode the whole particle.
	MPI_Datatype particle_types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Datatype mpi_particle_type;
    int block_counts[4] = {4, 4, 4, 1};
    MPI_Aint offsets[4], extent;
    offsets[0] = 0;
    MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[1] = 4 * extent;
	MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[2] = 4 * extent;
	MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[3] = 4 * extent;
    MPI_Type_struct(4, block_counts, offsets, particle_types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);

	// global universe variables
	int global_num_particles, global_width, global_height, time_steps;
	float density, alpha, beta, gamma;

	// instantiate io parser
	IOParser io_parser;
    if (rank == 0) {
		io_parser("settings.json", "out_display.bin", "out_log.bin");

		// read in json data: TODO: THIS IS PROBLEMATIC ATM, EITHER NEED TO STORE THESE GLOBALLY OR READ IN AND THEN DISTRIBUTE TO EACH PROCESS
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

    // init MPI grid communicator
    MPI_Comm grid_comm;
    MPI_Request request[4];
    MPI_Status status[4];  // I BELIEVE 4 IS BECAUSE OF THE left, top, right, down --> if we do diagonal this needs to be 8

	// init grid variables
    int grid_rows, grid_cols, reorder = 1, neighbours[4], send_counts[4]={0}, dims[2], my_coords[2], periods[2] = {1, 1};

	// create a division of processes in 2-D cartesian grid.
	calculate_grid_layout(num_procs, &grid_rows, &grid_cols);
    dims[0] = grid_rows;
    dims[1] = grid_cols;
    MPI_Dims_create(num_procs, 2, dims); 

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm);
    MPI_Cart_get(grid_comm, 2, dims, periods, my_coords); // myCoords[0] = row, myCoords[0] = col.

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

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    if (rank == 0) {
		io_parser.OpenOutFile();
		io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());
    }

	// instantiate universe (contains properties to be used by each miniverse)
  	Universe universe(global_num_particles, global_width, global_height, density, alpha, beta, gamma);

  	// instantiate miniverse, run by individual process
  	Miniverse miniverse(num_procs, rank, local_box, local_box_width, local_box_height, grid_comm, local_num_particles, universe);

	printf("Hello from process %d out of %d\n", rank, num_procs);

    // run simulation for all time steps
    for (int i = 0; i < time_steps; i++) {

		// need to run miniverse internal step, then sharing step, then position update step

		/*
		//Send neighbor's data to the corresponding neighbor.
        for(j = 0; j < 4; j++){
          MPI_Isend(sendBuffer[j], sendCounts[j], mpiPartType, neighbours[j], 0, gridComm, &(request[j]));
        }
        //Receive own data from neighbor.
        for(j = 0; j < 4; j++){
            MPI_Recv(recvBuffer[j], COMM_BUFFER_SIZE, mpiPartType, MPI_ANY_SOURCE, 0, gridComm, &(status[j]));
        }
        MPI_Waitall(4, request, MPI_STATUS_IGNORE); //Wait for non-blocking send completion.
		*/

		// Timer text_time;
		// universe.Step();

		// write output state to file
		// io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());
    }

    // output average time
    // timeAvg = Totaltime / (float) time_steps;
    // std::cout<<timeAvg<<"ms\n";

	// clean particles
	miniverse.Clean();

    if (rank == 0) {
		// close out file
		io_parser.CloseOutFile();

    	// call Python to plot data 
		// std::string command = "python3 display.py";
		// SystemCommandCall(command);
    }

    return 0;
}
