#include <cstring>
#include "universe.h"
#include "ioparser.h"
#include "nlohmann/json.hpp"
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

int main(int argc, char *argv[]) {

  	// get the current time, for benchmarking
	auto start_time = std::chrono::high_resolution_clock::now();

    // instantiate io parser
    IOParser io_parser("settings.json", "out_display.bin", "out_log.bin");

    // read in json data
    std::ifstream json_settings_file("settings.json");
    json settings = json::parse(json_settings_file);
    json_settings_file.close();
    int num_particles = settings["num_particles"].get<int>();
    int width = settings["width"].get<int>();
    int height = settings["height"].get<int>();
    float density = settings["density"].get<float>();
    float alpha = settings["alpha"].get<float>();
    float beta = settings["beta"].get<float>();
    float gamma = settings["gamma"].get<float>();
    int time_steps = settings["time_steps"].get<int>();

    // define seed
    //unsigned int seed = (unsigned int)time(0);
    unsigned int seed = 100; // for runtime analyses

    // instantiate universe
    Universe universe(num_particles, width, height, density, alpha, beta, gamma, seed);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    io_parser.OpenOutFile();
    io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());

  	// time recorded for initialisation
	auto finish_init_time = std::chrono::high_resolution_clock::now();

    // run simulation for all time steps
    auto parallel_time = std::chrono::microseconds::zero();
    for (int i = 0; i < time_steps; i++) {


		auto step_start_time = std::chrono::high_resolution_clock::now();
		
		// run PPS algorithm update step
		universe.Step();
		auto step_finish_time = std::chrono::high_resolution_clock::now();
		parallel_time += std::chrono::duration_cast<std::chrono::microseconds>(step_finish_time - step_start_time);
		
		// write state output to binary file
		io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());
    }

    // clean up memory
    io_parser.CloseOutFile();
    universe.Clean();

    // record final finish time
    auto finish_time = std::chrono::high_resolution_clock::now();

    // get time information
	auto time_spent_in_init = std::chrono::duration_cast<std::chrono::microseconds>(finish_init_time - start_time);
	auto total_time_spent = std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time);
	auto average_time_per_step = parallel_time.count() / (float) time_steps;
	auto serial_time_spent = total_time_spent.count() - parallel_time.count();

    // // for recording data [Num Particles, Num Procs, Total Time, Total Serial Time, Time Steps, Time Per Step, Total Parallel Time]
	std::cout << num_particles << "," << total_time_spent.count() << "," << serial_time_spent << "," 
	<< time_steps << "," << average_time_per_step << "," << parallel_time.count() << "\n";

  	// call Python to plot data 
    std::string command = "python3 display.py";
    SystemCommandCall(command);

    return 0;
}
