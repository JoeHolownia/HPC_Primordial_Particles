#include <cstring>
#include "Universe.h"
#include "IOParser.h"
#include "nlohmann/json.hpp"
using json = nlohmann::json;

void SystemCommandCall(std::string command)
{
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

    // instantiate universe
    Universe universe(num_particles, width, height, density, alpha, beta, gamma);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    io_parser.OpenOutFile();
    io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());

    // run simulation for all time steps
    for (int i = 0; i < time_steps; i++) {
      // Timer text_time;
      universe.Step();
      io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());
    }

    // output average time
    // timeAvg = Totaltime / (float) time_steps;
    // std::cout<<timeAvg<<"ms\n";

    // clean up memory
    io_parser.CloseOutFile();
    universe.Clean();

    // call Python to plot data 
    std::string command = "python display.py";
    SystemCommandCall(command);

    return 0;
}
