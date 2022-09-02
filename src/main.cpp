#include <cstring>
#include "Universe.h"
#include "IOParser.h"

void SystemCommandCall(std::string command)
{
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


  // CHECK OS TO RUN THIS ON LINUX:
  // int commLen = command.length();
  // char *commandchar = new char[commLen + 1];          // declaring character array
  // commandchar = strcpy(commandchar, command.c_str()); // copying the contents of the string to char array
  // system(commandchar);                                // Creates the directory incase it does not exist
  // delete[] commandchar;
}


int main(int argc, char *argv[]) {

    // instantiate io parser
    IOParser io_parser("universe_settings.json", "out_display.bin", "out_log.bin");

    // instantiate settings from file
    // json settings = *io_parser.ReadSettingsFile();
    // int num_particles = settings["num_particles"];
    // float width = settings["width"];
    // float height = settings["height"];
    // float radius = settings["radius"];
    // float close_radius = settings["close_radius"];
    // float alpha = settings["alpha"];
    // float beta = settings["beta"];
    // float velocity = settings["velocity"];
    //int time_steps = settings["time_steps"];
    int time_steps = 2000;

    // instantiate universe
    //Universe universe(num_particles, width, height, radius, close_radius, alpha, beta, velocity); TODO: this!!
    Universe universe(500, 50, 50, 5.0f, 1.3f, 180.0f, 17.0f, 0.67f);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    io_parser.OpenOutFile();
    io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());

    // Particle* state = universe.GetCurrentState();
    // printf("===ORIGINAL STATE===\n");
    // for (size_t i = 0; i < universe.GetNumParticles(); i++) {
    //     std::cout << "x: " << state[i].x << "  y: " << state[i].y << "  colour: " << state[i].colour << '\n';
    // }

    // run simulation for all time steps
    for (int i = 0; i < time_steps; i++) {
      universe.Step();
      io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());
    }

    // state = universe.GetCurrentState();
    // std::cout << "===STATE AFTER STEP===\n";
    // for (size_t i = 0; i < universe.GetNumParticles(); i++) {
    //     std::cout << "x: " << state[i].x << "  y: " << state[i].y << "  colour: " << state[i].colour << '\n';
    // }
    // std::cout.flush();

    // clean up memory
    io_parser.CloseOutFile();
    universe.Clean();

    //Call Python to plot data 
    std::string command = "python display.py";
    SystemCommandCall(command);

    return 0;
}
