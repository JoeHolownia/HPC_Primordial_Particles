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


    // instantiate universe (currently just set to replicate Nature paper)
    Universe universe(5000, 250, 250, 5, 1.3, 180, 17, 0.67);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    io_parser.OpenOutFile();
    io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());

    // Particle* state = universe.GetCurrentState();
    // printf("===ORIGINAL STATE===\n");
    // for (size_t i = 0; i < universe.GetNumParticles(); i++) {
    //     std::cout << "x: " << state[i].x << "  y: " << state[i].y << "  colour: " << state[i].colour << '\n';
    // }

    // run simulation for all time steps
    for (int i = 0; i < 100; i++) {
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
