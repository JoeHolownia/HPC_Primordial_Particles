#include <cstring>
#include "Universe.h"
#include "IOParser.h"

void SystemCommandCall(std::string command)
{
  int commLen = command.length();
  char *commandchar = new char[commLen + 1];          // declaring character array
  commandchar = strcpy(commandchar, command.c_str()); // copying the contents of the string to char array
  system(commandchar);                                // Creates the directory incase it does not exist
  delete[] commandchar;
}

int main(int argc, char *argv[]) {

    // instantiate settings from file
    IOParser io_parser("universe_settings.json", "out_display.bin", "out_log.bin");

    // instantiate universe (currently just set to replicate Nature paper)
    Universe universe(5000, 250, 250, 5, 180, 17, 0.67);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    io_parser.OpenOutFile();
    io_parser.WriteStateToOutFile(universe.GetCurrentState(), universe.GetNumParticles());

    // run simulation for all time steps
    //universe.Step();

    // clean up memory
    io_parser.CloseOutFile();
    universe.Clean();

    //Call Python to plot data 
    std::string command = "python debug/display.py";
    SystemCommandCall(command);

    return 0;
}
