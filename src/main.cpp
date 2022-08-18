#include <iostream>

#include "Universe.h"

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;

    // instantiate settings from file


    // instantiate universe (currently just set to replicate Nature paper)
    Universe universe(5000, 0.08, 250, 250, 5, 180, 17, 0.67);

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    //std::array currentState = universe.GetCurrentState();
    // IOReader.write(currentState)


    // run simulation for all time steps

    universe.Step();

    // clean up memory
    universe.Clean();


    return 0;
}
