#include <iostream>

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;

    // instantiate settings from file


    // instantiate universe
    Universe Universe();

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
    std::array currentState = universe.GetCurrentState();
    // IOReader.write(currentState)


    // run simulation for all time steps

    universe.Step();

    //


    return 0;
}
