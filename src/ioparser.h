//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_IOParser_H
#define PROJECT_IOParser_H
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include "particle.h"

void print_array(const double* A, int size);

struct Settings {
    int num_particles;
    int width;
    int height;
    float radius;
    float close_radius;
    float alpha;
    float beta;
    float velocity;
    int time_steps;
};

class IOParser {
public:

    // constructor and core functions
    IOParser(std::string out_disp_fpath, std::string out_log_fpath);
    void OpenOutFile();
    void WriteStateToOutFile(std::list<particle_type*> state, int n);
    void CloseOutFile();

private:
        std::string io_out_disp_fpath;
        std::ofstream io_out_disp_file;
        std::string io_out_log_fpath;
        std::ofstream io_out_log_file; 
};


#endif //PROJECT_IOParser_H
