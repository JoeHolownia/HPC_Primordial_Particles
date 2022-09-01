//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_IOParser_H
#define PROJECT_IOParser_H
#include <iostream>
#include <fstream>
#include <string>
#include "nlohmann/json.hpp"
#include "Particle.h"
using json = nlohmann::json;

struct Settings {
    int num_particles;
    float width;
    float height;
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
    IOParser(std::string settings_fpath, std::string out_disp_fpath, std::string out_log_fpath);
    json* ReadSettingsFile();
    void OpenOutFile();
    void WriteStateToOutFile(Particle* state, int n);
    void CloseOutFile();

private:
        std::string io_settings_fpath;

        std::string io_out_disp_fpath;
        std::ofstream io_out_disp_file; 

        std::string io_out_log_fpath;
        std::ofstream io_out_log_file; 
};


#endif //PROJECT_IOParser_H
