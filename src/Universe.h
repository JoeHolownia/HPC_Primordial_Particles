//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_UNIVERSE_H
#define PROJECT_UNIVERSE_H
#include <algorithm>
#include <random>
#include <math.h> 
#include <time.h>
#include <cmath>
#include "Particle.h"

class Universe {
public:

    // constructor and core functions
    Universe(int num_particles, int width, int height, double density, double close_radius, double a, double b, double g);
    void InitState();
    void Step();
    void Clean();

    // getters and setters
    Particle* GetCurrentState();
    int GetNumParticles();

private:
    int u_num_particles;  // number of particles
    int u_width;
    int u_height;
    double u_density;  // radius
    double u_close_radius;  // smaller radius for special colouring
    double u_a;  // alpha
    double u_b;  // beta
    double u_g;  // gamma = fixed velocity / radius
    double u_radius;  // radius
    double u_radius_sqrd;
    double u_velocity;  // velocity
    Particle *u_state;  // universe state, i.e. array of all particles
    std::mt19937 u_rand_gen;  // seed for random number generation
};


#endif //PROJECT_UNIVERSE_H
