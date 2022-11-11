//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_UNIVERSE_H
#define PROJECT_UNIVERSE_H
#include <algorithm>
#include <random>
#include <math.h> 
#include <cmath>
#include <chrono>
#include <iostream>
#include "particle.h"

class Universe {
public:

    // constructor and core functions
    Universe(int num_particles, int width, int height, float density, float a, float b, float g, unsigned int seed);
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
    float u_density;  // radius
    float u_close_radius;  // smaller radius for special colouring
    float u_a;  // alpha
    float u_b;  // beta
    float u_g;  // gamma = fixed velocity / radius
    float u_radius;  // radius
    float u_radius_sqrd;
    float u_close_radius_sqrd;
    float u_velocity;  // velocity
    Particle *u_state;  // universe state, i.e. array of all particles
    std::mt19937 u_rand_gen;  // random number generation
    unsigned int u_seed;
};


#endif //PROJECT_UNIVERSE_H
