//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_UNIVERSE_H
#define PROJECT_UNIVERSE_H
#include <random>
#include "Particle.h"

class Universe {
public:

    // constructor and core functions
    Universe(int num_particles, int width, int height, int radius, int close_radius, float a, float b, float velocity);
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
    int u_radius;  // radius
    int u_close_radius;  // smaller radius for special colouring
    float u_a;  // alpha
    float u_b;  // beta
    float u_velocity;  // fixed-velocity
    Particle *u_state;  // universe state, i.e. array of all particles
    std::mt19937 u_rand_gen;  // seed for random number generation
};


#endif //PROJECT_UNIVERSE_H
