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
#include <chrono>
#include <iostream>
#include "particle.h"
#include <list>
#include <mpi.h>

#define COMM_BUFFER_SIZE 2000

class Universe {
public:

    // constructor and core functions
    Universe(int num_particles, int width, int height, float density, float a, float b, float g);

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
    std::mt19937 u_rand_gen;  // seed for random number generation
};

typedef class Universe universe_type;

class Miniverse {
public:

    // constructor and core functions
    Miniverse(int num_proc, int rank, box_coord_type box_coords, int width, int height, MPI_Comm grid_comm, 
              int num_particles, Universe* universe);
    void InitState();
    void Step();
    void Clean();
    std::list<particle_type*> GetParticleList();

    // getters and setters
    particle_type* GetCurrentState();
    int GetNumParticles();

private:
    int m_num_proc; // MPI num procs
    int m_rank; // MPI rank
    box_coord_type m_box; // box coords
    int m_width;
    int m_height;
    MPI_Comm m_grid_comm; // MPI grid communicator
    int m_num_particles;
    Universe* m_universe; // properties of overall universe

    int m_send_counts[4]={0, 0, 0, 0};

    // version of universe properties for convenience
    int u_num_particles;
    int u_width;
    int u_height;
    float density;
    float alpha;
    float beta;
    float gamma;
    float radius;
    float velocity;
    float radius_sqrd;
    float close_radius;
    float close_radius_sqrd;

    // particles for miniverse as a doubly linked list of pointers to particles
    std::list<particle_type*> m_particle_list;

    // instantiate grid communication buffers
    particle_type send_buffer[4][COMM_BUFFER_SIZE];
    particle_type recv_buffer[4][COMM_BUFFER_SIZE];
};

#endif //PROJECT_UNIVERSE_H
