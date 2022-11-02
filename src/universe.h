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
#include <cstring>
#include "particle.h"
#include "ioparser.h"
#include <list>
#include <mpi.h>
#include <array>

#define COMM_BUFFER_SIZE 1000
#define NUM_NEIGHBOURS 4
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#define UP_RIGHT 4
#define DOWN_RIGHT 5
#define DOWN_LEFT 6
#define UP_LEFT 7

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
              int neighbours[NUM_NEIGHBOURS], int num_particles, Universe* universe);
    void InitState();
    void Step();
    std::array<int, 2> check_particle_edge_contact(particle_type* p);
    int check_particle_box(particle_type* p);
    void WriteLocalParticlesToOutFile();
    void Clean();

    // getters and setters
    std::list<particle_type*> GetParticleList();
    int GetNumParticles();

private:

    // miniverse properties
    int m_num_proc; // MPI num procs
    int m_rank; // MPI rank
    box_coord_type m_box; // box coords
    int m_width;
    int m_height;
    MPI_Comm m_grid_comm; // MPI grid communicator
    int m_num_particles;
    Universe* m_universe; // properties of overall universe

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

    // io parser
    IOParser* io_parse_obj;

    // grid communication variables and buffers
    int* m_neighbours;
    int m_send_counts[NUM_NEIGHBOURS]={0};
    particle_type send_buffer[NUM_NEIGHBOURS][COMM_BUFFER_SIZE];
    particle_type recv_buffer[NUM_NEIGHBOURS][COMM_BUFFER_SIZE];
};

#endif //PROJECT_UNIVERSE_H
