//
// Created by joeho on 13/08/2022.
//

#include "Universe.h"
#include "Particle.h"

// DEFINE MACROS HERE!
#define CIRCLE_DEGREES 360

Universe::Universe(int num_particles, int width, int height, int radius, float a, float b, float velocity) {
    /*
     * Instantiate Universe.
     */

    // initialize universe parameters
    u_num_particles = num_particles;
    u_width = width;
    u_height = height;
    u_radius = radius;
    u_a = a;
    u_b = b;
    u_velocity = velocity;

    // initialise the start state using the given parameters
    InitState();
}

void::Universe::InitState() {
    /**
     * @brief Place all particles randomly in the Universe box using the given density and box size.
     * 
     */

    // instantiate state with particles
    u_state = new Particle[u_num_particles];

    // get seeded uniform random distribution
    u_rand_gen.seed(10);  // (unsigned int)time(0) --> can use this to be different each time
    std::uniform_real_distribution<float> uniform_rand(0.0f, 1.0f);

    // initialise all particles with random positions and headings
    for (int i = 0; i < u_num_particles; i++) {

        // NOTE* WE EVENTUALLY PROBABLY WANT TO CONSIDER IN WHAT ORDER WE ASSIGN THESE RANDOMLY (WHEN DOING MPI)
        // (THIS CAN SUPPORT CONTIGUOUS MEMORY ACCESS)
        Particle& p = u_state[i];
        p.colour = green;
        p.x = uniform_rand(u_rand_gen) * u_width;
        p.y = uniform_rand(u_rand_gen) * u_height;
        p.heading = uniform_rand(u_rand_gen) * CIRCLE_DEGREES;
    }
}

void Universe::Step() {
    /**
     * @brief Perform a time step update of the Universe simulation.
     * 
     */

    // create new state to assign to
    Particle* new_state = new Particle[u_num_particles];

    // naive O(n^2) pairwise calculation
    for (int i = 0; i < u_num_particles; i++) {

        // current particle
        Particle &p1 = u_state[i];

        // interactions
        // for (int j = 0; j < u_num_particles; j++) {

        //     // other particle
        //     const Particle &p2 = u_state[j];

        //     // get deltas
        //     float dx = p1.x - p2.x;
        //     float dy = p1.y - p2.y;

        //     // first check if particle in square (simpler calculation)

               // now check if particle in circle
               // (x - center_x)² + (y - center_y)² < radius²

               // then simply check if point is above or below in y to determine the points in each half

        // }

        Particle &new_p1 = new_state[i];
        new_p1.x += 1;
        new_p1.y += 1;
        new_p1.colour = blue;
    }


    // delete memory from old state u_state, and then set pointer to new memory state
    delete[] u_state;
    u_state = new_state;
}

void Universe::Clean() {
    /**
     * @brief Clear memory for all particles, i.e. clear the state.
     * 
     */

    delete[] u_state; 
}

Particle* Universe::GetCurrentState() {
    /**
     * @brief Returns the current state of the universe.
     * 
     */
    return u_state;
}

int Universe::GetNumParticles() {
    return u_num_particles;
}

