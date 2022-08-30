//
// Created by joeho on 13/08/2022.
//

#include <math.h> 
#include "Universe.h"
#include "Particle.h"

// DEFINE MACROS HERE!
#define CIRCLE_DEGREES 360
#define RADIANS_TO_DEGREES = 57.2957795131
#define PI 3.14159265

int sign(int x) {
    /**
     * @brief calculates the sign(x) of input x.
     */
    if (x < 0) {
        return -1;
    } else if (x == 0) {
        return 0;
    } else {
        return 1;
    }
}

Colour get_colour(int n, int n_close) {
    /**
     * @brief determines colour based on an input number of neighbours, n.
     * 
     *  The order of checks is optimized to allow the least checks for most particles,
     *  i.e. it is expected that most particles are green at any time step.
     */

    // check common, lower cases
    if (n < 13) {
        return green;

    } else if (13 <= n && n <= 15) {
        return brown;

    } else if (15 < n_close) {
        return magenta;

    } else if (15 < n && n <= 35) {
        return blue;

    } else {
        return yellow;
    }
}


Universe::Universe(int num_particles, int width, int height, int radius, int close_radius, float a, float b, float velocity) {
    /**
     * @brief Universe constructor. 
     */

    // initialize universe parameters
    u_num_particles = num_particles;
    u_width = width;
    u_height = height;
    u_radius = radius;
    u_close_radius = close_radius;
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

    // DUMMY INITSTATE FOR DEBUG
    // int colour = green;
    // for (int i = 0; i < u_num_particles; i++) {
    //     Particle& p = u_state[i];
    //     p.colour = (Colour) colour;
    //     colour++;
    //     p.x = (float) i + 1.0;
    //     p.y = (float) i + 1.0;
    //     p.heading = uniform_rand(u_rand_gen) * CIRCLE_DEGREES;
    // }

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

        // current particle, and corresponding new memory address
        Particle &p1 = u_state[i];
        Particle &new_p1 = new_state[i];

        // counts of particles within left and right semi-circle
        int l = 0;
        int r = 0;
        int n_close = 0;
        
        // interactions
        for (int j = 0; j < u_num_particles; j++) {

            // other particle
            Particle &p2 = u_state[j];

            // get deltas
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;

            // TODO: first check if particle in square (simpler calculation) (OPTIMIZATION)
            if (abs(dx) <= u_radius && abs(dy) <= u_radius) {

                // check if particle in u_radius circle
                float lhs = dx * dx + dy * dy;
                if (lhs < u_radius * u_radius) {

                    // check if point is to the left or right in x to determine the points in each half
                    if (dx <= 0) {l++;} else {r++;}

                    // also check if particle in  smaller radius circle
                    if (lhs < u_close_radius * u_close_radius) {
                        n_close++;
                    }
                }
            }
        }

        // total num within circle
        int n = l + r;

        // set colour
        new_p1.colour = get_colour(n, n_close);

        // set change in heading direction
        float d_phi = u_a + u_b * n * sign(r - l);
        new_p1.heading = p1.heading + d_phi;

        // apply force to get new x and y
        // TODO: OPTIMISATION HERE IS TO FILL A LOOKUP TABLE WITH VALUES!
        new_p1.x = p1.x + cos(new_p1.heading) * u_velocity;
        new_p1.y = p1.y + sin(new_p1.heading) * u_velocity;

        // // DUMMY TEST        
        // new_p1.x = p1.x + 1.0;
        // new_p1.y = p1.y + 1.0;
        // new_p1.colour = (Colour) (((int) p1.colour) + 1);
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

