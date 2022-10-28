//
// Created by joeho on 13/08/2022.
//

#include "universe.h"
#include "particle.h"
#include <mpi.h>

#define DEGREES_TO_RADIANS 0.0174532925199432957692369076848861271344287188854172545609719144017
#define PI 3.141592653589793115997963468544185161590576171875
#define TAU 2 * PI
#define CLOSE_RADIUS_RATIO 3.846153846


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


Universe::Universe(int num_particles, int width, int height, float density, float a, float b, float g) {
    /**
     * @brief Universe constructor, describes the properties of the overall system. 
     */

    // initialize universe parameters
    u_num_particles = num_particles;
    u_width = width;
    u_height = height;
    u_density = density;
    u_a = a * DEGREES_TO_RADIANS;  // convert to radians
    u_b = b * DEGREES_TO_RADIANS;  // convert to radians
    u_g = g;

    // set radii and velocity based on box size
    u_radius = sqrt((u_width * u_height * u_density) / (u_num_particles * PI));
    u_velocity = u_g * u_radius; 
    u_radius_sqrd = u_radius * u_radius;
    u_close_radius = u_radius / CLOSE_RADIUS_RATIO;
    u_close_radius_sqrd = u_close_radius * u_close_radius;
}


Miniverse::Miniverse(int num_proc, int rank, box_coord_type box_coords, int width, int height, MPI_Comm grid_comm, int num_particles, Universe universe) {
    /**
     * @brief Miniverse constructor, a local square grid cell of the universe managing its own set of particles independently. 
     */

    // initialize miniverse parameters
    m_num_proc = num_proc;
    m_box = box_coords;
    m_width = width;
    m_height = height;
    m_grid_comm = grid_comm;
    m_num_particles = num_particles;
    m_universe = universe;

    // unpack universe settings
    

    // initialise the start state using the given parameters
    InitState();
}

void::Miniverse::InitState() {
    /**
     * @brief Place all particles randomly in the local miniverse box.
     */

    // get seeded uniform random distribution
    // m_universe.u_rand_gen.seed((unsigned int)time(0));
    m_universe.u_rand_gen.seed(100); // FOR TIMING!
    std::uniform_real_distribution<float> uniform_rand(0.0f, 1.0f);

    // initialise all particles with random positions and headings
    for (int i = 0; i < m_num_particles; i++) {
        Particle p;
        p.colour = green;

        // TODO: ENSURE THIS IS WITHIN THE BOUNDS OF OUR BOX!!!
        p.x = m_box.x0 + uniform_rand(m_universe.u_rand_gen) * m_width;
        p.y = m_box.y0 + uniform_rand(m_universe.u_rand_gen) * m_height;

        // InitParticleList(&particlesList);
        // int i, j, k, dir;
        // bool bFlag, hasBigPart;

        // // Initiate particles.
        // for(i = 0; i < lNumParts; i++)
        // {
        //     float x = FloatRand(localBox.x0, localBox.x1);
        //     float y = FloatRand(localBox.y0, localBox.y1);
        //     float r = FloatRand(0, max_vel);
        //     float theta = FloatRand(0, 2 * PI);
        //     float vx = r * cos(theta);
        //     float vy = r * sin(theta);
        //     InsertPartListFront(&particlesList, CreateParticleListItem(x, y, vx, vy));
        //     lTemp += (double)(r * r) / 2;   // Here r = v = velocity.
        // }
        p.heading = uniform_rand(m_universe.u_rand_gen) * TAU;

        // add particle to linked list
        m_particle_list.push_front(p);
    }
}

void Miniverse::Step() {
    /**
     * @brief Perform a time step update of the Universe simulation.
     * 
     */

    // O(n^2) pairwise calculation
    for (int i = 0; i < m_num_particles; i++) {

        // get particle
        Particle &p1 = u_state[i];
        Particle p1_new;

        // counts of particles within left and right semi-circle
        int l = 0;
        int r = 0;
        int n_close = 0;
        
        // interactions
        for (int j = 0; j < u_num_particles; j++) {

            // exclude self from check
            if (i == j) {
                continue;
            }

            // other particle
            Particle &p2 = u_state[j];

            // get deltas
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;

            // wrap deltas, as mirror image about box mid-point
            if (dx > u_width * 0.5f) {
                dx -= u_width;
            } else if (dx < -u_width * 0.5f) {
                dx += u_width;
            }
            if (dy > u_height * 0.5f) {
                dy -= u_height;
            } else if (dy < -u_height * 0.5f) {
                dy += u_height;
            }

            // first check if particle in square (simpler calculation)
            if (abs(dx) <= u_radius && abs(dy) <= u_radius) {

                // check if particle in u_radius circle
                float lhs = dx * dx + dy * dy;
                if (lhs <= u_radius_sqrd) {

                    // check if point is to the left or right of heading to determine the points in each half
                    if (dx * sin(p1.heading) - dy * cos(p1.heading) < 0) {
                        r++;
                    }  else  {
                        l++;
                    }

                    // also check if particle in smaller radius circle, just for colouring
                    if (lhs < u_close_radius_sqrd) {
                        n_close++;
                    }
                }
            }
        }

        // total num within circle
        int n = l + r;

        // set colour
        p1.colour = get_colour(n, n_close);

        // set change in heading direction
        p1.heading += u_a + u_b * n * sign(r - l);

        // apply force to get new x and y
        p1.x += cos(p1.heading) * u_velocity;
        p1.y += sin(p1.heading) * u_velocity;

        // wrap x
        if (p1.x < 0) {
            p1.x += u_width;
        } else if (p1.x >= u_width) {
            p1.x -= u_width;
        }

        // wrap y
        if (p1.y < 0) {
            p1.y += u_height;
        } else if (p1.y >= u_height) {
            p1.y -= u_height;
        }
    }
}

void Miniverse::Clean() {
    /**
     * @brief Clear memory for the particle list.
     * 
     */

    m_particle_list.clear; 
}

std::list<particle_type> Miniverse::GetParticleList() {
    /**
     * @brief Returns the current state of the universe.
     * 
     */
    return m_particle_list;
}

int Miniverse::GetNumParticles() {
    /**
     * @brief Returns the number of particles.
     */
    return m_num_particles;
}

