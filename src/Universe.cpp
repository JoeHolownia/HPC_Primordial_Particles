//
// Created by joeho on 13/08/2022.
//

#include "Universe.h"
#include "Particle.h"

// DEFINE MACROS HERE!
#define CIRCLE_DEGREES 360
#define DEGREES_TO_RADIANS 0.0174532925199432957692369076848861271344287188854172545609719144017
#define PI 3.141592653589793115997963468544185161590576171875
#define TAU 2 * PI


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

int* permuted_array(int n) {
    /**
     * @brief creates an array filled with a random permutation of values 0...n.
     * 
     */

    int* p_array = new int[n];

    // fill array
    for (int i = 0; i < n; i++) {
        p_array[i] = i;
    }

    // permute array
    std::random_shuffle(p_array, p_array + n);

    return p_array;
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


Universe::Universe(int num_particles, int width, int height, double density, double close_radius, double a, double b, double g) {
    /**
     * @brief Universe constructor. 
     */

    // initialize universe parameters
    u_num_particles = num_particles;
    u_width = width;
    u_height = height;
    u_density = density;
    u_close_radius = close_radius;
    u_a = a * DEGREES_TO_RADIANS;  // convert to radians
    u_b = b * DEGREES_TO_RADIANS;  // convert to radians
    u_g = g;

    // initialise the start state using the given parameters
    InitState();
}

void::Universe::InitState() {
    /**
     * @brief Place all particles randomly in the Universe box, and determine radius and velocity
     * based on density and box dimensions.
     * 
     */

    // instantiate state with particles
    u_state = new Particle[u_num_particles];

    // get seeded uniform random distribution
    u_rand_gen.seed((unsigned int)time(0));
    std::uniform_real_distribution<float> uniform_rand(0.0f, 1.0f);

    // initialise all particles with random positions and headings
    for (int i = 0; i < u_num_particles; i++) {
        Particle& p = u_state[i];
        p.colour = green;
        p.x = uniform_rand(u_rand_gen) * u_width;
        p.y = uniform_rand(u_rand_gen) * u_height;
        p.heading = uniform_rand(u_rand_gen) * CIRCLE_DEGREES;
    }

    // set radius and velocity based on box size
    u_radius = sqrt((u_width * u_height * 1 * u_density) / (u_num_particles * PI));
    u_velocity = u_g * u_radius; 
    u_radius_sqrd = u_radius * u_radius;

    // DUMMY INITSTATE FOR DEBUG
    // double xs[4] = {5.0, 5.0, 4.0, 5.0};
    // double ys[4] = {5.0, 6.0, 5.0, 4.0};
    // Colour colours[4] = {magenta, green, green, green};
    // double headings[4] = {0, 0, 0, 0};
    // for (int i = 0; i < u_num_particles; i++) {
    //     Particle& p = u_state[i];
    //     p.x = xs[i];
    //     p.y = ys[i];
    //     p.colour = colours[i];
    //     p.heading = headings[i];
    // }
}

void Universe::Step() {
    /**
     * @brief Perform a time step update of the Universe simulation.
     * 
     */

    // array of random particle move orders for asynchronous update in this step
    int* move_order = permuted_array(u_num_particles);

    // O(n^2) pairwise calculation
    for (int i = 0; i < u_num_particles; i++) {

        // choose current particle randomly
        // int rand_index = move_order[i];
        //printf("Rand Index P1:  %d\n", rand_index);
        //Particle &p1 = u_state[rand_index];
        //printf("Particle P1:  %d\n", i);
        Particle &p1 = u_state[i];

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
            //printf("Particle P2:  %d\n", j);

            // get deltas
            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;

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
            //if (abs(dx) <= u_radius && abs(dy) <= u_radius) {

            // check if particle in u_radius circle
            double lhs = dx * dx + dy * dy;
            if (lhs <= u_radius_sqrd) {

                // check if point is to the left or right of heading to determine the points in each half
                // double p2_angle = atan2(dy, dx);
                // //printf("Relative P2 Angle:  %f\n", p2_angle);
                // if (p2_angle >= 0) {
                //     l++;
                // } 
                // else {
                //     r++;
                // }
                if (dx * sin(p1.heading) - dy * cos(p1.heading) < 0) {
                    r++;  // Particle j is to the right of i
                }  else  {
                    l++;
                }
                // also check if particle in smaller radius circle
                // if (lhs < u_radius_sqrd) {
                //     n_close++;
                // }
            }
            //}
        }

        // total num within circle
        int n = l + r;
        //printf("P1 N: %d\n", n);
        //printf("P1 L: %d\n", l);
        //printf("P1 R: %d\n", r);

        // set colour
        p1.colour = get_colour(n, n_close);

        // set change in heading direction
        double d_phi = u_a + u_b * n * sign(r - l);
        //printf("Original P1 heading:  %f\n", p1.heading);
        p1.heading = std::fmod(p1.heading + d_phi, TAU);
        
        //printf("Delta P1 Heading:  %f\n", d_phi);
        //printf("New P1 heading:  %f\n", p1.heading);

        // apply force to get new x and y
        // TODO: OPTIMISATION HERE IS TO FILL A LOOKUP TABLE WITH VALUES!
        double x_change = cos(p1.heading) * u_velocity;
        double y_change = sin(p1.heading) * u_velocity;
        //printf("Original X Pos:  %f\n", p1.x);
        //printf("Original X Pos:  %f\n", p1.y);
        p1.x += x_change;
        p1.y += y_change;
        
        //printf("X Pos Delta:  %f\n", x_change);
        //printf("Y Pos Delta:  %f\n", y_change);
        //printf("New X Pos:  %f\n", p1.x);
        //printf("New Y Pos:  %f\n", p1.y);

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

    // delete memory from old state u_state, and then set pointer to new memory state
    delete[] move_order;
    // u_state = new_state;
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

