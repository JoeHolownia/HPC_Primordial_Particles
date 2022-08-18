//
// Created by joeho on 13/08/2022.
//

#include "Universe.h"
#include "Particle.h"

Universe::Universe(int num_particles, float density, int width, int height, int radius, float a, float b, float velocity) {
    /*
     * Instantiate Universe.
     */

    // initialize universe parameters
    u_num_particles = num_particles;
    u_density = density;
    u_width = width;
    u_height = height;
    u_radius = radius;
    u_a = a;
    u_b = b;
    u_velocity = velocity;


    //
    // initialise the start state using the given parameters.
    InitState();
    //m_rand_gen.seed((unsigned int)time(0));

    //SetPopulation(num_types, num_particles);
    //SetSize(float(width), float(height));

}

void::Universe::InitState() {
    /**
     * @brief Place all particles randomly in the Universe box using the given density and box size.
     * 
     */
 
  std::uniform_real_distribution<float> rand_uni(0.0f, 1.0f);
  std::normal_distribution<float>       rand_norm(0.0f, 1.0f);
  for (size_t i = 0; i < m_particles.size(); ++i) {
    Particle& p = m_particles[i];
    p.type = uint8_t(rand_type(m_rand_gen));
    p.x = (rand_uni(m_rand_gen)*0.5f + 0.25f) * m_width;
    p.y = (rand_uni(m_rand_gen)*0.5f + 0.25f) * m_height;
    p.vx = rand_norm(m_rand_gen) * 0.2f;
    p.vy = rand_norm(m_rand_gen) * 0.2f;
  }
}

void Universe::Step() {
    /**
     * @brief Perform a time step update of the Universe simulation.
     * 
     */

    // naive O(n^2) pairwise calculation
    for (int i = 0; i < u_num_particles; i++) {

        // current particle
        Particle &p = u_state[i];

        // interactions
        for (int j = 0; j < u_num_particles; j++) {

            // other particle
            const Particle &q = u_state[j];


        }
    }
}

void Universe::Clean() {
    /**
     * @brief Clear memory for all particle, i.e. clear the state.
     * 
     */

    // iterate through state array 

}

