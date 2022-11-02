//
// Created by joeho on 13/08/2022.
//

#include "universe.h"
#include "particle.h"
#include "ioparser.h"
#include <mpi.h>
#include <string> 
#include <sstream>

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


Miniverse::Miniverse(int num_proc, int rank, box_coord_type box_coords, int width, int height, MPI_Comm grid_comm, int num_particles, Universe* universe) {
    /**
     * @brief Miniverse constructor, a local square grid cell of the universe managing its own set of particles independently. 
     */

    // initialize miniverse parameters
    m_rank = rank;
    m_num_proc = num_proc;
    m_box = box_coords;
    m_width = width;
    m_height = height;
    m_grid_comm = grid_comm;
    m_num_particles = num_particles;
    m_universe = universe;

    // unpack universe settings
    u_num_particles = m_universe -> u_num_particles;
    u_width = m_universe -> u_width;
    u_height = m_universe -> u_height;
    density = m_universe -> u_density;
    alpha = m_universe -> u_a;
    beta = m_universe -> u_b;
    gamma = m_universe -> u_g;
    radius = m_universe -> u_radius;
    velocity = m_universe -> u_velocity; 
    radius_sqrd = m_universe -> u_radius_sqrd;
    close_radius = m_universe -> u_close_radius;
    close_radius_sqrd = m_universe -> u_close_radius_sqrd;

    // initialise the start state using the given parameters
    InitState();

    // instante IO (i.e. class which handles keeping log file stream open, formatting for reading/writing)
	// to write to a seperate file
    std::ostringstream out_file_name;
    out_file_name << "../results/out_display_" << std::to_string(m_rank) << ".bin";
	io_parse_obj = new IOParser(out_file_name.str(), "out_log.bin");
	io_parse_obj -> OpenOutFile();

    // write initial state to out file
    WriteLocalParticlesToOutFile();

    printf("Hello from process %d out of %d, I have Height: %d, Width: %d, Local Num Parts: %d, x0: %lf, x1: %lf, y0: %lf, y1: %lf\n", 
	m_rank, m_num_proc, m_width, m_height, m_num_particles, m_box.x0, m_box.x1, m_box.y0, m_box.y1);
}

void::Miniverse::InitState() {
    /**
     * @brief Place all particles randomly in the local miniverse box.
     */

    // get seeded uniform random distribution
    // m_universe.u_rand_gen.seed((unsigned int)time(0));
    m_universe -> u_rand_gen.seed(100); // FOR TIMING!
    std::uniform_real_distribution<float> uniform_rand(0.0f, 1.0f);

    // initialise all particles with random positions and headings
    for (int i = 0; i < m_num_particles; i++) {
        //particle_type* p = new Particle();
        particle_type* p = new Particle();
        p->id = i + m_rank * m_num_particles;
        p->colour = green;

        // TODO: ENSURE THIS IS WITHIN THE BOUNDS OF OUR BOX!!!
        p->x = m_box.x0 + uniform_rand(m_universe -> u_rand_gen) * m_width;
        p->y = m_box.y0 + uniform_rand(m_universe -> u_rand_gen) * m_height;
        p->heading = uniform_rand(m_universe -> u_rand_gen) * TAU;

        // add particle to linked list
        m_particle_list.push_front(p);
    }
}

void Miniverse::Step() {
    /**
     * @brief Perform a time step update of the Universe simulation.
     * 
     */

    // reset send counts for this iteration
    for(int j = 0; j < 4; j++) {
        m_send_counts[j] = 0;
    }

    // Step 1: determine interactions between particles within own grid cell
    // O((n/p)^2) pairwise calculation
    //particle_type* p1 =  m_particle_list.front();
    for (particle_type* &p1 : m_particle_list) {

        // new particle?
        // Particle p1_new;

        // counts of particles within left and right semi-circle
        int l = 0;
        int r = 0;
        int n_close = 0;
        
        // interactions
        //particle_type* p2 = m_particle_list.front();
        for (particle_type* &p2 : m_particle_list) {

            // exclude self from check
            if (p1->id == p2->id) {
                continue;
            }

            // get deltas
            float dx = p2->x - p1->x;
            float dy = p2->y - p1->y;

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
            if (abs(dx) <= radius && abs(dy) <= radius) {

                // check if particle in u_radius circle
                float lhs = dx * dx + dy * dy;
                if (lhs <= radius_sqrd) {

                    // check if point is to the left or right of heading to determine the points in each half
                    if (dx * sin(p1->heading) - dy * cos(p1->heading) < 0) {
                        r++;
                    }  else  {
                        l++;
                    }

                    // also check if particle in smaller radius circle, just for colouring
                    if (lhs < close_radius_sqrd) {
                        n_close++;
                    }
                }
            }
        }

        // total num within circle
        int n = l + r;

        // store n in particle!!

        // TODO: ALSO HERE ADD PARTICLE TO SENDBUFFER LISTS IF IT IS IN RADIUS OF GRID CELL EDGE!

        // TODO: THE FOLLOWING UPDATES NEED TO BE MOVED TO A LATER FINAL STEP ONCE ALL CHECKS HAVE BEEN DONE!!!

        // set colour
        p1->colour = get_colour(n, n_close);

        // set change in heading direction
        p1->heading += alpha + beta * n * sign(r - l);

        // apply force to get new x and y
        p1->x += cos(p1->heading) * velocity;
        p1->y += sin(p1->heading) * velocity;

        // wrap x
        if (p1->x < 0) {
            p1->x += u_width;
        } else if (p1->x >= u_width) {
            p1->x -= u_width;
        }

        // wrap y
        if (p1->y < 0) {
            p1->y += u_height;
        } else if (p1->y >= u_height) {
            p1->y -= u_height;
        }
    }

    // Step 2: send and receive particles which were in contact with walls, and process interactions between them, adding
    // count to n (may need to faff with some pointer stuff here for updating particles between lists)

        /*
    //Send neighbor's data to the corresponding neighbor.
    for(j = 0; j < 4; j++){
        MPI_Isend(sendBuffer[j], sendCounts[j], mpiPartType, neighbours[j], 0, gridComm, &(request[j]));
    }
    //Receive own data from neighbor.
    for(j = 0; j < 4; j++){
        MPI_Recv(recvBuffer[j], COMM_BUFFER_SIZE, mpiPartType, MPI_ANY_SOURCE, 0, gridComm, &(status[j]));
    }
    MPI_Waitall(4, request, MPI_STATUS_IGNORE); //Wait for non-blocking send completion.
    */

    // Step 3: now we have final counts, update colours velocities and headings, for all particles and then update their
    // positions, passing to other processes if they go outside of the local box.



}

void Miniverse::WriteLocalParticlesToOutFile() {
    /**
    *  @brief Writes the particle list to the io parser out file.
    */

    io_parse_obj->WriteStateToOutFile(m_particle_list, m_particle_list.size());
}

void Miniverse::Clean() {
    /**
     * @brief Clear memory for the particle list and io parser.
     * 
     */

    
    // clear particle list
    for (auto& particle : m_particle_list) {
        delete particle;
    }
    m_particle_list.clear(); 

    // close io parser
	io_parse_obj -> CloseOutFile();
	delete io_parse_obj;
}

std::list<particle_type*> Miniverse::GetParticleList() {
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

