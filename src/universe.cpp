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


void add_to_send_buffer(particle_type* send_buffer, particle_type* p) {
    /**
    *  @brief Copy particle to given send buffer, to be sent to another process.
    */
    memcpy(send_buffer, p, sizeof(particle_type));
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


Miniverse::Miniverse(int num_proc, int rank, box_coord_type box_coords, int width, int height, MPI_Comm grid_comm, 
                     int neighbours[NUM_NEIGHBOURS], int num_particles, Universe* universe) {
    /**
     * @brief Miniverse constructor, a local square grid cell of the universe managing its own set of particles independently. 
     */

    // instantiate MPI parameters
    m_rank = rank;
    m_num_proc = num_proc;
    m_grid_comm = grid_comm;
    m_neighbours = neighbours;

    // create MPI particle type: id/x_coord/y_coord/heading/colour
	int num_part_fields = 5;
	MPI_Datatype particle_types[num_part_fields] = {MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
    int block_counts[num_part_fields] = {1, 4, 4, 4, 1};
    MPI_Aint offsets[num_part_fields], extent;
    offsets[0] = 0;
	MPI_Type_extent(MPI_INT, &extent);
    offsets[1] = 1 * extent;
    MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[2] = 4 * extent;
	MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[3] = 4 * extent;
	MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[4] = 4 * extent;
    MPI_Type_struct(num_part_fields, block_counts, offsets, particle_types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);

    // initialize miniverse parameters
    m_box = box_coords;
    m_width = width;
    m_height = height;
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
    m_universe -> u_rand_gen.seed((unsigned int)time(0));
    //m_universe -> u_rand_gen.seed(100); // FOR TIMING!
    std::uniform_real_distribution<float> uniform_rand(0.0f, 1.0f);

    // initialise all particles with random positions and headings
    for (int i = 0; i < m_num_particles; i++) {

        // instantiate new particle
        particle_type* p = new Particle();
        p->id = i + m_rank * m_num_particles;
        p->colour = green;
        //p->colour = (Colour) m_rank;

        // randomly place particles in box with random headings
        p->x = m_box.x0 + uniform_rand(m_universe->u_rand_gen) * m_width;
        p->y = m_box.y0 + uniform_rand(m_universe->u_rand_gen) * m_height;
        p->heading = uniform_rand(m_universe->u_rand_gen) * TAU;

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
    for(int j = 0; j < NUM_NEIGHBOURS; j++) {
        m_edge_send_counts[j] = 0;
        m_send_counts[j] = 0;
    }

    // Step 1: determine interactions between particles within own grid cell
    // O((n/p)^2) pairwise calculation
    for (particle_type* &p1 : m_particle_list) {

        // reset counts stored in particle to 0
        p1->l = 0;
        p1->r = 0;
        p1->n_close = 0;
        
        // interactions
        for (particle_type* &p2 : m_particle_list) {

            // exclude self from check
            if (p1->id == p2->id) {
                continue;
            }

            // update neighbour counts for p1
            CheckIfNeighbours(p1, p2);
        }

        // add particle to send buffer and edge particles list if in contact with grid cell edge
        std::array<int, 2> grid_cell_dir = CheckParticleEdgeContact(p1);
        int x_dir = grid_cell_dir[0], y_dir = grid_cell_dir[1];

        // check x direction
        if (x_dir != -1) {
            m_edge_lists[x_dir].push_front(p1);
            AddToSendBuffer(&(m_send_buffer[x_dir][m_edge_send_counts[x_dir]++]), p1);
        }

        // check y direction
        if (y_dir != -1) {
            m_edge_lists[y_dir].push_front(p1);
            AddToSendBuffer(&(m_send_buffer[y_dir][m_edge_send_counts[y_dir]++]), p1);
        }
    }


    // Step 2: send and receive particles which were in contact with walls, and process interactions between them
    SendRecvParticles(m_edge_send_counts);
    
    // iterate through received data and edge list and determine interactions
    // between particles on each grid border
    for(int i = 0; i < NUM_NEIGHBOURS; i++) {

        int num_particles_recvd;
        MPI_Get_count(&(m_status[i]), mpi_particle_type, &num_particles_recvd);

        for (particle_type* &p1 : m_particle_list) {

            for(int j = 0; j < num_particles_recvd; j++) {

                // get received particle
                particle_type recvd_p = m_recv_buffer[i][j];

                // calculate neighbour counts for p1
                CheckIfNeighbours(p1, &recvd_p);
            }
        }
    }

    // Step 3: now we have final counts, update colours, velocities and headings, for all particles and then update their
    // positions, passing to other processes and removing from local list if they go outside of the local box.
    for (std::list<particle_type*>::const_iterator iter = m_particle_list.begin(), end = m_particle_list.end(); iter != end; ++iter) {

        particle_type* p = *iter;

        // total num within circle
        int n = p->l + p->r;

        // set colour
        p->colour = get_colour(n, p->n_close);
        // p->colour = (Colour) m_rank;

        // set change in heading direction
        p->heading += alpha + beta * n * sign(p->r - p->l);

        // apply force to get new x and y
        p->x += cos(p->heading) * velocity;
        p->y += sin(p->heading) * velocity;

        // wrap x
        if (p->x < 0) {
            p->x += u_width;
        } else if (p->x >= u_width) {
            p->x -= u_width;
        }

        // wrap y
        if (p->y < 0) {
            p->y += u_height;
        } else if (p->y >= u_height) {
            p->y -= u_height;
        }

        // check which grid cell new coords are in, add to send buffer and remove from list
        // if they are in a different cell
        int neighbour = CheckParticleBox(p);
        if (neighbour != -1) {

            // copy particle to be sent to other grid cell
            AddToSendBuffer(&(m_send_buffer[neighbour][m_send_counts[neighbour]++]), p);

            // remove particle from current list
            delete p;
            iter = m_particle_list.erase(iter);
            --iter;
        }
    }

    // send and receive particles to neighbouring processes
    SendRecvParticles(m_send_counts);
    //printf("WE FINISHED SECOND SEND RECEIVE DAWG!!!...\n");

    // add received particles to local list of particles
    for(int i = 0; i < NUM_NEIGHBOURS; i++) {

        int num_particles_recvd;
        MPI_Get_count(&(m_status[i]), mpi_particle_type, &num_particles_recvd);
        for(int j = 0; j < num_particles_recvd; j++) {

            // instantiate new particle
            particle_type* new_p = new Particle();
            particle_type recvd_p = m_recv_buffer[i][j];
            new_p->id = recvd_p.id;
            new_p->x = recvd_p.x;
            new_p->y = recvd_p.y;
            new_p->heading = recvd_p.heading;
            new_p->colour = recvd_p.colour;

            // add particle to linked list
            m_particle_list.push_front(new_p);
        }
    }

    // clear edge particle lists at end of each step
    for (int i = 0; i < NUM_NEIGHBOURS; i++) {

        // DEBUG: set particle colours on edge to magenta
        // for (particle_type* &p : m_edge_lists[i]) {
        //     p->colour = magenta;
        // }

        m_edge_lists[i].clear(); 
    }
}


void Miniverse::CheckIfNeighbours(particle_type* p1, particle_type* p2) {
    /**
    *  @brief Determine whether p2 is neighbouring p1, and update p1 neighbour
              counts accordingly.
    */
    
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
                p1->r++;
            }  else  {
                p1->l++;
            }

            // also check if particle in smaller radius circle, just for colouring
            if (lhs < close_radius_sqrd) {
                p1->n_close++;
            }
        }
    }
}


void Miniverse::SendRecvParticles(int* send_counts) {
    /**
    *  @brief Send all particles added to send buffer, and receive corresponding particles from neighbours
              in receive buffer.
    */

    for(int i = 0; i < NUM_NEIGHBOURS; i++){

        // send this process' particles to corresponding neighbor
        MPI_Isend(m_send_buffer[i], send_counts[i], mpi_particle_type, m_neighbours[i], 0, m_grid_comm, &(m_request[i]));

        //printf("Direction %d send counts: %d\n", i, m_send_counts[i]);

        //Receive own data from neighbor
        MPI_Recv(m_recv_buffer[i], COMM_BUFFER_SIZE, mpi_particle_type, MPI_ANY_SOURCE, 0, m_grid_comm, &(m_status[i]));
    }

    // wait for non-blocking send on all processes.
    MPI_Waitall(NUM_NEIGHBOURS, m_request, MPI_STATUS_IGNORE);
}


std::array<int, 2> Miniverse::CheckParticleEdgeContact(particle_type* p) {
    /**
    *  @brief Check if the radius about a particle is in contact 
              with the edges of a grid cell, so as to consider it for 
              checking neighbour counts in other processes. 
    **/

    // x, y edge contacts
    std::array<int, 2> edge_contacts = {-1, -1};

    // get deltas between particle and four edges
    float dx_left = m_box.x0 - p->x;
    float dx_right = m_box.x1 - p->x;
    float dy_down = m_box.y0 - p->y;
    float dy_up = m_box.y1 - p->y;

    // check x edges
    if (abs(dx_left) <= radius) {
        edge_contacts[0] = LEFT;

    } else if (abs(dx_right) <= radius) {
        edge_contacts[0] = RIGHT;
    }

    // check y edges
    if (abs(dy_down) <= radius) {
        edge_contacts[1] = DOWN;

    } else if (abs(dy_up) <= radius) {
        edge_contacts[1] = UP;
    }

    return edge_contacts;
}


int Miniverse::CheckParticleBox(particle_type* p) {
    /**
    *  @brief Check if particle has moved into one of the surrounding grid cells. 
    **/

    if(p->x < m_box.x0)  return LEFT;
    if(p->x >= m_box.x1)  return RIGHT;
    if(p->y < m_box.y0) return DOWN;
    if(p->y >= m_box.y1) return UP;

    // TODO: CHECK DIAGONALS AS WELL!!!

    // return flag that the particle is still in the same box
    return -1;
}

void Miniverse::AddToSendBuffer(particle_type *send_buffer, particle_type *p) {
    /**
    *  @brief Copies the memory of a particle to a specified send buffer, to
              be sent to another grid cell process.
    **/

    memcpy(send_buffer, p, sizeof(particle_type));
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

    printf("IMA BOUT TO CLEAN THIS MOTHER\n");
    // clear particle list
    for (auto& particle : m_particle_list) {
        delete particle;
    }
    m_particle_list.clear(); 

    // close io parser
	io_parse_obj -> CloseOutFile();
	delete io_parse_obj;

    // clear MPI particle type
    MPI_Type_free(&mpi_particle_type);
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

