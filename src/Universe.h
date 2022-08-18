//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_UNIVERSE_H
#define PROJECT_UNIVERSE_H


class Universe {
public:

    Universe(int num_particles, float density, int width, int height, int radius, float a, float b, float velocity);
    void Step();
    void InitState();

private:
    int u_num_particles;  // number of particles
    float u_density;  // particle placement density
    int u_width;
    int u_height;
    int u_radius;  // radius
    float u_a;  // alpha
    float u_b;  // beta
    float u_velocity;  // fixed-velocity
    //std::array<Particle> u_state;
    Particle u_state [];
};


#endif //PROJECT_UNIVERSE_H
