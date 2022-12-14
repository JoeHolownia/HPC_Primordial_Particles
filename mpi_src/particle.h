//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_PARTICLE_H
#define PROJECT_PARTICLE_H

enum Colour {
    green,
    brown,
    magenta,
    blue,
    yellow
};

struct Particle {
    int id;  // unique particle id
    int l;  // left neighbour counts
    int r;  // right neighbour counts
    int n_close;  // close neighbour counts
    float x, y;  // x and y coordinates
    float heading;  // direction particle is facing (+ is to the right, - to the left)
    Colour colour;  // current colour

    // bool operator==(const Particle* rhs) const {
    //     return id == rhs->id;
    // }
};

struct BoxCoords {
    float x0;
    float x1;
    float y0;
    float y1;
} ;

// struct part_cord {
//     float x;
//     float y;
//     float vx;
//     float vy;
// } ;

typedef struct BoxCoords box_coord_type;
typedef struct Particle particle_type;
// typedef struct part_cord pcord_t ;

#endif //PROJECT_PARTICLE_H
