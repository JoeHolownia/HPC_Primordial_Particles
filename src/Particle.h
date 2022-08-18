//
// Created by joeho on 13/08/2022.
//

#ifndef PROJECT_PARTICLE_H
#define PROJECT_PARTICLE_H

enum colours {
    green,
    brown,
    magenta,
    blue,
    yellow
};

struct Particle {
    float x, y;  // x and y coordinates
    float heading;  // direction particle is facing (+ is to the right, - to the left)
    uint8_t colour;  // current colour
};

#endif //PROJECT_PARTICLE_H
