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
    double x, y;  // x and y coordinates
    double heading;  // direction particle is facing (+ is to the right, - to the left)
    Colour colour;  // current colour
};

#endif //PROJECT_PARTICLE_H
