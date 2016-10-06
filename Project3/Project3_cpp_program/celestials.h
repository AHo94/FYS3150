#ifndef CELESTIALS_H
#define CELESTIALS_H
#include <vec3.h>

class Celestials
{
public:
    vec3 position;
    vec3 velocity;
    vec3 acceleration;
    double mass;

    Celestials(vec3 position, vec3 velocity, double mass);
    Celestials(double x, double y, double z, double vx, double vy, double vz, double mass);
};

#endif // CELESTIALS_H
