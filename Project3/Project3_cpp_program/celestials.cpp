#include "celestials.h"

Celestials::Celestials(vec3 pos, vec3 vel, double in_mass)
{
    position = pos;
    velocity = vel;
    mass = in_mass;
}

Celestials::Celestials(double x, double y, double z, double vx, double vy, double vz, double in_mass)
{
    position = vec3(x,y,z);
    velocity = vec3(vx, vy, vz);
    mass = in_mass;
}

