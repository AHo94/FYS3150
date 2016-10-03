#include "celestials.h"

Celestials::Celestials()
{
    mass = 0;
    distance = 0;
}

void Celestials::set_properties(double M, double R)
{
    // Initializes the properties of the object
    mass = M;
    distance = R;
}

double Celestials::get_mass(void)
{
    // Returns the mass of the object
    return mass;
}

double Celestials::get_distance(void)
{
    // Returns the distance to from the object to the Sun
    return distance;
}
