#include <iostream>
#include <cmath>
#include "vec3.h"
#include "celestials.h"
#include <fstream>

using namespace std;

int main()
{
    // Masses of the celestial bodies in the solar system
    double M_sun, M_earth, M_jupiter, M_mars;
    M_sun = 2*pow(10,30);
    M_earth = 6*pow(10,24);
    M_jupiter = 1.9*pow(10,27);
    M_mars = 6.6*pow(10,23);

    double AU = 1.5*pow(10,11); // Defining one astronomical unit in meters
    double R_earth = 1*AU;

    Celestials Earth;
    Earth.set_properties(M_earth, R_earth);
    cout << "Mass of Earth = " << Earth.get_mass() << endl;;


    return 0;
}
