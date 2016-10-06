#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include "celestials.h"
#include "vec3.h"
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
class SolarSystem
{
public:
    SolarSystem();
    double pi_m = acos(-1);
    Celestials &createCelestialBody(vec3 position, vec3 velocity, double mass);

    int NumberofBodies();
    void CalculateAccelerationAndEnergy();
    void write_file(std::string filename);

    std::vector<Celestials> &bodies();
private:
    std::vector<Celestials> m_bodies;
    std::ofstream m_file;
    double m_kin_energy;
    double m_pot_energy;
    double new_tot_energy;
    double old_tot_energy;
    vec3 angular_momentum;
};

#endif // SOLARSYSTEM_H
