#include "solarsystem.h"
#include <iostream>

using namespace std;

SolarSystem::SolarSystem()
{
    new_tot_energy = 0;
    old_tot_energy = 0;
}

Celestials &SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass)
{
    m_bodies.push_back( Celestials(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

int SolarSystem::NumberofBodies()
{
    return m_bodies.size();
}

void SolarSystem::CalculateAccelerationAndEnergy(){
    /*
     * This function calculates the new accelerations for the objects.
     * Will also calculate the total energy of the system.
     * An if test will be used to see if the total energy of the system is conserved.
     * If total energy is not conserved, the program stops.
     */
    m_kin_energy = 0;
    m_pot_energy = 0;
    old_tot_energy = new_tot_energy;

    // Reset forces/acceleration on all bodies
    for (Celestials &body : m_bodies){
        body.acceleration.zeros();
    }

    double four_pi2 = 4*pi_m*pi_m;  // Defines 4pi^2, used multiple times in for-loop
    for (int i=0; i<NumberofBodies(); i++){
        Celestials &body1 = m_bodies[i];
        for (int j=i+1; j<NumberofBodies(); j++){
            Celestials &body2 = m_bodies[j];
            vec3 dRvector = body2.position - body1.position;
            double dR = dRvector.length();

            vec3 factor = -four_pi2*dRvector/(dR*dR*dR);
            body1.acceleration = factor*body2.mass;
            body2.acceleration = factor*body1.mass;

            m_pot_energy -= four_pi2*body1.mass*body2.mass;
        }
    m_kin_energy += 0.5*body1.mass*body1.velocity.dot(body1.velocity);
    }
    new_tot_energy = m_kin_energy + m_pot_energy;   // New total energy
    if (old_tot_energy != 0){
        if (fabs(new_tot_energy - old_tot_energy) > 1e-6){
            cout << "Total energy not conserved, stopping program" << endl;
            terminate();
        }
    }
}

void SolarSystem::write_file(string filename){
    /*
     * Writes the positions of the celestial bodies to a file.
     */
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }
    for(Celestials &body : m_bodies){
        m_file << body.position.x() << " " << body.position.y() << " " << body.position.z() << " ";
    }
    m_file << "\n";
}

std::vector<Celestials> &SolarSystem::bodies()
{
    return m_bodies;
}
