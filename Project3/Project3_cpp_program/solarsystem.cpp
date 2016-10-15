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
     * Also checks if the angular momentum is conserved.
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
            body1.acceleration += factor*body2.mass;
            body2.acceleration += factor*body1.mass;

            m_pot_energy -= four_pi2*body1.mass*body2.mass;
        }
    m_kin_energy = 0.5*body1.mass*body1.velocity.dot(body1.velocity);
    //angular_momentum = body1.mass*(body1.position.cross(body1.velocity));
    }
    new_tot_energy = m_kin_energy + m_pot_energy;   // New total energy
    if (old_tot_energy != 0){
        if (fabs(new_tot_energy - old_tot_energy) > 1e-4){
            //cout << fabs(new_tot_energy - old_tot_energy) << endl;
            cout << "Total energy not conserved, stopping program" << endl;
            //terminate();
        }
    }
}

void SolarSystem::CalculateAccelerationAndEnergy_GR(){
    /*
     * Function used to test the perihelion precession of Mercury using general relativity.
     * This assumes only two celestial bodies, the Sun and Mercury.
     * Also assumes that the Sun is the first celestial in m_bodies
     */
   // m_kin_energy = 0;
   // m_pot_energy = 0;
   // old_tot_energy = new_tot_energy;

    // Reset forces/acceleration on all bodies
    for (Celestials &body : m_bodies){
        body.acceleration.zeros();
    }

    double four_pi2 = 4*pi_m*pi_m;  // Defines 4pi^2, used multiple times in for-loop

    Celestials &Sun = m_bodies[0];
    Celestials &Mercury = m_bodies[1];
    vec3 dRvector = Mercury.position - Sun.position;
    double dR = dRvector.length();
    vec3 orbital_angular_momentum = Sun.position.cross(Mercury.velocity);
    //vec3 orbital_angular_momentum = dRvector.cross(Mercury.velocity);
    double l = orbital_angular_momentum.lengthSquared();
    vec3 factor = -four_pi2*dRvector/(dR*dR*dR);
    Sun.acceleration += factor*Mercury.mass*(1 + 3*l/(dR*dR*c*c));
    Mercury.acceleration += factor*Sun.mass*(1 + 3*l/(dR*dR*c*c));
    //cout << 3*l/(dR*dR*c*c) << endl;


    /*
    for (int i=0; i<NumberofBodies(); i++){
        Celestials &body1 = m_bodies[i];
        for (int j=i+1; j<NumberofBodies(); j++){
            Celestials &body2 = m_bodies[j];
            vec3 dRvector = body2.position - body1.position;
            double dR = dRvector.lengthSquared();

            vec3 factor = -four_pi2*dRvector/dR;
            dRvector.cross(body2.velocity)
            body1.acceleration += factor*body2.mass;
            body2.acceleration += factor*body1.mass;
            */
}

void SolarSystem::write_file(string filename, string NumSteps, string dt){
    // Writes the positions of the celestial bodies to a file.
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }
    cout << "Wrote shait " << endl;
    m_file << NumSteps << dt;
    for(Celestials &body : m_bodies){
        m_file << body.position.x() << " " << body.position.y() << " " << body.position.z() << " ";
        body.position.print("POSITIN   ");
    }
    m_file << "\n";
}

std::vector<Celestials> &SolarSystem::bodies()
{
    return m_bodies;
}
