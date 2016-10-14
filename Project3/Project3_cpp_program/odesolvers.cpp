#include "odesolvers.h"
#include "solarsystem.h"

ODEsolvers::ODEsolvers(double dt):
m_dt(dt)
{

}

void ODEsolvers::Euler_step(SolarSystem &Ssystem){
    // Solving Euler's method
    Ssystem.CalculateAccelerationAndEnergy();

    for (Celestials &body : Ssystem.bodies()){
        body.position += body.velocity*m_dt;
        body.velocity += body.acceleration*m_dt;
    }
}

void ODEsolvers::EulerCromer(SolarSystem &Ssystem){
    // Solving Euler-Cromer's method
    Ssystem.CalculateAccelerationAndEnergy();

    for (Celestials &body : Ssystem.bodies()){
        body.velocity += body.acceleration*m_dt;
        body.position += body.velocity*m_dt;
    }
}

void ODEsolvers::Verlet(SolarSystem &Ssystem){
    // Solving Verlet's method
    Ssystem.CalculateAccelerationAndEnergy();

    double dt2_2 = m_dt*m_dt/2.0;
    double dt_2 = m_dt/2.0;

    for (Celestials &body : Ssystem.bodies()){
        body.position += body.velocity*m_dt + dt2_2*body.acceleration;
        vec3 old_acc = body.acceleration;
        Ssystem.CalculateAccelerationAndEnergy();
        body.velocity += dt_2*(body.acceleration + old_acc);
    }
}

void ODEsolvers::Verlet_GR(SolarSystem &Ssystem){
    // Verlet's method, used to test general relativity.

    Ssystem.CalculateAccelerationAndEnergy_GR();
    double dt2_2 = m_dt*m_dt/2.0;
    double dt_2 = m_dt/2.0;

    for (Celestials &body : Ssystem.bodies()){
        body.position += body.velocity*m_dt + dt2_2*body.acceleration;
        vec3 old_acc = body.acceleration;
        Ssystem.CalculateAccelerationAndEnergy_GR();
        body.velocity += dt_2*(body.acceleration + old_acc);
    }

}
