#include "eulermethod.h"
#include "solarsystem.h"
eulermethod::eulermethod(double dt):
    m_dt(dt)
{
}

void eulermethod::Euler3D_step(SolarSystem &system)
{
    system.CalculateAccelerationAndEnergy();

    for (Celestials &body : system.bodies()){
        body.position += body.velocity*m_dt;
        body.velocity += body.acceleration*m_dt;
    }
}
