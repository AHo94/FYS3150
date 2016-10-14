#ifndef ODESOLVERS_H
#define ODESOLVERS_H

class ODEsolvers
{
public:
    ODEsolvers(double dt);
    double m_dt;
    void Euler_step(class SolarSystem &system);
    void EulerCromer(class SolarSystem &system);
    void Verlet(class SolarSystem &system);
    void Verlet_GR(class SolarSystem &system);
};

#endif // ODESOLVERS_H
