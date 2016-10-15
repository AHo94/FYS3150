#ifndef ODESOLVERS_H
#define ODESOLVERS_H

class ODEsolvers
{
public:
    ODEsolvers(double dt);
    double m_dt;
    void Euler_step(class SolarSystem &Ssystem);
    void EulerCromer(class SolarSystem &Ssystem);
    void Verlet(class SolarSystem &Ssystem);
    void Verlet_GR(class SolarSystem &Ssystem);
};

#endif // ODESOLVERS_H
