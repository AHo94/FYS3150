#ifndef EULERMETHOD_H
#define EULERMETHOD_H

class eulermethod
{
public:
    eulermethod(double dt);
    double m_dt;
    void Euler3D_step(class SolarSystem &system);
};

#endif // EULERMETHOD_H
