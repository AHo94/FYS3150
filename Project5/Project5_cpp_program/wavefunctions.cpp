#include "wavefunctions.h"

Wavefunctions::Wavefunctions()
{

}

double Wavefunctions::Wavefunction_T1(vec3 r1, vec3 r2, double alpha, double omega)
{
    // Function for first trial wavefunction
    return exp(-0.5*alpha*omega*(r1.lengthSquared() + r2.lengthSquared()));
}

double Wavefunctions::Wavefunction_T2(vec3 r1, vec3 r2, double alpha, double beta, double omega)
{
    // Function for the second trial wavefunction
    double r_12 = (r1-r2).length();
    double Wavefunc1 = Wavefunction_T1(r1, r2, alpha, omega);
    return Wavefunc1*exp(r_12/(2*(1+beta*r_12)));
}
