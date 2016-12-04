#include "wavefunctions.h"
#include <iostream>
Wavefunctions::Wavefunctions(int WaveFuncChoice):
FuncFactor(WaveFuncChoice)
{
    if (FuncFactor != 0 && FuncFactor != 1){
        std::cout << "Instance call value invalid. Current value = " << FuncFactor << std::endl;
        std::cout << "Try value: 0 or 1" << std::endl;;
        exit(1);
    }
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
    //double Wavefunc1 = Wavefunction_T1(r1, r2, alpha, omega);
    double Wavefunc1 = exp(-0.5*alpha*omega*(r1.lengthSquared() + r2.lengthSquared()));
    return Wavefunc1*exp(r_12/(2*(1+beta*r_12)));
}
