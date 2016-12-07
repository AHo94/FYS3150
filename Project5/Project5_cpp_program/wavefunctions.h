#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H
#include "vec3.h"
#include <cmath>
class Wavefunctions
{
public:
    Wavefunctions(int WaveFuncChoice);
    int FuncFactor;  // If FuncFactor = 0, returns first trial wave function. If FuncFactor = 1, returns second wavefunc.
    double operator()(vec3 r1, vec3 r2, double omega, double alpha, double beta=0)\
                {return ((1-FuncFactor)*Wavefunction_T1(r1, r2, alpha, omega) \
                + FuncFactor*Wavefunction_T2(r1, r2, alpha, beta, omega));}
    double Wavefunction_T1(vec3 r1, vec3 r2, double alpha, double omega);
    double Wavefunction_T2(vec3 r1, vec3 r2, double alpha, double beta, double omega);
};

#endif // WAVEFUNCTIONS_H
