#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H
#include "vec3.h"
#include <cmath>

class Wavefunctions
{
public:
    Wavefunctions();
    double Wavefunction_T1(vec3 r1, vec3 r2, double alpha, double omega);
    double Wavefunction_T2(vec3 r1, vec3 r2, double alpha, double beta, double omega);
};

#endif // WAVEFUNCTIONS_H
