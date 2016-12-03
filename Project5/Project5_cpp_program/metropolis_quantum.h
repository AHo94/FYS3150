#ifndef METROPOLIS_QUANTUM_H
#define METROPOLIS_QUANTUM_H
#include "vec3.h"
#include "wavefunctions.h"
#include <fstream>
#include <string>
#include <chrono>
#include <random>
#include <iomanip>

class Metropolis_Quantum
{
public:
    Metropolis_Quantum();
    double CalculateStepLength(vec3 r1, vec3 r2, double alpha, double omega, double *s);
    double LaplaceAnalytic(vec3 r1, vec3 r2, double alpha, double omega);
    double LaplaceOperator(Wavefunctions &WaveFunc, vec3 r1, vec3 r2, double alpha, double omega);

    void Metropolis_T1(int MC_cycles, Wavefunctions &WaveFunc, double *ExpectationValues, double alpha,\
                       double omega, int CoulombInt, int Analytic = 0);

    void Metropolis_T2(int MC_cycles, Wavefunctions &WaveFunc, double *ExpectationValues\
                       , double alpha, double beta, double omega);

private:
};

#endif // METROPOLIS_QUANTUM_H
