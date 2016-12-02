#ifndef METROPOLIS_QUANTUM_H
#define METROPOLIS_QUANTUM_H
#include "vec3.h"
#include "wavefunctions.h"

class Metropolis_Quantum
{
public:
    Metropolis_Quantum();
    double CalculateStepLength(vec3 r1, vec3 r2, double alpha, double omega, double *s);
    double LaplaceAnalytic(vec3 r1, vec3 r2, double alpha, double omega);
    double LaplaceOperator(vec3 r1, vec3 r2, double alpha, double omega);
    void Metropolis_T1(int MC_cycles, double alpha, double omega, int Analytic = 0);
    void Metropolis_T2(int MC_cycles, double alpha, double beta, double omega);
    void Write_file(double alpha, double omega);

private:
    double EnergyExpectation;
    double EnergyExpectationSquared;
    double MeanDistanceExpectation;
};

#endif // METROPOLIS_QUANTUM_H
