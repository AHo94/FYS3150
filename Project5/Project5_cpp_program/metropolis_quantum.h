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
    void Metropolis_T1(int MC_cycles, Wavefunctions &WaveFunc, double alpha, double omega, int CoulombInt, int Analytic = 0);
    void Metropolis_T2(int MC_cycles, Wavefunctions &WaveFunc, double alpha, double beta, double omega);
    void Write_file(int MC_cycles, std::string filename, double omega, double alpha ,double beta=0);

private:
    std::fstream ofile_global;
    double EnergyExpectation;
    double EnergyExpectationSquared;
    double MeanDistanceExpectation;
};

#endif // METROPOLIS_QUANTUM_H
