#include "metropolis_quantum.h"
#include <chrono>
#include <random>
#include <iostream>
using namespace std;
Metropolis_Quantum::Metropolis_Quantum()
{
    EnergyExpectation = 0;
    EnergyExpectationSquared = 0;
    MeanDistanceExpectation = 0;
}

double Metropolis_Quantum::CalculateStepLength(vec3 r1, vec3 r2, double alpha, double omega, double *s){
    // Calculate step length
    double a = 0;
    for (int i = 0; i<6; i++){
        a += s[i]*s[i];
    }
    double b = 2*(r1[0]*s[0] + r1[1]*s[1] + r1[2]*s[2] + r2[0]*s[3] + r2[1]*s[4] + r2[2]*s[5]);
    double c = log(0.5)/(alpha*omega);
    return (-b + sqrt(b*b - 4*a*c))/(2*a);
}

double Metropolis_Quantum::LaplaceAnalytic(vec3 r1, vec3 r2, double alpha, double omega)
{
    // Function that calculates the analytical expression of the Laplace operator
    double AlphaOmega = alpha*omega;
    return  -0.5*(pow(AlphaOmega,2)*r1.lengthSquared() - 3*AlphaOmega) \
            - 0.5*(pow(AlphaOmega,2)*r2.lengthSquared() - 3*AlphaOmega);
}

double Metropolis_Quantum::LaplaceOperator(vec3 r1, vec3 r2, double alpha, double omega)
{
    // Function that calculates the Laplace operator numerically
    double SecondDerivative = 0;
    double dr = 1e-5;
    double wavefunc = WaveInstance.Wavefunction_T1(r1, r2, alpha, omega);
    for (int i=0; i<3; i++){
        vec3 rchange(0,0,0);
        rchange[i] = dr;
        SecondDerivative -= (WaveInstance.Wavefunction_T1(r1+rchange, r2, alpha, omega) -\
                2*wavefunc + WaveInstance.Wavefunction_T1(r1-rchange, r2, alpha, omega));
        SecondDerivative -= (WaveInstance.Wavefunction_T1(r1, r2+rchange, alpha, omega) -\
                2*wavefunc + WaveInstance.Wavefunction_T1(r1, r2-rchange, alpha, omega));
    }
    return 0.5*SecondDerivative/(wavefunc*(dr*dr));
}

void Metropolis_Quantum::Metropolis_T1(int MC_cycles, double alpha, double omega, int Analytic){
    /* Function that solves the Metropolis method for the first trial function.
    * Argument "Analytic" is set to zero by default. Acts like an optional argument.
    * If Analytic = 1, function uses analytic Laplace operator.
    * If Analytic = 0 (or an arbitrary number), uses numerical Laplace operator.
    */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distrUniform(0, 1.0);
    std::uniform_real_distribution<double> distr(-1.0, 1.0);
    // Random initial starting position for both electrons
    vec3 r1(distr(generator),distr(generator),distr(generator));
    vec3 r2(distr(generator),distr(generator),distr(generator));

    double E_local = 0;
    double EnergySum = 0;
    double EnergySquaredSum = 0;
    double MeanDistance = 0;
    double NewWavefuncSquared = 0;
    double omega2 = omega*omega;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);

    int counter = 0;
    // Runs the Monte Carlo cycles
    double OldWavefuncSquared = pow(WaveInstance.Wavefunction_T1(r1, r2, alpha, omega), 2);
    for (int cycle=0; cycle<MC_cycles; cycle++){
        // Running Monte Carlo cycles
        vec3 r1_new(0,0,0);
        vec3 r2_new(0,0,0);
        for (int j=0; j<3; j++){
            rDistr[j] = distr(generator);
            rDistr[j+3] = distr(generator);
            r1_new[j] = r1[j] +step_length*rDistr[j];
            r2_new[j] = r2[j] +step_length*rDistr[j+3];
        }
        NewWavefuncSquared = pow(WaveInstance.Wavefunction_T1(r1_new, r2_new, alpha, omega), 2);
        if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);
        if (Analytic == 1){
            E_local = LaplaceAnalytic(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        }
        else{
            E_local = LaplaceOperator(r1, r2, alpha, omega) + \
            0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        }
        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
        MeanDistance += (r1-r2).length();
    }

    // Adding the energies and mean distance in their arrays
    EnergyExpectation = EnergySum;
    EnergyExpectationSquared = EnergySquaredSum;
    MeanDistanceExpectation = MeanDistance;
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Energy = "<< EnergyExpectation/(MC_cycles) << endl;
    cout << "Variance = " << EnergyExpectationSquared/MC_cycles - \
            EnergyExpectation*EnergyExpectation/MC_cycles/MC_cycles << endl;;
    cout << "Accepted configs = " << (double)counter/MC_cycles << endl;

}


