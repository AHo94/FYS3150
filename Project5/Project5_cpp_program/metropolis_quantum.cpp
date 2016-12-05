#include "metropolis_quantum.h"
#include <iostream>
using namespace std;

Metropolis_Quantum::Metropolis_Quantum()
{
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

double Metropolis_Quantum::LaplaceOperator(Wavefunctions &WaveFunc, vec3 r1, vec3 r2, double alpha, double omega)
{
    // Function that calculates the Laplace operator numerically
    double SecondDerivative = 0;
    double dr = 1e-5;
    double wavefunc = WaveFunc(r1, r2, omega, alpha);
    for (int i=0; i<3; i++){
        vec3 rchange(0,0,0);
        rchange[i] = dr;
        SecondDerivative -= (WaveFunc(r1+rchange, r2 , omega, alpha) -\
                2*wavefunc + WaveFunc(r1-rchange, r2, omega, alpha));
        SecondDerivative -= (WaveFunc(r1, r2+rchange, omega, alpha) -\
                2*wavefunc + WaveFunc(r1, r2-rchange, omega, alpha));
    }
    return 0.5*SecondDerivative/(wavefunc*(dr*dr));
}

void Metropolis_Quantum::Metropolis_T1(int MC_cycles, Wavefunctions &WaveFunc, double *ExpectationValues, double alpha,
                                       double omega, int ExactEne, int Analytic){
    /* Function that solves the Metropolis method for the first trial function.
    * Argument "Analytic" is set to zero by default. Acts like an optional argument.
    * If Analytic = 1, function uses analytic Laplace operator.
    * If Analytic = 0 (or an arbitrary number), uses numerical Laplace operator.
    */
    if (ExactEne != 0 && ExactEne != 1){
        cout << "ExactEne not set properly. Currently = " << ExactEne << endl;
        cout << "Try: ExactEne = 0 or 1" << endl;
        exit(1);
    }
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
    double alpha2 = alpha*alpha;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);

    int counter = 0;
    // Runs the Monte Carlo cycles

    double OldWavefuncSquared = pow(WaveFunc(r1, r2, omega, alpha), 2);
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
        NewWavefuncSquared = pow(WaveFunc(r1_new, r2_new, omega, alpha), 2);
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
            if (ExactEne == 1){
                E_local = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha2) + 3*alpha*omega;
            }
            else{
                E_local = LaplaceOperator(WaveFunc, r1, r2, alpha, omega) + \
                    0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
            }
        }
        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
        MeanDistance += (r1-r2).length();
    }

    // Adding the energies and mean distance in their arrays
    ExpectationValues[0] = EnergySum;
    ExpectationValues[1] = EnergySquaredSum;
    ExpectationValues[2] = MeanDistance;
    ExpectationValues[3] = counter;
}

void Metropolis_Quantum::Metropolis_T2(int MC_cycles, Wavefunctions &WaveFunc, double *ExpectationValues
                                       , double alpha, double beta, double omega)
{
    // Function that solves the Metropolis method for the second trial wave function.

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
    double r_12 = 0;
    double omega2 = omega*omega;
    double alpha2 = alpha*alpha;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);

    int counter = 0;
    double OldWavefuncSquared = pow(WaveFunc(r1, r2, omega, alpha, beta), 2);
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

        NewWavefuncSquared = pow(WaveFunc(r1_new, r2_new, omega, alpha,  beta), 2);
        if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);
        r_12 = (r1-r2).length();
        E_local = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha2) + 3*alpha*omega + 1/r_12\
                + (1/(2*(1+beta*r_12)*(1+beta*r_12)))*(alpha*omega*r_12 - 1/(2*(1+beta*r_12)*(1+beta*r_12))\
                - 2/r_12 + 2*beta/(1+beta*r_12));

        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
        MeanDistance += r_12;
    }
    // Saving calculated values to an array
    ExpectationValues[0] = EnergySum;
    ExpectationValues[1] = EnergySquaredSum;
    ExpectationValues[2] = MeanDistance;
    ExpectationValues[3] = counter;
}

void Metropolis_Quantum::Metropolis_Virial(int MC_cycles, Wavefunctions &WaveFunc, double *ExpectationValues
                                           , double alpha, double beta, double omega, int CoulombInt)
{
    if (CoulombInt != 1 && CoulombInt != 0){
        cout << "CoulombInt = " << CoulombInt << endl;
        cout << "Invalid value, try CoulombInt = 1 or CoulombInt = 0" << endl;
        exit(1);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distrUniform(0, 1.0);
    std::uniform_real_distribution<double> distr(-1.0, 1.0);
    // Random initial starting position for both electrons
    vec3 r1(distr(generator),distr(generator),distr(generator));
    vec3 r2(distr(generator),distr(generator),distr(generator));
    double E_pot = 0;
    double E_tot = 0;
    double KineticSum = 0;
    double KineticSquaredSum = 0;
    double PotentialSum = 0;
    double PotentialSquaredSum = 0;
    double MeanDistance = 0;
    double NewWavefuncSquared = 0;
    double r_12 = 0;
    double omega2 = omega*omega;
    double alpha2 = alpha*alpha;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);

    int counter = 0;
    double OldWavefuncSquared = pow(WaveFunc(r1, r2, omega, alpha, beta), 2);
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

        NewWavefuncSquared = pow(WaveFunc(r1_new, r2_new, omega, alpha,  beta), 2);
        if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        step_length = CalculateStepLength(r1, r2, alpha, omega, rDistr);
        r_12 = (r1-r2).length();
        E_tot = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha2) + 3*alpha*omega \
                + CoulombInt*(1.0/r_12)\
                + (1/(2*(1+beta*r_12)*(1+beta*r_12)))*(alpha*omega*r_12 - 1/(2*(1+beta*r_12)*(1+beta*r_12))\
                - 2/r_12 + 2*beta/(1+beta*r_12));
        E_pot = + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared()) + CoulombInt*(1.0/r_12);

        KineticSum += E_tot;
        KineticSquaredSum += E_tot*E_tot;
        PotentialSum += E_pot;
        PotentialSquaredSum += E_pot*E_pot;
        MeanDistance += r_12;
    }

    // Adding the energies and mean distance in their arrays
    ExpectationValues[0] = KineticSum - PotentialSum;
    ExpectationValues[1] = KineticSquaredSum - PotentialSquaredSum;
    ExpectationValues[2] = PotentialSum;
    ExpectationValues[3] = PotentialSquaredSum;
    ExpectationValues[4] = MeanDistance;
    ExpectationValues[5] = counter;
}

