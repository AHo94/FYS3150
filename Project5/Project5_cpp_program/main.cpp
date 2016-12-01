#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
#include "vec3.h"
using namespace std;

ofstream ofile_global;

double TrialWavefunc_T1(vec3 r1, vec3 r2, double alpha, double omega){
    // Function that calculates the first trial wave function
    return exp(-0.5*alpha*omega*(r1.lengthSquared() + r2.lengthSquared()));
}

double TrialWavefunc_T2(vec3 r1, vec3 r2, double alpha, double beta, double omega){
    // Function that calculates the second trial wave function
    double r_12 = (r1-r2).length();
    double Wavefunc1 = TrialWavefunc_T1(r1, r2, alpha, omega);
    return Wavefunc1*exp(r_12/(2*(1+beta*r_12)));
}

double StepLength(vec3 r1, vec3 r2, double alpha, double omega, double *r){
    // Calculates the optimal step length based on the alpha value
    double a = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double b = r[0]*(r1[0]+r2[0]) + r[1]*(r1[1]+r1[1]) + r[2]*(r1[2]+r2[2]);
    double c = log(0.5)/(2*alpha*omega);
    return (-b + sqrt(b*b - 4*a*c))/(2*a);
    /*
    double vecSum = r1.length() + r2.length();
    return (-(vecSum) + sqrt(pow(vecSum,2) - 2*log(0.5)/(alpha*omega)))/(2.0*r);
    */
}

double StepLength2(vec3 r1, vec3 r2, double alpha, double omega, double *r){
    // Calculates the optimal step length based on the alpha value
    double b = r[0]*(r1[0]+r2[0]) + r[1]*(r1[1]+r2[1]) + r[2]*(r[2]+r[2]);
    double a = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double c = log(0.5)/(2*alpha*omega);
    return (-b+sqrt(b*b - 4*a*c))/(2*a);
}

double LaplaceOperator(vec3 r1, vec3 r2, double alpha, double omega){
    // Function that calculates the Laplace operator numerically
    double SecondDerivative = 0;
    double dr = 1e-5;
    double wavefunc =  TrialWavefunc_T1(r1, r2, alpha, omega);
    for (int i=0; i<3; i++){
        vec3 rchange(0,0,0);
        rchange[i] = dr;
        SecondDerivative -= (TrialWavefunc_T1(r1+rchange, r2, alpha, omega) -\
                2*wavefunc + TrialWavefunc_T1(r1-rchange, r2, alpha, omega));
        SecondDerivative -= (TrialWavefunc_T1(r1, r2+rchange, alpha, omega) -\
                2*wavefunc + TrialWavefunc_T1(r1, r2-rchange, alpha, omega));
    }
    return 0.5*SecondDerivative/(wavefunc*(dr*dr));
}

double LaplaceOperator_T2(vec3 r1, vec3 r2, double alpha, double beta, double omega){
    // Function that calculates the Laplace operator numerically. Used for the second trial wave function
    double SecondDerivative = 0;
    double dr = 1e-5;
    double wavefunc =  TrialWavefunc_T2(r1, r2, alpha, beta, omega);
    for (int i=0; i<3; i++){
        vec3 rchange(0,0,0);
        rchange[i] = dr;
        SecondDerivative -= (TrialWavefunc_T2(r1+rchange, r2, alpha, beta, omega) -\
                2*wavefunc + TrialWavefunc_T2(r1-rchange, r2, alpha, beta, omega));
        SecondDerivative -= (TrialWavefunc_T2(r1, r2+rchange, alpha, beta, omega) -\
                2*wavefunc + TrialWavefunc_T2(r1, r2-rchange, alpha, beta, omega));
    }
    return 0.5*SecondDerivative/(wavefunc*(dr*dr));
}

double LaplaceAnalytic(vec3 r1, vec3 r2, double alpha, double omega){
    // Function that calculates the analytical expression of the Laplace operator
    double AlphaOmega = alpha*omega;
    return  -0.5*(pow(AlphaOmega,2)*r1.lengthSquared() - 3*AlphaOmega) \
            - 0.5*(pow(AlphaOmega,2)*r2.lengthSquared() - 3*AlphaOmega);
}

void Metropolis_method_T1(int MC_cycles, double omega, double alpha, double beta, double *ExpectationValues\
                       , int Analytic = 0){
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
    double r_12 = 0;
    double omega2 = omega*omega;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = StepLength(r1, r2, alpha, omega, rDistr);
    int counter = 0;
    // Runs the Monte Carlo cycles
    double OldWavefuncSquared = pow(TrialWavefunc_T1(r1, r2, alpha, omega), 2);
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
        NewWavefuncSquared = pow(TrialWavefunc_T1(r1_new, r2_new, alpha, omega), 2);
        if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        step_length = StepLength(r1, r2, alpha, omega, rDistr);
        if (Analytic == 1){
            E_local = LaplaceAnalytic(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        }
        else{
            E_local = LaplaceOperator(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        }

        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
        MeanDistance += (r1-r2).length();
    }

    // Adding the energies and mean distance in their arrays
    ExpectationValues[0] = EnergySum;
    ExpectationValues[1] = EnergySquaredSum;
    ExpectationValues[2] = MeanDistance;
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Energy = "<< ExpectationValues[0]/(MC_cycles) << endl;
    cout << "Variance = " << EnergySquaredSum/MC_cycles - EnergySum*EnergySum/MC_cycles/MC_cycles << endl;;
    cout << "Accepted configs = " << counter << endl;
}

void Metropolis_method_T2(int MC_cycles, double omega, double alpha, double beta, double *ExpectationValues\
                          , int IncludeCoulomb){
    /* Function that solves the Metropolis method for the second trial wave function.
     * Argument IncludeCoulomb
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
    double r_12 = 0;
    double omega2 = omega*omega;

    double *rDistr = new double[6];
    for (int i=0; i<6; i++){
        rDistr[i] = distr(generator);
    }
    double step_length = StepLength(r1, r2, alpha, omega, rDistr);

    int counter = 0;
    double OldWavefuncSquared = pow(TrialWavefunc_T2(r1, r2, alpha, beta, omega), 2);
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

        NewWavefuncSquared = pow(TrialWavefunc_T2(r1_new, r2_new, alpha,  beta, omega), 2);
        if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        step_length = StepLength(r1, r2, alpha, omega, rDistr);
        r_12 = (r1-r2).length();
        E_local = LaplaceOperator_T2(r1, r2, alpha, beta, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared()) \
                + (1.0/r_12)*IncludeCoulomb;
        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
        MeanDistance += r_12;
    }

    // Adding the energies and mean distance in their arrays
    ExpectationValues[0] = EnergySum;
    ExpectationValues[1] = EnergySquaredSum;
    ExpectationValues[2] = MeanDistance;
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "E numeric = "<< ExpectationValues[0]/(MC_cycles) << endl;
    cout << "Variance = " << EnergySquaredSum/MC_cycles - EnergySum*EnergySum/MC_cycles/MC_cycles << endl;;
    cout << "Accepted configs = " << counter << endl;
}

void write_file(int MC_cycles, double *H_expectation, double alpha, double omega){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Energy = H_expectation[0]/MC_cycles;
    double Variance = H_expectation[1]/MC_cycles - Energy*Energy;
    double MeanDistance = H_expectation[2]/MC_cycles;
    ofile_global << setw(15) << alpha;
    ofile_global << setw(15) << Energy;
    ofile_global << setw(15) << Variance;
    ofile_global << setw(15) << MeanDistance;
    ofile_global << setw(15) << omega << endl;
}

void initialize_outfile(){
    ofile_global << setw(15) << "Alpha";
    ofile_global << setw(15) << "Energy";
    ofile_global << setw(15) << "Variance";
    ofile_global << setw(15) << "Mean Distance";
    ofile_global << setw(15) << "Omega" << endl;
}

int main(int argc, char *argv[])
{
    double *ExpectationValues;
    ExpectationValues = new double [3];
    int MC_cycles = 100000;
    double omega, alpha, step_length, beta;
    omega = 1;
    alpha = 1;
    beta = 1;
    step_length = 1;
    Metropolis_method_T1(MC_cycles, omega, alpha, beta, ExpectationValues);

    /*
    ofile_global.open("Testing3.txt");
    initialize_outfile();
    for (alpha = 0.2; alpha<=2; alpha +=0.1){
        Metropolis_method_T1(MC_cycles, omega, alpha, beta, ExpectationValues);
        write_file(MC_cycles, ExpectationValues, alpha, omega);
    }
    ofile_global.close();
    */
    /*
    string filename = "Energy_Alpha_Mdistance_omega_";
    double omegas[] = {0.01, 0.5, 1};
    for (int i=0; i<3; i++){
        string fileout = filename;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile();
        for (alpha = 0.5; alpha<=2.3; alpha +=0.1){
            Metropolis_method_T1(MC_cycles, omegas[i], alpha, beta, ExpectationValues);
            write_file(MC_cycles, ExpectationValues, alpha, omegas[i]);
        }
        ofile_global.close();
    }
    */
    return 0;
}
