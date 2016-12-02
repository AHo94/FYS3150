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
#include "metropolis_quantum.h"
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

double StepLength(vec3 r1, vec3 r2, double alpha, double omega, double *s){
    // Calculates the optimal step length based on the alpha value
    double a = 0;
    for (int i = 0; i<6; i++){
        a += s[i]*s[i];
    }
    double b = 2*(r1[0]*s[0] + r1[1]*s[1] + r1[2]*s[2] + r2[0]*s[3] + r2[1]*s[4] + r2[2]*s[5]);
    double c = log(0.5)/(alpha*omega);
    return (-b + sqrt(b*b - 4*a*c))/(2*a);
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

void Metropolis_method_T1(int MC_cycles, double omega, double alpha, double *ExpectationValues\
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
            E_local = LaplaceOperator(r1, r2, alpha, omega) + \
                    0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
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
    cout << "Variance = " << EnergySquaredSum/MC_cycles - EnergySum*EnergySum/MC_cycles/MC_cycles << endl;
    cout << "Accepted configs (percentage) = " <<(double)counter/MC_cycles << endl;
    cout << "\n";
}

void Metropolis_method_T2(int MC_cycles, double omega, double alpha, double beta, double *ExpectationValues){
    /* Function that solves the Metropolis method for the second trial wave function.
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
        //E_local = LaplaceOperator_T2(r1, r2, alpha, beta, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared()) \
                + (1.0/r_12);
        E_local = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha*alpha) + 3*alpha*omega + 1/r_12\
                + (1/(2*(1+beta*r_12)*(1+beta*r_12)))*(alpha*omega*r_12 - 1/(2*(1+beta*r_12)*(1+beta*r_12))\
                - 2/r_12 + 2*beta/(1+beta*r_12));

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
    cout << "Accepted configs (percentage) = " << (double) counter/MC_cycles << endl;
}

void Metropolis_Virial(int MC_cycles, double alpha, double beta, double omega, double *ExpectationValues,
                       int CoulombInt=1){
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

    double E_kin = 0;
    double E_pot = 0;
    double KineticSum = 0;
    double KineticSquaredSum = 0;
    double PotentialSum = 0;
    double PotentialSquaredSum = 0;
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
        E_kin = LaplaceOperator_T2(r1, r2, alpha, beta, omega);
        E_pot = + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared()) + CoulombInt*(1.0/r_12);

        KineticSum += E_kin;
        KineticSquaredSum += E_kin*E_kin;
        PotentialSum += E_pot;
        PotentialSquaredSum += E_pot*E_pot;
    }

    // Adding the energies and mean distance in their arrays
    ExpectationValues[0] = KineticSum;
    ExpectationValues[1] = KineticSquaredSum;
    ExpectationValues[2] = PotentialSum;
    ExpectationValues[3] = PotentialSquaredSum;
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectationValues[0]/(MC_cycles) << endl;
    cout << "Potential numeric = "<< ExpectationValues[2]/(MC_cycles) << endl;
    cout << "Accepted configs (percentage) = " <<(double) counter/MC_cycles << endl;
    cout << "\n" << endl;
}


void write_file(int MC_cycles, double *ExpectationValues, double alpha, double beta, double omega){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Energy = ExpectationValues[0]/MC_cycles;
    double Variance = ExpectationValues[1]/MC_cycles - Energy*Energy;
    double MeanDistance = ExpectationValues[2]/MC_cycles;
    ofile_global << setw(15) << alpha;
    ofile_global << setw(15) << beta;
    ofile_global << setw(15) << Energy;
    ofile_global << setw(15) << Variance;
    ofile_global << setw(15) << MeanDistance;
    ofile_global << setw(15) << omega << endl;
}

void initialize_outfile(){
    ofile_global << setw(15) << "Alpha";
    ofile_global << setw(15) << "Beta";
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
    double omega, alpha, beta;
    omega = 1;
    alpha = 1;
    beta = 1;
    Metropolis_method_T1(MC_cycles, omega, alpha, ExpectationValues);
    Metropolis_method_T1(MC_cycles, omega, alpha, ExpectationValues, 1);
    //Metropolis_method_T2(MC_cycles, omega, alpha, beta, ExpectationValues);


    string filename = "Energy_Alpha_Mdistance_omega_";
    double omegas[] = {0.01, 0.5, 1};
    /*
    for (int i=0; i<3; i++){
        string fileout = filename;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile();
        for (double alphas = 0.5; alphas<=2.3; alphas +=0.1){
            Metropolis_method_T1(MC_cycles, omegas[i], alphas, ExpectationValues);
            write_file(MC_cycles, ExpectationValues, alphas, beta, omegas[i]);
        }
        ofile_global.close();
    }
    */
/*
    ofile_global.open("Beta_test2.txt");
    initialize_outfile();
    for (double alphas = 0.8; alphas <= 1.3; alphas += 0.02){
        for (double betas = 0.2; betas <= 1.3; betas += 0.02){
            Metropolis_method_T2(MC_cycles, omega, alphas, betas, ExpectationValues);
            write_file(MC_cycles, ExpectationValues, alphas, betas, omega);
        }
    }
    ofile_global.close();
    */

    double *VirialExpect = new double[4];
    double AlphaOptimal = 1.06;
    double BetaOptimal = 0.56;
    Metropolis_Virial(MC_cycles, AlphaOptimal, BetaOptimal, omega, VirialExpect, 0);
    /*
    string filename2 = "Virial_omega_";
    for (int i=0; i<3; i++){
        string fileout2 = filename2;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout2.append(argument);
        fileout2.append(".txt");
        ofile_global.open(fileout2);

        Metropolis_Virial(MC_cycles, AlphaOptimal, BetaOptimal, omegas[i], VirialExpect, 0);
        write_file(MC_cycles, VirialExpect, AlphaOptimal, BetaOptimal, omegas[i]);
    }
    */
    cout << "CLASS TEST" << endl;
    Metropolis_Quantum MSolver;
    MSolver.Metropolis_T1(MC_cycles, alpha, omega);
    MSolver.Metropolis_T1(MC_cycles, alpha, omega, 1);
    return 0;
}
