#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
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

double StepLength(vec3 r1, vec3 r2, double alpha, double omega, double r){
    // Calculates the optimal step length based on the alpha value
    double vecSum = r1.length() + r2.length();
    return (-(vecSum) + sqrt(pow(vecSum,2) - 2*log(0.5)/(alpha*omega)))/(2.0*r);
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

void Metropolis_method(int MC_cycles, double omega, double alpha, double beta, double *H_expect, int TrialFunc){
    // A function that uses the Metropolis method
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

    double step_length = StepLength(r1, r2, alpha, omega, distr(generator));
    double r_12 = (r1-r2).length();
    double CoulombInt = 0;
    if (TrialFunc == 2){
        // Adds Coulomb interaction for second trial wave function
        CoulombInt = 1.0/r_12;
    }
    // Analytic expressions as a comparison
    double E_an = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha*alpha) + 3*alpha*omega + CoulombInt;
    double E_an_T2 = E_an + 1.0/(2*pow(1+beta*r_12, 2))*(alpha*omega*r_12 - 1.0/(2*pow(1+beta*r_12, 2)) - 2.0/r_12 \
                                                         + 2.0*beta/(1+beta*r_12));
    int counter = 0;

    if (TrialFunc == 1){
        double OldWavefuncSquared = pow(TrialWavefunc_T1(r1, r2, alpha, omega), 2);
        for (int cycle=0; cycle<MC_cycles; cycle++){
            // Running Monte Carlo cycles
            vec3 r1_new(0,0,0);
            vec3 r2_new(0,0,0);
            double s = distr(generator);
            /*
            for (int j=0; j<3; j++){
                r1_new[j] = r1[j] +s;// step_length*distr(generator);
                r2_new[j] = r2[j] +s;// step_length*distr(generator);
            }
            */
            r1_new = r1 + s;
            r2_new = r2 + s;
            NewWavefuncSquared = pow(TrialWavefunc_T1(r1_new, r2_new, alpha, omega), 2);
            if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
                r1 = r1_new;
                r2 = r2_new;
                OldWavefuncSquared = NewWavefuncSquared;
                counter += 1;
            }
            step_length = StepLength(r1, r2, alpha, omega, s);
            //E_local = LaplaceAnalytic(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
            E_local = LaplaceOperator(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
            EnergySum += E_local;
            EnergySquaredSum += E_local*E_local;
            MeanDistance += (r1-r2).length();
        }
    }
    else if(TrialFunc == 2){
        double OldWavefuncSquared = pow(TrialWavefunc_T2(r1, r2, alpha, beta, omega), 2);
        for (int cycle=0; cycle<MC_cycles; cycle++){
            // Running Monte Carlo cycles
            vec3 r1_new(0,0,0);
            vec3 r2_new(0,0,0);
            double s = distr(generator);
            /*
            for (int j=0; j<3; j++){
                 r1_new[j] = r1[j] +s;// step_length*distr(generator);
                 r2_new[j] = r2[j] +s;// step_length*distr(generator);
            }
            */
            r1_new = r1 + s;
            r2_new = r2 + s;
            NewWavefuncSquared = pow(TrialWavefunc_T2(r1_new, r2_new, alpha,  beta, omega), 2);
            if (distrUniform(generator) <= NewWavefuncSquared/OldWavefuncSquared){
                r1 = r1_new;
                r2 = r2_new;
                OldWavefuncSquared = NewWavefuncSquared;
                counter += 1;
            }
            step_length = StepLength(r1, r2, alpha, omega, s);
            r_12 = (r1-r2).length();
            E_local = LaplaceOperator_T2(r1, r2, alpha, beta, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared()) \
                    + 1.0/r_12;
            EnergySum += E_local;
            EnergySquaredSum += E_local*E_local;
            MeanDistance += r_12;
        }
    }
    else{
        cout << "TrialFunc = " << TrialFunc << endl;
        cout << "Use TrialFunc = 1 or TrialFunc = 2. Stopping program" << endl;
        exit(1);
    }
    H_expect[0] = EnergySum;
    H_expect[1] = EnergySquaredSum;
    H_expect[2] = MeanDistance;
    cout << "E numeric = "<< H_expect[0]/(MC_cycles) << endl;
    if (TrialFunc == 1){cout << "E analytic = " << E_an << endl;}
    else if (TrialFunc == 2){cout << "E analytic = " << E_an_T2 << endl;}
    cout << "Variance = " << EnergySquaredSum/MC_cycles - EnergySum*EnergySum/MC_cycles/MC_cycles << endl;;
    cout << "Accepted configs = " << counter << endl;
}

void write_file(int MC_cycles, double *H_expectation, double alpha, double omega){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Energy = H_expectation[0]/MC_cycles;
    double MeanDistance = H_expectation[2]/MC_cycles;
    ofile_global << setw(15) << alpha;
    ofile_global << setw(15) << Energy;
    ofile_global << setw(15) << MeanDistance;
    ofile_global << setw(15) << omega << endl;
}

void initialize_outfile(){
    ofile_global << setw(15) << "alpha";
    ofile_global << setw(15) << "Energy";
    ofile_global << setw(15) << "Mean Distance";
    ofile_global << setw(15) << "Omega" << endl;
}

int main(int argc, char *argv[])
{
    double *H_expect;
    H_expect = new double [3];
    int MC_cycles = 100000;
    double omega, alpha, step_length, beta;
    omega = 1;
    alpha = 1;
    beta = 1;
    step_length = 1;
    Metropolis_method(MC_cycles, omega, alpha, beta, H_expect, 2);
    /*
    ofile_global.open("Energy_Alpha_Mdistance.txt");
    initialize_outfile();
    for (alpha = 0; alpha<=2; alpha +=0.1){
        Metropolis_method(MC_cycles, omega, alpha, beta, step_length, H_expect, 1);
        write_file(MC_cycles, H_expect, alpha, omega);
    }
    ofile_global.close();
*/

    return 0;
}
