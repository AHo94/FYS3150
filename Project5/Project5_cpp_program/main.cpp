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
    // Function that calculates the trial wave function when r_12 -> 0
    return exp(-0.5*alpha*omega*(r1.lengthSquared() + r2.lengthSquared()));
}

double StepLength(vec3 r1, vec3 r2, double alpha, double omega, double r){
    double vecSum = r1.length() + r2.length();
    return (-(vecSum) + sqrt(pow(vecSum,2) - 2*log(0.5)*r/(alpha*omega)))/(2.0*r);
}

double LaplaceOperator(vec3 r1, vec3 r2, double alpha, double omega){
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

double LaplaceAnalytic(vec3 r1, vec3 r2, double alpha, double omega){
    double AlphaOmega = alpha*omega;
    return  -0.5*(pow(AlphaOmega,2)*r1.lengthSquared() - 3*AlphaOmega) \
            - 0.5*(pow(AlphaOmega,2)*r2.lengthSquared() - 3*AlphaOmega);
}

double LocalEnergy(vec3 r1, vec3 r2, double alpha, double omega){

}

void Metropolis_method(int MC_cycles, double omega, double alpha, double step_length, double *H_expect){
    // A function that uses the Metropolis method
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    // Random initial starting position for both electrons
    vec3 r1(distr(generator),distr(generator),distr(generator));
    vec3 r2(distr(generator),distr(generator),distr(generator));

    double E_local = 0;
    double EnergySum = 0;
    double EnergySquaredSum = 0;
    double NewWavefuncSquared = 0;
    double omega2 = omega*omega;
    double OldWavefuncSquared = pow(TrialWavefunc_T1(r1, r2, alpha, omega), 2);

    //double step_length = StepLength(r1, r2, alpha, omega, distr(generator));

    double E_an = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha*alpha) + 3*alpha*omega;
    int counter = 0;
    for (int cycle=0; cycle<MC_cycles; cycle++){
        vec3 r1_new(0,0,0);
        vec3 r2_new(0,0,0);
        for (int j=0; j<3; j++){
             r1_new[j] = r1[j] + step_length*distr(generator);
             r2_new[j] = r2[j] + step_length*distr(generator);
        }
        NewWavefuncSquared = pow(TrialWavefunc_T1(r1_new, r2_new, alpha, omega), 2);
        if (distr(generator) <= NewWavefuncSquared/OldWavefuncSquared){
            r1 = r1_new;
            r2 = r2_new;
            OldWavefuncSquared = NewWavefuncSquared;
            counter += 1;
        }
        //step_length = StepLength(r1, r2, alpha, omega, distr(generator));
        //E_local = LaplaceAnalytic(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        E_local = LaplaceOperator(r1, r2, alpha, omega) + 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared());
        EnergySum += E_local;
        EnergySquaredSum += E_local*E_local;
    }
    H_expect[0] = EnergySum;
    H_expect[1] = EnergySquaredSum;
    cout << "E = "<< H_expect[0]/(MC_cycles) << endl;
    cout << "E analytic = " << E_an << endl;
    cout << "Variance = " << EnergySquaredSum/MC_cycles - EnergySum*EnergySum/MC_cycles/MC_cycles << endl;;
    cout << "Accepted configs = " << counter << endl;
}

void write_file(int MC_cycles, double *H_expectation, double alpha){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Energy = H_expectation[0]/MC_cycles;
    ofile_global << setw(15) << alpha ;
    ofile_global << setw(15) << Energy << endl;;
}

void initialize_outfile(){
    ofile_global << setw(15) << "alpha";
    ofile_global << setw(15) << "Energy" << endl;
}

int main(int argc, char *argv[])
{
    double *H_expect;
    H_expect = new double [2];
    int MC_cycles = 100000;
    double omega, alpha, step_length;
    omega = 1;
    //alpha = 1;
    step_length = 0.005;
    ofile_global.open("Energy_alpha.txt");
    initialize_outfile();
    for (alpha = 0.2; alpha<=10; alpha +=0.1){
        Metropolis_method(MC_cycles, omega, alpha, step_length, H_expect);
        write_file(MC_cycles, H_expect, alpha);
    }
    ofile_global.close();
    return 0;
}
