#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "vec3.h"
#include "metropolis_quantum.h"
#include "wavefunctions.h"
using namespace std;

ofstream ofile_global;

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
        //E_kin = LaplaceOperator_T2(r1, r2, alpha, beta, omega);
        E_tot = 0.5*omega2*(r1.lengthSquared() + r2.lengthSquared())*(1-alpha*alpha) + 3*alpha*omega \
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
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectationValues[0]/(MC_cycles) << endl;
    cout << "Potential numeric = "<< ExpectationValues[2]/(MC_cycles) << endl;
    cout << "Accepted configs (percentage) = " <<(double) counter/MC_cycles << endl;
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

void write_file_virial(int MC_cycles, double *ExpectationValues, double alpha, double beta, double omega){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Kinetic = ExpectationValues[0]/MC_cycles;
    double KineticVariance = ExpectationValues[1]/MC_cycles - Kinetic*Kinetic;
    double Potential = ExpectationValues[2]/MC_cycles;
    double PotentialVariance = ExpectationValues[3]/MC_cycles - Potential*Potential;
    double MeanDistance = ExpectationValues[4]/MC_cycles;
    ofile_global << setw(15) << alpha;
    ofile_global << setw(15) << beta;
    ofile_global << setw(15) << Kinetic;
    ofile_global << setw(15) << KineticVariance;
    ofile_global << setw(15) << Potential;
    ofile_global << setw(15) << PotentialVariance;
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

void initialize_outfile_virial(){
    ofile_global << setw(15) << "Alpha";
    ofile_global << setw(15) << "Beta";
    ofile_global << setw(15) << "Kinetic";
    ofile_global << setw(15) << "Kin Variance";
    ofile_global << setw(15) << "Potential";
    ofile_global << setw(15) << "Pot Variance";
    ofile_global << setw(15) << "Mean Distance";
    ofile_global << setw(15) << "Omega" << endl;
}

void Find_Optimal_AlphaBeta(int MC_cycles, int N_omegas, double *omegas, double *OptimalAlphas, double *OptimalBetas){
    double *ExpectationValues;
    ExpectationValues = new double [3];

    double alpha_max = 1.0;
    double beta_max = 0.7;
    double alpha_min = 0.5;
    double beta_min = 0.2;
    double StepSize = 0.02;

    int N_size = (int)((alpha_max - alpha_min)/(StepSize))*((beta_max-beta_min)/(StepSize));

    double *MinEnergies = new double[N_size+1];
    double *AlphasCalc = new double[N_size+1];
    double *BetasCalc = new double[N_size+1];
    int indexMin = 0;
    double minimum = 0;
    int indexCounter = 0;
    for (int i=0; i<N_omegas; i++){
        for (double alphas = 0.5; alphas <= 1.0; alphas += 0.02){
            for (double betas = 0.2; betas <= 0.7; betas += 0.02){
                Metropolis_method_T2(MC_cycles, omegas[i], alphas, betas, ExpectationValues);
                MinEnergies[indexCounter] = ExpectationValues[0]/MC_cycles;
                AlphasCalc[indexCounter] = alphas;
                BetasCalc[indexCounter] = betas;
                indexCounter += 1;
            }
        }
        indexCounter = 0;
        minimum = MinEnergies[0];
        for (int j=1; j<N_size; j++){
            if (MinEnergies[j] < minimum){
                minimum = MinEnergies[j];
                indexMin = j;
            }
        }
        OptimalAlphas[i] = AlphasCalc[indexMin];
        OptimalBetas[i] = BetasCalc[indexMin];
        indexMin = 0;
    }
}

int main(int argc, char *argv[])
{
    double *ExpectValues;
    ExpectValues = new double [4];
    int MC_cycles = 100000;
    double omega, alpha, beta;
    omega = 1;
    alpha = 1;
    beta = 1;

    cout << "CLASS TEST" << endl;
    Wavefunctions FirstTrialFunc(0);
    Wavefunctions SecondTrialFunc(1);
    Metropolis_Quantum MSolver;

    // Testing the algorithm
    MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alpha, omega, 0, 1);
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectValues[0]/(MC_cycles) << endl;
    cout << "Variance = "<< ExpectValues[1]/(MC_cycles) - ExpectValues[0]*ExpectValues[0]/MC_cycles/MC_cycles << endl;
    cout << "Accepted configs (percentage) = " << ExpectValues[3]/MC_cycles << endl;

    MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alpha, omega, 0, 0);
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectValues[0]/(MC_cycles) << endl;
    cout << "Variance = "<< ExpectValues[1]/(MC_cycles) - ExpectValues[0]*ExpectValues[0]/MC_cycles/MC_cycles << endl;
    cout << "Accepted configs (percentage) = " << ExpectValues[3]/MC_cycles << endl;
    /*
    // Find optimal alpha
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
        for (double alphas = 0.5; alphas<=1.5; alphas +=0.02){
            MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectationValues, alphas, omegas[i], 0, 1);
            write_file(MC_cycles, ExpectationValues, alphas, beta, omegas[i]);
        }
        ofile_global.close();
    }

    // Find optimal alpha, beta
    string filename_optimal = "Optimal_AlphaBeta_omega_";
    for (int i=0; i<3; i++){
        string fileout = filename_optimal;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);
        initialize_outfile();
        for (double alphas = 0.5; alphas <= 1.0; alphas += 0.02){
            for (double betas = 0.2; betas <= 0.7; betas += 0.02){
                MSolver.Metropolis_T2(MC_cycles, SecondTrialFunc, ExpectationValues, alphas, betas, omegas[i]);
                write_file(MC_cycles, ExpectationValues, alphas, betas, omegas[i]);
            }
        }
        ofile_global.close();
    }

    // Use optimal alpha, beta
    double AlphaOptimal[] = {0.58, 0.98, 0.98};
    double BetaOptimal[] = {0.2, 0.3, 0.52};
    ofile_global.open("Optimal_Energy_SecondTrialWaveFunc.txt");
    initialize_outfile();
    for (int i=0; i<3; i++){
        MSolver.Metropolis_T2(MC_cycles, SecondTrialFunc, ExpectationValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
        write_file(MC_cycles, ExpectationValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
    }
    ofile_global.close();

    ofile_global.open("Optimal_Energy_FirstTrialWaveFunc.txt");
    initialize_outfile();
    for (int i=0; i<3; i++){
        MSolver.Metropolis_T2(MC_cycles, FirstTrialFunc, ExpectationValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
        write_file(MC_cycles, ExpectationValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
    }
    ofile_global.close();
    */
    // Virial test

    double AlphaOptimal[] = {0.58, 0.98, 0.98};
    double BetaOptimal[] = {0.2, 0.36, 0.48};
    double *VirialExpect = new double[6];
    // New omegas
    int N_omegas = 50;
    double *omegas2 = new double[N_omegas];
    double omega_step = (1.0-0.01)/(N_omegas-1);
    for (int i=0; i<N_omegas; i++){
        omegas2[i] = 0.01 + i*omega_step;
    }

    string filename2 = "Virial_data_V";
    for (int j=0; j<3; j++){
        string fileout = filename2;
        stringstream stream;
        stream << fixed << j;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile_virial();
        for (int i=0; i<N_omegas; i++){
            MSolver.Metropolis_Virial(MC_cycles, SecondTrialFunc, VirialExpect,\
                                      AlphaOptimal[j], BetaOptimal[j], omegas2[i], 1);
            write_file_virial(MC_cycles, VirialExpect, AlphaOptimal[j], BetaOptimal[j], omegas2[i]);
        }
        ofile_global.close();
    }

    string filename3 = "Virial_NoCoulomb_data_V";
    for (int j=0; j<3; j++){
        string fileout = filename3;
        stringstream stream;
        stream << fixed << j;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile_virial();
        for (int i=0; i<N_omegas; i++){
            MSolver.Metropolis_Virial(MC_cycles, SecondTrialFunc, VirialExpect,\
                                      AlphaOptimal[j], BetaOptimal[j], omegas2[i], 0);
            write_file_virial(MC_cycles, VirialExpect, AlphaOptimal[j], BetaOptimal[j], omegas2[i]);
        }

        ofile_global.close();
    }

    /*
    string filename2 = "Virial_data_V";
    //ofile_global.open(filename2);
    for (int j=0; j<3; j++){

        string fileout = filename2;
        stringstream stream;
        stream << fixed << j;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile_virial();
        for (int i=0; i<N_omegas; i++){
            Metropolis_Virial(MC_cycles, AlphaOptimal[j], BetaOptimal[j], omegas2[i], VirialExpect, 1);
            write_file_virial(MC_cycles, VirialExpect, AlphaOptimal[j], BetaOptimal[j], omegas2[i]);
        }
        ofile_global.close();
    }
    string filename3 = "Virial_NoCoulomb_data_V";
    //ofile_global.open(filename3);
    for (int j=0; j<3; j++){
        string fileout = filename3;
        stringstream stream;
        stream << fixed << j;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile_virial();
        for (int i=0; i<N_omegas; i++){
            Metropolis_Virial(MC_cycles, AlphaOptimal[j], BetaOptimal[j], omegas2[i], VirialExpect, 0);
            write_file_virial(MC_cycles, VirialExpect, AlphaOptimal[j], BetaOptimal[j], omegas2[i]);
        }

        ofile_global.close();
    }
    */
    return 0;
}
