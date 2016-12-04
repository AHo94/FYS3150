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
//                Metropolis_method_T2(MC_cycles, omegas[i], alphas, betas, ExpectationValues);
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

int main()
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
    /*
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
    return 0;
}
