#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
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
    ofile_global << setw(15) << omega;
    ofile_global << setw(15) << ExpectationValues[3]/MC_cycles << endl;
}

void write_file_virial(int MC_cycles, double *ExpectationValues, double alpha, double beta, double omega){
    // Writes out data to the output file. Function made specific for the parallellization part.
    double Kinetic = ExpectationValues[0]/MC_cycles;
    double KineticVariance = ExpectationValues[1]/MC_cycles - Kinetic*Kinetic;
    double Potential = ExpectationValues[2]/MC_cycles;
    double PotentialVariance = ExpectationValues[3]/MC_cycles - Potential*Potential;
    //double Kinetic = ExpectationValues[0]/MC_cycles - ExpectationValues[2]/MC_cycles;
    //double KineticSquared = ExpectationValues[1]/MC_cycles - ExpectationValues[3]/MC_cycles \
                ;//- 2*Potential*Kinetic;
    //double KineticVariance = KineticSquared - Kinetic*Kinetic;

    double MeanDistance = ExpectationValues[4]/MC_cycles;
    ofile_global << setw(15) << alpha;
    ofile_global << setw(15) << beta;
    ofile_global << setw(15) << Kinetic;
    ofile_global << setw(15) << KineticVariance;
    ofile_global << setw(15) << Potential;
    ofile_global << setw(15) << PotentialVariance;
    ofile_global << setw(15) << MeanDistance;
    ofile_global << setw(15) << omega;
    ofile_global << setw(15) << ExpectationValues[5]/MC_cycles << endl;
}
void initialize_outfile(){
    ofile_global << setw(15) << "Alpha";
    ofile_global << setw(15) << "Beta";
    ofile_global << setw(15) << "Energy";
    ofile_global << setw(15) << "Variance";
    ofile_global << setw(15) << "Mean Distance";
    ofile_global << setw(15) << "Omega";
    ofile_global << setw(15) << " Accepted configs (%)" << endl;
}

void initialize_outfile_virial(){
    ofile_global << setw(15) << "Alpha";
    ofile_global << setw(15) << "Beta";
    ofile_global << setw(15) << "Kinetic";
    ofile_global << setw(15) << "Kin Variance";
    ofile_global << setw(15) << "Potential";
    ofile_global << setw(15) << "Pot Variance";
    ofile_global << setw(15) << "Mean Distance";
    ofile_global << setw(15) << "Omega ";
    ofile_global << setw(15) << " Accepted configs (%)" << endl;
}

void write_file_stability(int MC_cycles, double *energies, double *MCrange){
    for (int i = 0; i<= MC_cycles; i++){
        ofile_global << energies[i]  << endl;
    }
}

int main()
{
    clock_t start, finish;
    double *ExpectValues;
    ExpectValues = new double [4];
    int MC_cycles = 1000000;
    double omega, alpha, beta;
    omega = 1;
    alpha = 1;
    beta = 1;
    Wavefunctions FirstTrialFunc(0);
    Wavefunctions SecondTrialFunc(1);
    Metropolis_Quantum MSolver;

    // Testing the algorithm with laplace operators
    /*
    start = clock();
    MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alpha, omega, 0, 1);
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectValues[0]/(MC_cycles) << endl;
    cout << "Variance = "<< ExpectValues[1]/(MC_cycles) - ExpectValues[0]*ExpectValues[0]/MC_cycles/MC_cycles << endl;
    cout << "Accepted configs (percentage) = " << ExpectValues[3]/MC_cycles << endl;
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

    start = clock();
    MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alpha, omega, 0, 0);
    cout << "Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Kinetic numeric = "<< ExpectValues[0]/(MC_cycles) << endl;
    cout << "Variance = "<< ExpectValues[1]/(MC_cycles) - ExpectValues[0]*ExpectValues[0]/MC_cycles/MC_cycles << endl;
    cout << "Accepted configs (percentage) = " << ExpectValues[3]/MC_cycles << endl;
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    */
    //

    // Find optimal alpha
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
        for (double alphas = 0.5; alphas<=1.5; alphas +=0.02){
            MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alphas, omegas[i], 1);
            write_file(MC_cycles, ExpectValues, alphas, beta, omegas[i]);
        }
        ofile_global.close();
    }
    */
    // Testing stability 1
    /*
    string filename_stability = "Stability_test_MC_";
    for (int MCrun = 100; MCrun <= 1000001; MCrun *=10){
        string fileout = filename_stability;
        stringstream stream;
        stream << fixed  << MCrun;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);
        initialize_outfile();
        for (double alphas = 0.5; alphas<=1.5; alphas +=0.02){
            MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, alphas, 1, 1);
            write_file(MC_cycles, ExpectValues, alphas, beta, 1);
        }
        ofile_global.close();
    }

    */
    /*
    // Calculate mean distance with optimal alpha
    for (int i=0; i<3; i++){
       MSolver.Metropolis_T1(MC_cycles, FirstTrialFunc, ExpectValues, 1, omegas[i], 1);
       cout << "Omega = " << omegas[i] << endl;
       cout << "Mean distance = " << ExpectValues[2]/MC_cycles << endl;
    }
    */
    /*
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
                MSolver.Metropolis_T2(MC_cycles, SecondTrialFunc, ExpectValues, alphas, betas, omegas[i]);
                write_file(MC_cycles, ExpectValues, alphas, betas, omegas[i]);
            }
        }
        ofile_global.close();
    }
    */
    /*
    // Use optimal alpha, beta
    double AlphaOptimal[] = {0.58, 0.98, 0.98};
    double BetaOptimal[] = {0.2, 0.3, 0.52};
    ofile_global.open("Optimal_Energy_SecondTrialWaveFunc.txt");
    initialize_outfile();
    for (int i=0; i<3; i++){
        MSolver.Metropolis_T2(MC_cycles, SecondTrialFunc, ExpectValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
        write_file(MC_cycles, ExpectValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
    }
    ofile_global.close();

    ofile_global.open("Optimal_Energy_FirstTrialWaveFunc.txt");
    initialize_outfile();
    for (int i=0; i<3; i++){
        MSolver.Metropolis_T2(MC_cycles, FirstTrialFunc, ExpectValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
        write_file(MC_cycles, ExpectValues, AlphaOptimal[i], BetaOptimal[i], omegas[i]);
    }
    ofile_global.close();
    */
    // Virial test

    double AlphaOptimal[] = {0.6, 0.98, 0.98};
    double BetaOptimal[] = {0.2, 0.36, 0.52};
    double *VirialExpect = new double[6];
    // New omegas
    int N_omegas = 50;
    double *omegas2 = new double[N_omegas];
    double omega_step = (1.0-0.01)/(N_omegas-1);
    for (int i=0; i<N_omegas; i++){
        omegas2[i] = 0.01 + i*omega_step;
    }
    /*
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
    */
    // Testing stability 2, using virial problem
    string filename_virial_stability = "Virial_stability_MC_";
    for (int MCrun = 10000; MCrun <= MC_cycles; MCrun *= 10){
        string fileout = filename_virial_stability;
        stringstream stream;
        stream << fixed << MCrun;
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        ofile_global.open(fileout);

        initialize_outfile_virial();
        for (int i=0; i<N_omegas; i++){
            MSolver.Metropolis_Virial(MC_cycles, SecondTrialFunc, VirialExpect,\
                                      0.98, 0.48, omegas2[i], 1);
            write_file_virial(MC_cycles, VirialExpect, 0.98, 0.48, omegas2[i]);
        }
        ofile_global.close();
    }
    return 0;
}
