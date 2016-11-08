#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

using namespace std;

inline int periodic(int i, int limit, int add){
    // Funcion that ensures periodic boundary conditions
    return (i+limit+add) % (limit);
}

void initialize_system(int L, double **Spin_matrix, int random_state = 0){
    /* Initializes the lattice by setting spin states.
     * The argument random_state acts as an optional argument. Set to zero by default.
     * If random_state = 0, we set up a random state.
     * If random_state = 1, all spins will be up-spins.
     * If random_state = -1, all spins will be down-spins.
     * Else, print out an error.
    */
    std::random_device rd; // obtain a random number from hardware
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distr(0,1); // define the range

    if (random_state == 0){
        // Sets up a random microstate
        for (int x=0; x < L; x++){
            for (int y=0; y<L; y++){
                int spin = distr(generator);
                if (spin == 0){
                    Spin_matrix[x][y] = -1.0;
                }
                else{
                    Spin_matrix[x][y] = 1.0;
                }
            }
        }
    }
    else if(random_state == 1){
        // All spins will be up-spins
        for (int x=0; x < L; x++){
            for (int y=0; y<L; y++){
                Spin_matrix[x][y] = 1;
            }
        }
    }
    else if (random_state == -1){
        // All spins will be down-spins
        for (int x=0; x < L; x++){
            for (int y=0; y<L; y++){
                Spin_matrix[x][y] = -1;
            }
        }
    }
    else{
        cout << "Argument random_state not set correctly. By default, random_state = 0."
                " Try random_state = 1 or random_state = -1" << endl;
        cout << "Current input for random_state = " << random_state << endl;
        exit(1);
    }
}

double Calculate_E(double **Spin_matrix, int L){
    // Calculates the energy
    double Energy = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            Energy -= Spin_matrix[i][j]*(Spin_matrix[periodic(i, L, -1)][j] +
                                         Spin_matrix[i][periodic(j, L, -1)]);
        }
    }
    return Energy;
}

double Calculate_M(double **Spin_matrix, int L){
    // Calculates the magnetic moment
    double Magnetic_moment = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            Magnetic_moment += Spin_matrix[i][j];
        }
    }
    return Magnetic_moment;
}

void Metropolis_method(int L, int MC_cycles, double Temperature, double *Expectation_values, int &accepted_flip,
                       double *Mean_Energies, double *Mean_Mag_moments, int random_state = 0){
    // A function that uses the Metropolis method
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    double **Spin_matrix;
    Spin_matrix = new double*[L];
    for (int i=0; i<L; i++){
        Spin_matrix[i] = new double[L];
    }
    double EnergySum=0;
    double MSum=0;
    double EnergySquaredSum=0;
    double MSquaredSum=0;
    double fabsMSum = 0;
    initialize_system(L, Spin_matrix, random_state);

    double currentEnergy = Calculate_E(Spin_matrix, L);
    double currentM = Calculate_M(Spin_matrix, L);
    for (int cycle=0; cycle<MC_cycles; cycle++){
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                int ix = distr(generator)*L;
                int iy = distr(generator)*L;
                Spin_matrix[ix][iy] *= -1.0;    // Flip a random spin
                double E_t = Calculate_E(Spin_matrix, L);   // Calculate new energy after spin flip
                double dE = E_t - currentEnergy;    // Calculate energy difference
                /*
                int dE =  2*Spin_matrix[ix][iy]*
                      (Spin_matrix[ix][periodic(iy,L,-1)]+
                       Spin_matrix[periodic(ix,L,-1)][iy] +
                       Spin_matrix[ix][periodic(iy,L,1)] +
                       Spin_matrix[periodic(ix,L,1)][iy]);
                if (distr(generator) <= exp(-dE/Temperature)){
                    Spin_matrix[ix][iy] *= -1.0;
                    currentM += 2*Spin_matrix[ix][iy];//Calculate_M(Spin_matrix, L);
                    currentEnergy += dE;
                }
                */

                if (dE <= 0){
                    // Accept new state, use the energy for this state
                    currentEnergy = Calculate_E(Spin_matrix, L);
                    currentM = Calculate_M(Spin_matrix, L);
                    accepted_flip += 1;
                } else {
                    if (distr(generator) <= exp(-dE/Temperature)){
                        // Accept new state, use new energy
                        currentEnergy = Calculate_E(Spin_matrix, L);
                        currentM = Calculate_M(Spin_matrix, L);
                        accepted_flip += 1;
                    } else {
                        // Do not accept new state. Flip back spin and use old energy
                        Spin_matrix[ix][iy] *= -1.0;
                    }
                }
            }
        }
        // Sums up the samplings
        EnergySum += currentEnergy;
        EnergySquaredSum += currentEnergy*currentEnergy;
        MSum += currentM;
        fabsMSum += fabs(currentM);
        MSquaredSum += currentM*currentM;

        Mean_Energies[cycle] = currentEnergy;
        Mean_Mag_moments[cycle] = fabs(currentM);
    }
    Expectation_values[0] = EnergySum;
    Expectation_values[1] = EnergySquaredSum;
    Expectation_values[2] = MSum;
    Expectation_values[3] = MSquaredSum;
    Expectation_values[4] = fabsMSum;
}

//void write_file(int L, int MC_cycles, double Temperature, double *Expectation_values){
    /* Function that writes data to the output file.
     * All expectation values are normalized with respect to the number of Monte Carlo cycles,
     * as well as the number of spins.
    */
    /*
    double E_expectation = Expectation_values[0]/MC_cycles;
    double E2_expectation = Expectation_values[1]/MC_cycles;
    double M_expectation = Expectation_values[2]/MC_cycles;
    double M2_expectation = Expectation_values[3]/MC_cycles;
    double Mabs_expectation = Expectation_values[4]/MC_cycles;

    double E_variance = (E2_expectation - E_expectation*E_expectation)/L/L;
    double M_variance = (M2_expectation - Mabs_expectation*Mabs_expectation)/L/L;
    double C_v = E_variance/Temperature/Temperature;
    double chi = M_variance/Temperature/Temperature;

    ofile << setw(15) << L;
    ofile << setw(15) << MC_cycles;
    ofile << setw(15) << setprecision(5) << Temperature;
    ofile << setw(15) << setprecision(5) << E_expectation/L/L;
    ofile << setw(15) << setprecision(5) << C_v;
    ofile << setw(15) << setprecision(5) << Mabs_expectation/L/L;
    ofile << setw(15) << setprecision(5) << chi;
    ofile << "\n";
    */
//}
/*
void Run_simulation(int L, double Temperature, int MC_max){
    clock_t start, finish;
    double *Expectation_values;
    cout << "Running for temperature: " << Temperature << endl;
    for (int MC_cycles = 100; MC_cycles <= MC_max; MC_cycles *= 2){
        start = clock();
        cout << "Using MC_cycles = " << MC_cycles << endl;
        Expectation_values = new double[5];
        //Metropolis_method(L, MC_cycles, Temperature, Expectation_values, 1);
        write_file(L, MC_cycles, Temperature, Expectation_values);
        finish = clock();
        cout << "Time elapsed for MC_cycles = " << MC_cycles << ":  " <<
                ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    }
}
*/
void set_initial_text(){
    /* Function that sets up an initial text for the output file.
     * Makes it easier to read off in the text file.
    */
    /*
    ofile << setw(15) << "# Spins" << setw(15) << "MC cycles" << setw(15) << "Temperature" << setw(15) << "<E>"
          << setw(15) << "C_v" << setw(15) << "<|M|>" << setw(15) << "Chi" << "\n";
    */
}

void write_file(int L, int T, int MC_cycles, int accepted_flip, double *Mean_energies, double *Mean_mag_moments,
                string filename_E, string filename_M){
    cout << "Saving data to file" << endl;
    ofstream ofileE, ofileM; // Output file
    ofileE.open(filename_E);
    ofileE << "<E>" << setw(15) << "MC Cycles" << setw(15);
    ofileE << "# Spins" << setw(15) << "Temperature" << setw(15) << "Accepted flips" << "\n";
    ofileE << setw(15) << MC_cycles << setw(15) << L << setw(15) << T << setw(15) << accepted_flip << "\n";

    ofileM.open(filename_M);
    ofileM << "<|M|>" << setw(15) << "MC Cycles" << setw(15);
    ofileM << "# Spins" << setw(15) << "Temperature" << "\n";
    ofileM << setw(15) << MC_cycles << setw(15) << L << setw(15) << T << "\n";
    int counter = 100;
    for (int i=0; i<MC_cycles; i++){
        if (counter == 100){
            ofileE << Mean_energies[i] << "\n";
            ofileM << Mean_mag_moments[i] << "\n";
            counter = 0;
        }
        counter += 1;
    }

    ofileE.close();
    ofileM.close();
}

int main()
{
    clock_t start, finish;
    double *Expectation_values;
    Expectation_values = new double[5];
    double T_init = 1.0;     // Temperature = 1.0 kT/J
    int L = 2;  // Number of spins
    int accepted_flip = 0;

    // Analytical expressions
    double AC_v = 64.0*(1+3*cosh(8.0/T_init))/(T_init*pow((cosh(8.0/T_init)+3), 2));
    double Achi = 8*(exp(8.0/T_init) + cosh(8.0/T_init) + 3.0/2.0)/(T_init*pow((cosh(8.0/T_init)+3), 2));

    int MC_cycles = 1000;
    double *Mean_energies, *Mean_mag_moments;
    Mean_energies = new double[MC_cycles];
    Mean_mag_moments = new double [MC_cycles];
    Metropolis_method(L, MC_cycles, T_init, Expectation_values, accepted_flip, Mean_energies, Mean_mag_moments);
    double C_v = (Expectation_values[1]/MC_cycles -
            Expectation_values[0]*Expectation_values[0]/MC_cycles/MC_cycles)/T_init/T_init;
    double Chi = (Expectation_values[3]/MC_cycles -
            Expectation_values[4]*Expectation_values[4]/MC_cycles/MC_cycles)/T_init/T_init;
    cout << "Number of Monte Carlo cycles = " << MC_cycles << endl;
    cout << "Analytic C_v = " << AC_v << ", Numerical C_v = " << C_v << endl;
    cout << "Analytic Chi = " << Achi << ", Numerical Chi = " << Chi << endl;
    cout << "\n" << "Running multiple times, using MC_cycles = "<< MC_cycles << endl;
    for (int i=0; i<=5; i++)
    {
        // Runs this example multiple times to showcase the stability of the algorithm
        Expectation_values = new double[5];
        Metropolis_method(L, MC_cycles, T_init, Expectation_values, accepted_flip, Mean_energies, Mean_mag_moments);
        double C_v = (Expectation_values[1]/MC_cycles -
                Expectation_values[0]*Expectation_values[0]/MC_cycles/MC_cycles)/T_init/T_init;
        double Chi = (Expectation_values[3]/MC_cycles -
                Expectation_values[4]*Expectation_values[4]/MC_cycles/MC_cycles)/T_init/T_init;
        cout << "Analytic C_v = " << AC_v << ", Numerical C_v = " << C_v << endl;
        cout << "Analytic Chi = " << Achi << ", Numerical Chi = " << Chi << endl;
        cout << " " << endl;
    }

    // 4c) Let now L = 20
    L = 20;
    cout << "Running L = 20. NOTE: This may take a while" << endl;
    MC_cycles = 100000;
    double T_final = 2.4;
    string filename_E = "Mean_E_T";
    string filename_M = "Mean_M_T";
    for (double Temperature = T_init; Temperature <= T_final; Temperature += 1.4){
        cout << "Running for T = " << Temperature << endl;
        start = clock();
        accepted_flip = 0;
        Expectation_values = new double[5];
        Mean_energies = new double[MC_cycles];
        Mean_mag_moments = new double [MC_cycles];
        Metropolis_method(L, MC_cycles, Temperature, Expectation_values, accepted_flip, Mean_energies, Mean_mag_moments);
        finish = clock();
        cout << "Time elapsed for MC_cycles = " << MC_cycles << ":  " <<
                ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

        string fileout_E = filename_E;
        string fileout_M = filename_M;
        stringstream stream;
        stream << fixed << setprecision(2) << Temperature;
        string argument = stream.str();
        fileout_E.append(argument);
        fileout_E.append(".txt");
        fileout_M.append(argument);
        fileout_M.append(".txt");
        write_file(L, Temperature, MC_cycles, accepted_flip, Mean_energies, Mean_mag_moments, fileout_E, fileout_M);
    }

    /*
    string fileout = "4c.txt";
    ofile.open(fileout);
    set_initial_text();
    double T_final = 2.4;
    int MC_max = 10000;

    for (double Temperature = T_init; Temperature <= T_final; Temperature += 1.4){
        Run_simulation(L, Temperature, MC_max);
    }
    ofile.close();
    */

    /*
    for (double Temperature = T_init; Temperature <= T_final; Temperature += 1.4){
        cout << "Running for temperature: " << Temperature << endl;
            for (MC_cycles = 100; MC_cycles <= 100000; MC_cycles *= 10){
            start = clock();
            cout << "Using MC_cycles = " << MC_cycles << endl;
            Expectation_values = new double[5];
            Metropolis_method(L, MC_cycles, Temperature, Expectation_values, 1);
            write_file(L, MC_cycles, Temperature, Expectation_values);
            finish = clock();
            cout << "Time elapsed for MC_cycles = " << MC_cycles << ":  " <<
                    ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
        }
    }
    ofile.close();
    */
    return 0;
}
