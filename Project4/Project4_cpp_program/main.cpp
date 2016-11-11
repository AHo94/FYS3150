#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <mpi.h>

using namespace std;
ofstream ofile_global;

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
    double CurrentE = 0;
    for (int i=0; i<L; i++){
         for (int j=0; j<L; j++){
             CurrentE -= Spin_matrix[i][j]*(Spin_matrix[periodic(i, L, -1)][j] +
                                       Spin_matrix[i][periodic(j, L, -1)]);
        }
    }
    return CurrentE;
}

double Calculate_M(double **Spin_matrix, int L){
    double CurrentM = 0;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            CurrentM += Spin_matrix[i][j];
        }
    }
    return CurrentM;
}

void Metropolis_method(int L, int MC_cycles, double Temperature, double *Expectation_values, double *accepted_config,
                       double *Energies_array, double *Mag_moments_array, int random_state = 0){
    // A function that uses the Metropolis method
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    double **Spin_matrix;
    Spin_matrix = new double*[L];
    for (int i=0; i<L; i++){
        Spin_matrix[i] = new double[L];
    }
    int accepted_flip = 0;
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
                int dE =  2*Spin_matrix[ix][iy]*
                      (Spin_matrix[ix][periodic(iy,L,-1)]+
                       Spin_matrix[periodic(ix,L,-1)][iy] +
                       Spin_matrix[ix][periodic(iy,L,1)] +
                       Spin_matrix[periodic(ix,L,1)][iy]);
                if (distr(generator) <= exp(-dE/Temperature)){
                    Spin_matrix[ix][iy] *= -1.0;
                    currentM += 2*Spin_matrix[ix][iy];//Calculate_M(Spin_matrix, L);
                    currentEnergy += dE;
                    accepted_flip += 1;
                }
            }
        }
        // Sums up the samplings
        EnergySum += currentEnergy;
        EnergySquaredSum += currentEnergy*currentEnergy;
        MSum += currentM;
        fabsMSum += fabs(currentM);
        MSquaredSum += currentM*currentM;

        Energies_array[cycle] = currentEnergy;
        Mag_moments_array[cycle] = fabs(currentM);
        accepted_config[cycle] = accepted_flip;
    }
    Expectation_values[0] = EnergySum;
    Expectation_values[1] = EnergySquaredSum;
    Expectation_values[2] = MSum;
    Expectation_values[3] = MSquaredSum;
    Expectation_values[4] = fabsMSum;
}

void Metropolis_parallelization(int L, double Temperature, double *Expectation_values, int random_state = 0){
    // Function that solves the Metropolis method specifically for the parallellization part.
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    double **Spin_matrix;
    Spin_matrix = new double*[L];
    for (int i=0; i<L; i++){
        Spin_matrix[i] = new double[L];
    }
    initialize_system(L, Spin_matrix, random_state);

    double currentEnergy = Calculate_E(Spin_matrix, L);
    double currentM = Calculate_M(Spin_matrix, L);
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            int ix = distr(generator)*L;
            int iy = distr(generator)*L;
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
        }
    }
    // Adds samplings
    Expectation_values[0] += currentEnergy;
    Expectation_values[1] += currentEnergy*currentEnergy;
    Expectation_values[2] += currentM;
    Expectation_values[3] += currentM*currentM;
    Expectation_values[4] += fabs(currentM);
}

void write_file(int L, double T, int MC_cycles, double *accepted_flip, double *Mean_energies, double *Mean_mag_moments,
                double *Expectation_values, string filename_E, string filename_M){
    cout << "Saving data to file" << endl;
    double norm = 1.0/MC_cycles;
    double Variance_E = Expectation_values[1]*norm - Expectation_values[0]*norm*Expectation_values[0]*norm;
    double Variance_M = Expectation_values[3]*norm - Expectation_values[4]*norm*Expectation_values[4]*norm;

    ofstream ofileE, ofileM; // Output file
    ofileE.open(filename_E);
    ofileE << "<E>" << setw(15) << "MC Cycles" << setw(15);
    ofileE << "# Spins" << setw(15) << "Temperature" << setw(15) << "Variance_E" << "\n";
    ofileE << setw(15) << MC_cycles << setw(15) << L << setw(15) << T << setw(15) << Variance_E << "\n";

    ofileM.open(filename_M);
    ofileM << "<|M|>" << setw(15) << "MC Cycles" << setw(15);
    ofileM << "# Spins" << setw(15) << "Temperature" << setw(15) << "Variance_M" << "\n";
    ofileM << setw(15) << MC_cycles << setw(15) << L << setw(15) << T << setw(15) << Variance_M <<"\n";
    int counter = 100;
    for (int i=0; i<MC_cycles; i++){
        if (counter == 100){
            // Saves every 100 point
            ofileE << Mean_energies[i]  << setw(15) << accepted_flip[i] << "\n";
            ofileM << Mean_mag_moments[i] << setw(15) << accepted_flip[i] << "\n";
            counter = 0;
        }
        counter += 1;
    }

    ofileE.close();
    ofileM.close();
}

void write_file_4d(int L, double T, int MC_cycles, double *Expectation_values, string filename){
    float norm = 1.0/MC_cycles;
    double variance = (Expectation_values[1]*norm -
            Expectation_values[0]*Expectation_values[0]*norm*norm)/L/L;
    ofile_global << Expectation_values[0]*norm/L/L << setw(15) << variance;
}

void write_file_parallellization(int L, double T, int MC_cycles, double *Total_expectation_values){
    double norm = 1.0/MC_cycles;
    double E_expect = Total_expectation_values[0]*norm;
    double E_expect_2 = Total_expectation_values[1]*norm;
    double M_expect_2 = Total_expectation_values[2]*norm;
    double M_abs_expect = Total_expectation_values[4]*norm;

    double E_variance = E_expect_2 - E_expect*E_expect;
    double M_variance = M_expect_2 - M_abs_expect*M_abs_expect;
    double C_v = E_variance/T;
    double Chi = M_variance/T;
    ofile_global << setw(15) << T;
    ofile_global << setw(15) << L;
    ofile_global << setw(15) << E_expect;
    ofile_global << setw(15) << M_abs_expect;
    ofile_global << setw(15) << C_v;
    ofile_global << setw(15) << Chi << endl;;

}

int main(int nargs, char*args[])
{

    clock_t start, finish;
    double *Expectation_values;
    Expectation_values = new double[5];
    int L = 2;  // Number of spins
    int MC_cycles = 0;  // Number of Monte Carlo cycles

    double *Energies_array, *Mag_moments_array, *accepted_config;
    if (nargs <= 1){
        // Analytical expressions
        double T_init = 1.0; // Temperature = 1.0 kT/J
        double AC_v = 64.0*(1+3*cosh(8.0/T_init))/(T_init*pow((cosh(8.0/T_init)+3), 2));
        double Achi = 8*(exp(8.0/T_init) + cosh(8.0/T_init) + 3.0/2.0)/(T_init*pow((cosh(8.0/T_init)+3), 2));
        int MC_cycles = 1000000;
        Energies_array = new double[MC_cycles];
        Mag_moments_array = new double [MC_cycles];
        accepted_config = new double[MC_cycles];

        start = clock();
        Metropolis_method(L, MC_cycles, T_init, Expectation_values, accepted_config, Energies_array, Mag_moments_array);
        double C_v = (Expectation_values[1]/MC_cycles -
            Expectation_values[0]*Expectation_values[0]/MC_cycles/MC_cycles)/T_init/T_init;
        double Chi = (Expectation_values[3]/MC_cycles -
            Expectation_values[4]*Expectation_values[4]/MC_cycles/MC_cycles)/T_init/T_init;
        cout << "Number of Monte Carlo cycles = " << MC_cycles << endl;
        cout << "Analytic C_v = " << AC_v << ", Numerical C_v = " << C_v << endl;
        cout << "Analytic Chi = " << Achi << ", Numerical Chi = " << Chi << endl;

        finish = clock();
        cout << "Time elapsed for MC_cycles = " << MC_cycles << ":  " <<
            ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;


        cout << "\n" << "Running multiple times, using MC_cycles = "<< MC_cycles << endl;
        for (int i=0; i<=5; i++)
        {
            // Runs this example multiple times to showcase the stability of the algorithm
            Expectation_values = new double[5];
            Metropolis_method(L, MC_cycles, T_init, Expectation_values, accepted_config, Energies_array, Mag_moments_array);
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
        MC_cycles = 1000000;
        double T_final = 2.4;
        string filename_E = "E_expect_T";
        string filename_M = "M_expect_T";
        // Runs for an arbritary state.
        cout << "Doing calculations with abritary intial state" << endl;
        for (double Temperature = T_init; Temperature <= T_final; Temperature += 1.4){
            cout << "Running for T = " << Temperature << endl;
            start = clock();
            Expectation_values = new double[5];
            Energies_array = new double[MC_cycles];
            Mag_moments_array = new double [MC_cycles];
            Metropolis_method(L, MC_cycles, Temperature, Expectation_values,
                               accepted_config, Energies_array, Mag_moments_array);
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
            write_file(L, Temperature, MC_cycles, accepted_config, Energies_array, Mag_moments_array,
                   Expectation_values, fileout_E, fileout_M);
        }


        string filename_E_allup = "E_expect_AllUpState_T";
        string filename_M_allup = "M_expect_AllUpState_T";
        // Running for an initial state where all spins point up
        cout << "Doing calculations with initial state with all spins up" << endl;
        for (double Temperature = T_init; Temperature <= T_final; Temperature += 1.4){
            cout << "Running for T = " << Temperature << endl;
            start = clock();
            Expectation_values = new double[5];
            Energies_array = new double[MC_cycles];
            Mag_moments_array = new double [MC_cycles];
            Metropolis_method(L, MC_cycles, Temperature, Expectation_values,
                               accepted_config, Energies_array, Mag_moments_array, 1);
            finish = clock();
            cout << "Time elapsed for MC_cycles = " << MC_cycles << ":  " <<
                ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

            string fileout_E = filename_E_allup;
            string fileout_M = filename_M_allup;
            stringstream stream;
            stream << fixed << setprecision(2) << Temperature;
            string argument = stream.str();
            fileout_E.append(argument);
            fileout_E.append(".txt");
            fileout_M.append(argument);
            fileout_M.append(".txt");
            write_file(L, Temperature, MC_cycles, accepted_config, Energies_array, Mag_moments_array,
                   Expectation_values, fileout_E, fileout_M);
        }
    }
    else{
        // If there are arguments in the input line, then run for the last two tasks.
        // Reads L and MC cycles from command line and sets up temperature limits.
        L = atoi(args[1]);
        MC_cycles = atoi(args[2]);
        double T_init = 2.0;
        double T_final = 2.3;
        double Temp_step = 0.02;

        // Initialize filename to fit the number of spins L
        string filename = "../build-Project4_cpp_program-Desktop_Qt_5_7_0_GCC_64bit-Debug/4e_data_L";
        stringstream stream;
        stream << fixed << setprecision(0) << L;
        string argument = stream.str();
        filename.append(argument);
        filename.append(".txt");

        // New empty arrays
        double *Total_expectation_values;
        Expectation_values = new double[5];
        Total_expectation_values = new double[5];

        int numprocs, my_rank;
        //  Initialize MPI
        cout << "Wait for MPI" << endl;
        MPI_Init(&nargs, &args);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        if (my_rank == 0){
            ofile_global.open(filename);
            // IF TO APPEND:  std::ios_base::app
        }

        // Determine number of intervals for the nodes
        int no_intervals = MC_cycles/numprocs;
        int loop_begin = my_rank*no_intervals + 1;
        int loop_end = (my_rank+1)*no_intervals;
        if ((my_rank == numprocs - 1) && (loop_end < MC_cycles)){
            loop_end = MC_cycles;
        }

        // Broadcast common variables to all nodes
        MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&T_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&T_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (my_rank == 0){
            cout << "Everything initialized, startig algorithm with:" << endl;
            cout << "L = " << L << endl;
            cout << "MC_cycles = "<< MC_cycles << endl;
            cout << "This may take a while..." << endl;
        }
        double Time_start, Time_end, Time_total;
        Time_start = MPI_Wtime();
        for (double temperature = T_init; temperature <= T_final; temperature += Temp_step){
            if (fabs(temperature - 2.2) <= 1e-7){
                // Changes temperature step length when T = 2.2. Interesting things happens between T = 2.2 and T = 2.3.
                if (my_rank == 0){
                    cout << "Changing temperature step length" << endl;
                }
                Temp_step = 0.01;
            }

            for (int cycles = loop_begin; cycles <= loop_end; cycles ++)
                Metropolis_parallelization(L, temperature, Expectation_values);

            for (int i = 0; i<5; i++){
                // Merges all values from the different nodes
                MPI_Reduce(&Expectation_values[i], &Total_expectation_values[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            if(my_rank == 0){
                // Write to output file
                write_file_parallellization(L, temperature, MC_cycles, Total_expectation_values);
            }
        }
        ofile_global.close();
        Time_end = MPI_Wtime();
        Time_total = Time_end - Time_start;
        if (my_rank == 0){
            cout << "Time = " << Time_total << " on number of processors: " << numprocs << endl;
        }
        //  End MPI
        MPI_Finalize ();
    }
    return 0;
}
/*
SELF NOTE:
Compile with
mpic++ -std=c++11 ./main.x main.cpp
Run with
mpiexec -n 4 ./main.x (arguments)
Example
mpiexec -n 4 ./main.x 20 10000
*/
