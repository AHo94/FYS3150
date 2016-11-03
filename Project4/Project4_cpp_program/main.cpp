#include <iostream>
#include <chrono>   // Used to seed random generator based on the time
#include <random>
#include <cmath>
#include <fstream>
#include <string>
//#include <armadillo>

using namespace std;

inline int periodic(int i, int limit, int add){
    // Funcion that ensures periodic boundary conditions
    return (i+limit+add) % (limit);
}

void initialize_system(int L, double **spins, double &energy, double &magnetic_moment){
    // Initializes the system, by setting spin states and calculates the energies and magnetizations.
    std::random_device rd; // obtain a random number from hardware
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> distr(0,1); // define the range
    for (int x=0; x < L; x++){
        for (int y=0; y<L; y++){
            int spin = distr(generator);
            if (spin == 0){
                spins[x][y] = -1;
            }
            else{
                spins[x][y] = 1;
            }
            magnetic_moment +=  spins[x][y];
        }
    }
    for (int x=0; x < L; x++){
        for (int y=0; y<L; y++){
            energy -= spins[x][y]*(spins[periodic(x,L,-1)][y] + spins[x][periodic(y,L,-1)]);
        }
    }
}

void Metropolis_method(int L, int MC_cycles, double Temperature, double *Expectation_values){
    // A function that uses the Metropolis method
    //std::random_device rd;
    //std::mt19937_64 generator(rd());
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Time dependent seed
    std::default_random_engine generator(seed);
    //std::uniform_int_distribution<int> distr(0,3); // define the range
    //std::normal_distribution<double> distr(0.0,1.0);
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    double **Spin_matrix;
    Spin_matrix = new double*[L];
    for (int i=0; i<L; i++){
        Spin_matrix[i] = new double[L];
    }
    double Energy=0;
    double Mag_moment=0;
    initialize_system(L, Spin_matrix, Energy, Mag_moment);

    for (int cycle=1; cycle<=MC_cycles; cycle++){
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                int ix = distr(generator)*L;
                int iy = distr(generator)*L;
                double E_b =2*Spin_matrix[ix][iy]*(
                            Spin_matrix[ix][periodic(iy, L, -1)]+
                            Spin_matrix[periodic(ix, L, -1)][iy] +
                            Spin_matrix[ix][periodic(iy, L, 1)]+
                            Spin_matrix[periodic(ix, L, 1)][iy]);

                Spin_matrix[ix][iy] *= -1.0;

                double E_t = 2*Spin_matrix[ix][iy]*(
                            Spin_matrix[ix][periodic(iy, L, -1)]+
                            Spin_matrix[periodic(ix, L, -1)][iy] +
                            Spin_matrix[ix][periodic(iy, L, 1)]+
                            Spin_matrix[periodic(ix, L, 1)][iy]);

                double dE = E_t - E_b;
                if (dE <= 0){
                    Energy += E_t;
                    Mag_moment += (Spin_matrix[0][0] + Spin_matrix[0][1] +
                            Spin_matrix[1][0] + Spin_matrix[1][1]);
                }
                else if(dE > 0){
                        if (distr(generator) <= exp(-dE/Temperature)){
                        Energy += E_t;
                        Mag_moment += (Spin_matrix[0][0] + Spin_matrix[0][1] +
                                Spin_matrix[1][0] + Spin_matrix[1][1]);
                        }

                    else{
                        Spin_matrix[ix][iy] *= -1.0;
                        Energy += E_b;
                        Mag_moment += (Spin_matrix[0][0] + Spin_matrix[0][1] +
                                Spin_matrix[1][0] + Spin_matrix[1][1]);
                    }
                }
            }
        }
    }
    /*
    for (int cycle=1; cycle<MC_cycles; cycle++){
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                int ix = (distr(generator)*L);
                int iy = (distr(generator)*L);
                int dE = 2*Spin_matrix[ix][iy]*(
                        Spin_matrix[ix][periodic(iy, L, -1)]+
                        Spin_matrix[periodic(ix, L, -1)][iy] +
                        Spin_matrix[ix][periodic(iy, L, 1)]+
                        Spin_matrix[periodic(ix, L, 1)][iy]);
                if (distr(generator) <= exp(-dE/Temperature)){
                    Spin_matrix[ix][iy] *= -1.0;
                    Mag_moment += 2*Spin_matrix[ix][iy];
                    Energy += dE;
                }
            }
        }
    }
    */
    Expectation_values[0] += Energy;
    Expectation_values[1] += Energy*Energy;
    Expectation_values[2] += Mag_moment;
    Expectation_values[3] += Mag_moment*Mag_moment;
    Expectation_values[4] += fabs(Mag_moment);
}

void write_file(int L, int MC_cycles, double Temperature, double *Expectation_values){
    double norm = 1.0/MC_cycles;
    double E_expectation = Expectation_values[0]*norm;
    double E2_expectation = Expectation_values[1]*norm;
    double M_expectation = Expectation_values[2]*norm;
    double M2_expectation = Expectation_values[3]*norm;
    double Mabs_expectation = Expectation_values[4]*norm;

    double E_variance = (E2_expectation - E_expectation*E_expectation)/L/L;
    double M_variance = (M2_expectation - M_expectation*M_expectation)/L/L;
}

int main()
{
    double **Spin_matrix, *Expectation_values;
    double T = 1.0;     // Temperature = 1.0 kT/J
    int L = 2;  // Number of spins
    Spin_matrix = new double*[L];
    Expectation_values = new double[5];
    for(int i=0; i<L; i++) {
        Spin_matrix[i] = new double[L];
    }

    int MC_cycles = 1000000;
    Metropolis_method(L, MC_cycles, T, Expectation_values);
    //initialize_system(N, lattice, spins, energies, magnetization);
    // Analytical expressions
    double AC_v = 64.0*(4+cosh(8))/(T*pow((cosh(8)+3), 2));
    double Achi = 32.0*(exp(8) + 1)/(T*(cosh(8) + 3));
    double norm = 1.0/(MC_cycles);
    double C_v = (Expectation_values[1] - Expectation_values[0]*Expectation_values[0])/L/L/T/T;
    double Chi = (Expectation_values[4] - Expectation_values[2]*Expectation_values[2])/L/L/T/T;

    cout << "Analytic C_v = " << AC_v << " Numerical C_v = " << C_v << endl;

    cout << "Analytic Chi = " << Achi << " Numerical Chi = " << Chi << endl;
    return 0;
}
