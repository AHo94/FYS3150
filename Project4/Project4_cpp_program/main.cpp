#include <iostream>
#include <chrono>
#include <random>
#include <cmath>

using namespace std;

void initialize_system(int N, int lattice, double **spins, double *energies){
    std::random_device rd; // obtain a random number from hardware
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<> distr(-1,1); // define the range
    for (int i=0; i<N; i++){
        for (int j=0; j<lattice; j++){
            if (distr(generator) == 0){
                spins[i][j] = -1;
            }
            else{
                spins[i][j] = 1;
            }
        }
    }

    for (int i=0; i<N; i++){
        energies[i] = 2*(spins[i][0]*spins[i][1] + spins[i][0]*spins[i][2]
                + spins[i][1]*spins[i][3] + spins[i][2]*spins[i][3]);
    }

}

int main()
{
    double **spins, *energies;
    double T = 1.0;     // Temperature = 1.0 kT/J
    double beta = 1/T;
    int lattice = 4;
    int N = pow(2, lattice);
    energies = new double[N];
    spins = new double*[N];
    for (int i=0; i<N; i++){
        spins[i] = new double[lattice];
    }
    initialize_system(N, lattice, spins, energies);
    double Z;
    for (int i=0; i<N; i++){
        Z += exp(-beta*energies[i]);
    }

    /*
    for (int i=0; i<N; i++){
        for (int j=0; j<lattice; j++){
            cout << spins[i][j];
        }
        cout << "\n";
    }
    */
    /*
    for (int i=0; i<N; i++){
        cout << energies[i] << " ";
    }
    */
    return 0;
}
