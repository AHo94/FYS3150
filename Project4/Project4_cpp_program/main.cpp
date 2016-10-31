#include <iostream>
#include <chrono>
#include <random>
#include <cmath>

using namespace std;

void initialize_system(int N, int lattice, double **spins, double *energies, double *magnetization){
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
        energies[i] = -2*(spins[i][0]*spins[i][1] + spins[i][0]*spins[i][2]
                + spins[i][1]*spins[i][3] + spins[i][2]*spins[i][3]);
        magnetization[i] = (spins[i][0] + spins[i][1] + spins[i][2] + spins[i][3]);
    }

}

int main()
{
    double **spins, *energies, *magnetization;
    double T = 1.0;     // Temperature = 1.0 kT/J
    double beta = 1/T;
    int lattice = 4;
    int N = pow(2, lattice);
    energies = new double[N];
    magnetization = new double[N];
    spins = new double*[N];
    for (int i=0; i<N; i++){
        spins[i] = new double[lattice];
    }

    initialize_system(N, lattice, spins, energies, magnetization);
    double Z, E_mean, E_mean2, M_mean, M_mean2, C_v, chi, totC_v, totChi;
    double totE, totM;

    // Analytical expressions
    double AC_v = 256.0/(T*pow(cosh(8)+3, 2));
    double Achi = 32.0*(exp(8) + 1)/(T*(cosh(8) + 3));

    for (int j=0; j<1000; j++){
    for (int i=0; i<N; i++){
        Z += exp(-beta*energies[i]);
    }
    for (int i=0; i<N; i++){
        E_mean += (1.0/Z)*(energies[i]*exp(-beta*energies[i]));
        E_mean2 += (1.0/Z)*(energies[i]*energies[i]*exp(-beta*energies[i]));
        M_mean += (1.0/Z)*(magnetization[i]*exp(-beta*energies[i]));
        M_mean2 += (1.0/Z)*(magnetization[i]*magnetization[i]*exp(-beta*energies[i]));
    }

    C_v = (E_mean2 - pow(E_mean, 2))/T;
    chi = (M_mean2 - pow(M_mean, 2))/T;

    totC_v += C_v;
    totChi += chi;
    if (fabs(totC_v/j - AC_v) < 1e-3 && fabs(totChi/j - Achi) < 1e-3){
        j = 999;
        cout << "AYY WORKED" << endl;
    }
    }
    cout << "Numerical C_v = " << totC_v/100.0 << ", Analytical C_v = " << AC_v << endl;
    cout << "Numerical chi = " << totChi/100.0 << ", Analytical chi = " << Achi << endl;

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
