#include <iostream>
#include <cmath>    // Math library
#include <fstream>  // Used to read files
#include <iomanip>  // setw identitation for output file
#include <string>
#include "eulermethod.h"
#include "vec3.h"
#include "celestials.h"

using namespace std;

int main()
{
    // Masses of the celestial bodies in the solar system
    eulermethod Euler_solver;
    double M_sun, M_earth, M_jupiter, M_mars;
    M_sun = 2*pow(10,30);
    M_earth = 6*pow(10,24);
    M_jupiter = 1.9*pow(10,27);
    M_mars = 6.6*pow(10,23);

    double AU = 1.5*pow(10,11); // Defining one astronomical unit in meters
    double R_earth = 1*AU;

    Celestials Earth;
    Earth.set_properties(M_earth, R_earth);
    cout << "Mass of Earth = " << Earth.get_mass() << endl;;

    vec3 Earth_pos;
    string earth_file = "../NASA_data/horizons_results.txt";
    std::ifstream infile(earth_file);
    std::string line;
    char c;
    int counter = 0;
    while (std::getline(infile, line)){
        if (line == "$$SOE"){
            counter += 1;
        }
        if (line == "$$EOE"){
            counter -= 1;
        }
        if (counter >= 1 && line != "$$SOE"){
            cout << line << endl;
        }
    }
    double **Epos, **Evel;
    int n = 1000;
    Epos = new double*[2];
    Evel = new double*[2];
    for (int i=0; i<2; i++){
        Epos[i] = new double[n];
        Evel[i] = new double[n];
    }
    Epos[0][0] = 1.479572465138224E+08;
    Epos[1][0] = 2.645009868848537E+07;
    Evel[0][0] = -5.659086230512699E+00;
    Evel[1][0] = 2.924889478278030E+01;
    Euler_solver.solve_euler_2D(Epos, Evel, n);
    ofstream datafile;
    datafile.open("TESTFILE.txt");
    for (int i=0; i<n; i++){
        datafile << Epos[0][i] << setw(15)
                 << Epos[1][i] << setw(15)
                 << Evel[0][i] << setw(15)
                 << Evel[1][i] << "\n";
    }
    datafile.close();

    return 0;
}
