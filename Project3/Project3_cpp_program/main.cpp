#include <iostream>
#include <cmath>    // Math library
#include "vec3.h"
#include "celestials.h"
#include "solarsystem.h"
#include "odesolvers.h"
using namespace std;

int main()
{
    // Masses of the celestial bodies in the solar system

    SolarSystem System;
    double M_sun, M_earth, M_jupiter, YrstoDays;
    YrstoDays = 1.0/365;  // Converts one year to 365 days
    M_sun = 1.0;
    M_earth = 6*pow(10,24)/(2*pow(10,30));
    M_jupiter = 1.9*pow(10,27)/(2*pow(10,30));

    System.createCelestialBody(vec3(0,0,0), vec3(0,0,0), M_sun);
    vec3 Earthpos (9.890331046925951E-01, 1.768079890757788E-01, -1.738715302893284E-04);
    vec3 Earthvel (-3.268395786841218E-03, 1.689265025904021E-02, -9.889230545039174E-07);
    Earthvel /= YrstoDays;
    System.createCelestialBody(Earthpos, Earthvel, M_earth);

    double dt = 0.01;
    int NumTimesteps = 1000;
    ODEsolvers solver(dt);
    for (int step=0; step<NumTimesteps; step++){
        System.write_file("TESTINGFILE.txt");
        //solver.Euler_step(System);
        //solver.EulerCromer(System);
        solver.Verlet(System);
    }


    //Celestials Earth(1,0,0,0,0,1);
    /*
    string earth_file = "../NASA_data/horizons_results.txt";

    std::ifstream infile(earth_file);
    std::string line;
    std::string delimeter = "  ";
    char c;
    int counter = 0;
    while (std::getline(infile, line)){
        std::stringstream ss(line);

        int n;
        std::vector<int> v;
        while (iss > n);
        {
            v.push_back(n);
        }

        if (line == "$$SOE"){
            counter += 1;
        }
        if (line == "$$EOE"){
            counter -= 1;
        }
        if (counter >= 1 && line != "$$SOE"){
            cout << std:: line << endl;
        }

    }
    */

    return 0;
}
