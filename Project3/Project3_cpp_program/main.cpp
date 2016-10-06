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
    M_jupiter = 0;//1.9*pow(10,27)/(2*pow(10,30));

    System.createCelestialBody(vec3(0,0,0), vec3(0,0,0), M_sun);
    vec3 Earthpos (9.890331046925951E-01, 1.768079890757788E-01, -1.738715302893284E-04);
    vec3 Earthvel (-3.268395786841218E-03, 1.689265025904021E-02, -9.889230545039174E-07);
    Earthvel /= YrstoDays;
    System.createCelestialBody(Earthpos, Earthvel, M_earth);

    vec3 Jupiterpos (-5.433468170028908E+00, -3.819061221110369E-01, 1.231004384238452E-01);
    vec3 Jupitervel (4.425651679847022E-04, -7.171108917491057E-03, 1.992744446163222E-05);
    Jupitervel /= YrstoDays;
    System.createCelestialBody(Jupiterpos, Jupitervel, M_jupiter);

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
