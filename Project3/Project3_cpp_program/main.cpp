#include <iostream>
#include <cmath>    // Math library
#include <string>
#include <fstream>
#include <sstream>
#include "vec3.h"
#include "celestials.h"
#include "solarsystem.h"
#include "odesolvers.h"
using namespace std;

void set_initial_cond(vec3 &pos, vec3 &vel, string planet_name){
    /*
     * Function that reads data from NASA and sets the coordinate/velocities as initial conditions.
     * Datafile is assumed to only contain positions in the first row and velocities in the second row.
     * Source of data: http://ssd.jpl.nasa.gov/horizons.cgi
     */
    string data_filename = "../NASA_data/";
    data_filename.append(planet_name);
    data_filename.append("_data.txt");
    ifstream infile (data_filename);
    string line;
    double YrstoDays = 1.0/365.0;   // Converts from years to days
    double x, y, z;
    int counter = 0;
    while (getline(infile, line)){
        istringstream ss(line);
        ss >> x >> y >> z;
        if (counter == 0)
        {
            pos(0) = x;
            pos(1) = y;
            pos(2) = z;
        }
        else if (counter == 1){
            vel(0) = x;
            vel(1) = y;
            vel(2) = z;
        }
        counter += 1;
    }
    vel /= YrstoDays;   // Converts velocity from AU/days to AU/yrs
}

int main(){

    SolarSystem System;
    // Masses of the celestial bodies in the solar system
    double M_sun, M_earth, M_jupiter, M_mercury, M_venus, M_mars, M_saturn, M_uranus, M_neptune, M_pluto;
    M_sun = 1.0;
    M_earth = 6*pow(10,24)/(2*pow(10,30));
    M_jupiter = 1.9*pow(10,27)/(2*pow(10,30));
    M_mercury = 2.4*pow(10,23)/(2*pow(10,30));
    M_venus = 4.9*pow(10,24)/(2*pow(10,30));
    M_mars = 6.6*pow(10,23)/(2*pow(10,30));
    M_saturn = 5.5*pow(10,26)/(2*pow(10,30));
    M_uranus = 8.8*pow(10,25)/(2*pow(10,30));
    M_neptune = 1.03*pow(10,26)/(2*pow(10,30));
    M_pluto = 1.31*pow(10,22)/(2*pow(10,30));


    // Adding the Sun
    System.createCelestialBody(vec3(0,0,0), vec3(0,0,0), M_sun);

    // Adding Earth
    vec3 Earthpos, Earthvel;
    set_initial_cond(Earthpos, Earthvel, "earth");
    System.createCelestialBody(Earthpos, Earthvel, M_earth);

    // Adding Jupiter
    vec3 Jupiterpos, Jupitervel;
    set_initial_cond(Jupiterpos, Jupitervel, "jupiter");
    System.createCelestialBody(Jupiterpos, Jupitervel, M_jupiter);

    // Adding Mercury
    vec3 Mercurypos, Mercuryvel;
    set_initial_cond(Mercurypos, Mercuryvel, "mercury");
    System.createCelestialBody(Mercurypos, Mercuryvel, M_mercury);

    // Adding Venus
    vec3 Venuspos, Venusvel;
    set_initial_cond(Venuspos, Venusvel, "venus");
    System.createCelestialBody(Venuspos, Venusvel, M_venus);

    // Adding Mars
    vec3 Marspos, Marsvel;
    set_initial_cond(Marspos, Marsvel, "mars");
    System.createCelestialBody(Marspos, Marsvel, M_mars);

    // Adding Saturn
    vec3 Saturnpos, Saturnvel;
    set_initial_cond(Saturnpos, Saturnvel, "saturn");
    System.createCelestialBody(Saturnpos, Saturnvel, M_saturn);

    // Adding Uranus
    vec3 Uranuspos, Uranusvel;
    set_initial_cond(Uranuspos, Uranusvel, "uranus");
    System.createCelestialBody(Uranuspos, Uranusvel, M_uranus);

    // Adding Neptune
    vec3 Neptunepos, Neptunevel;
    set_initial_cond(Neptunepos, Neptunevel, "neptune");
    System.createCelestialBody(Neptunepos, Neptunevel, M_neptune);

    // Adding pluto
    vec3 Plutopos, Plutovel;
    set_initial_cond(Plutopos, Plutovel, "pluto");
    System.createCelestialBody(Plutopos, Plutovel, M_pluto);

    // Solving system
    double dt = 0.001;
    int NumTimesteps = 10000;
    ODEsolvers solver(dt);
    int plot_counter = 10;
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 25){
            // Saves every 25 steps.
            System.write_file("Celestial_positions.txt");
            plot_counter = 0;
        }
        //solver.Euler_step(System);
        //solver.EulerCromer(System);
        solver.Verlet(System);
        plot_counter += 1;
    }


    return 0;
}
