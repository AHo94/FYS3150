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
    double DaystoYrs = 1.0/365.0;   // Converts from years to days
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
    vel /= DaystoYrs;   // Converts velocity from AU/days to AU/yrs
}

int main(){
    SolarSystem System;     // Initializes the solar system
    // Masses of the celestial bodies in the solar system
    double M_sun, M_earth, M_jupiter, M_mercury, M_venus, M_mars, M_saturn, M_uranus, M_neptune, M_pluto;
    M_sun = 1.0;
    double M_sun_real = 2*pow(10,30);
    M_earth = 6*pow(10,24)/M_sun_real;
    M_jupiter = 1.9*pow(10,27)/M_sun_real;
    M_mercury = 2.4*pow(10,23)/M_sun_real;
    M_venus = 4.9*pow(10,24)/M_sun_real;
    M_mars = 6.6*pow(10,23)/M_sun_real;
    M_saturn = 5.5*pow(10,26)/M_sun_real;
    M_uranus = 8.8*pow(10,25)/M_sun_real;
    M_neptune = 1.03*pow(10,26)/M_sun_real;
    M_pluto = 1.31*pow(10,22)/M_sun_real;

    // Creates arrays for names and masses for the celestials


    string Celestial_names[] = {"earth", "jupiter", "mercury", "venus", "mars", "saturn",
                            "uranus", "neptune", "pluto"};
    double Celestial_masses[] = {M_earth, M_jupiter, M_mercury, M_venus, M_mars, M_saturn,
                             M_uranus, M_neptune, M_pluto};

    // Hard coding the Sun into the system in the origin and at rest
    System.createCelestialBody(vec3(0,0,0), vec3(0,0,0), M_sun);

    // Adds all the celestials to the system
    for (int i=0; i<sizeof(Celestial_masses)/sizeof(*Celestial_masses); i++){
        vec3 position, velocity;
        set_initial_cond(position, velocity, Celestial_names[i]);
        System.createCelestialBody(position, velocity, Celestial_masses[i]);
    }

    // Solving system
    double dt = 0.0001;
    int NumTimesteps = 25000;
    ODEsolvers solver(dt);
    int plot_counter = 25;

    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 25){
            // Saves every 25 steps.
            System.write_file("TESTTT.txt");
            plot_counter = 0;
        }
        //solver.Euler_step(System);
        //solver.EulerCromer(System);
        solver.Verlet(System);
        plot_counter += 1;
    }

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

    /*
    SolarSystem MercurySys;
    MercurySys.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1);
    vec3 Mercpos(0.3075, 0, 0);
    vec3 Mercvel(0, 12.44, 0);
    MercurySys.createCelestialBody(Mercpos, Mercvel, M_mercury);
    dt = 0.0001;
    NumTimesteps = 100000;
    ODEsolvers GRsolver(dt);
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            MercurySys.write_file("Mercury_GR.txt");
            plot_counter = 0;
        }
        GRsolver.Verlet_GR(MercurySys);
        plot_counter += 1;
    }
    */

    return 0;
}
