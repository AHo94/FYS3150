#include <iostream>
#include <cmath>    // Math library
#include <string>
#include <fstream>
#include <sstream>
#include "vec3.h"
#include "celestials.h"
#include "solarsystem.h"
#include "odesolvers.h"
#include <time.h>
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

void solve_system(SolarSystem &Solar_system, double dt, double N, string filename, string method){
    /*
    ODEsolvers solver(dt);
    string SNsteps = to_string(N) + " ";
    string Sdt = to_string(dt) + "\n ";
    Solar_system.write_file(filename, SNsteps, Sdt);
    int plot_counter = 100;
    if (method == "verlet"){
        cout << "Running Verlet method" << endl;
        for (int step=0; step<N; step++){
            if (plot_counter == 100){
                // Saves every 100 steps.
                Solar_system.write_file(filename, SNsteps, Sdt);
                plot_counter = 0;
                SNsteps = "";
                Sdt = "";
            }
            solver.Verlet(Solar_system);
            plot_counter += 1;
        }
    }

    else if(method == "euler"){
        cout << "Running Euler method" << endl;
        for (int step=0; step<N; step++){
            if (plot_counter == 100){
                // Saves every 100 steps.
                Solar_system.write_file(filename, SNsteps, Sdt);
                plot_counter = 0;
                SNsteps = "";
                Sdt = "";;
            }
            solver.Euler_step(Solar_system);
            plot_counter += 1;
        }
    }
    else if (method == "eulercromer"){
        cout << "Running Euler Cromer method" << endl;
        for (int step=0; step<N; step++){
            if (plot_counter == 100){
                // Saves every 100 steps.
                Solar_system.write_file(filename, SNsteps, Sdt);
                plot_counter = 0;
                SNsteps = "";
                Sdt = "";
            }
            solver.EulerCromer(Solar_system);
            plot_counter += 1;
        }
    }
    else if (method == "verletGR"){
        cout << "Running Verlet method with GR correction" << endl;
        for (int step=0; step<N; step++){
            if (plot_counter == 100){
                // Saves every 100 steps.
                Solar_system.write_file(filename, SNsteps, Sdt);
                plot_counter = 0;
                SNsteps = "";
                Sdt = "";
            }
            solver.Verlet_GR(Solar_system);
            plot_counter += 1;
        }
    }
    else{
        cout << "Method input invalid, input = " << method << ", try the following: \n"
             << "verlet \n"
             << "euler \n"
             << "eulercromer \n"
             << "verletGR" <<endl;
        terminate();
    }
    */
}

int main(){
    clock_t start, finish;
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

    // Hard coding position and velocity of the Sun. Assuming at origin and at rest
    vec3 Sunpos(0,0,0);
    vec3 Sunvel(0,0,0);

    // Creates arrays for names and masses for the celestials
    string Celestial_names[] = {"earth", "jupiter", "mercury", "venus", "mars", "saturn",
                            "uranus", "neptune", "pluto"};
    double Celestial_masses[] = {M_earth, M_jupiter , M_mercury, M_venus, M_mars, M_saturn,
                             M_uranus, M_neptune, M_pluto};



    SolarSystem System;     // Initializes the solar system
    // Hard coding the Sun into the system in the origin and at rest
    System.createCelestialBody(Sunpos, Sunvel, M_sun);

    // Adds all the celestials to the system
    for (int i=0; i<sizeof(Celestial_masses)/sizeof(*Celestial_masses); i++){
        vec3 position, velocity;
        set_initial_cond(position, velocity, Celestial_names[i]);
        System.createCelestialBody(position, velocity, Celestial_masses[i]);
    }

    // Solving system
    start = clock();
    double dt = 0.001;
    int NumTimesteps = 100000;
    ODEsolvers solver(dt);
    string SNsteps = to_string(NumTimesteps) + " ";
    string Sdt = to_string(dt) + "\n ";
    //solve_system(System, dt, NumTimesteps, "Celestial_positions.txt", "verlet");
    finish = clock();
    cout << "Time elapsed for non interactive case: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

    // Solving Earth-Jupiter and Sun system
    SolarSystem EJ_SunSys;
    EJ_SunSys.createCelestialBody(Sunpos, Sunvel, M_sun);
    string EJ_names[] = {"earth", "jupiter"};
    double EJ_masses[] = {M_earth, M_jupiter};

    for (int i=0; i<sizeof(EJ_masses)/sizeof(*EJ_masses); i++){
        vec3 position, velocity;
        set_initial_cond(position, velocity, EJ_names[i]);
        EJ_SunSys.createCelestialBody(position, velocity, EJ_masses[i]);
    }
    string filename = "Earth_Sun_Jupiter.txt";
    int plot_counter = 100;
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            // Saves every 100 steps.
            EJ_SunSys.write_file(filename);//, SNsteps, Sdt);
            plot_counter = 0;
            SNsteps = "";
            Sdt = "";
        }
        solver.Verlet(EJ_SunSys);
        plot_counter += 1;
    }
    //solve_system(EJ_SunSys, dt, NumTimesteps, "Earth_Sun_Jupiter.txt", "verlet");
    /*
    ODEsolvers solver(dt);
    int plot_counter = 100;
    string init_text = to_string(NumTimesteps) + " " +  to_string(dt) + " POOP \n";
    string SNsteps = to_string(NumTimesteps) + " ";
    string Sdt = to_string(dt) + "\n";
    // Solving for Sun - Earth system - Euler method
    System.write_file("Testing.txt", SNsteps, Sdt);
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            // Saves every 25 steps.
            //cout << "write" << endl;
            System.write_file("Testing.txt", SNsteps, Sdt);
            plot_counter = 0;
            SNsteps = "";
            Sdt = "";
        }
        solver.Euler_step(System);
        plot_counter += 1;
    }


    // Solving for Sun - Earth system - Euler Cromer
    System.write_file("Earth_Sun_EulerCromer.txt", init_text);
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            // Saves every 25 steps.
            System.write_file("Earth_Sun_EulerCromer.txt", "");
            plot_counter = 0;
        }
        solver.EulerCromer(System);
        plot_counter += 1;
    }

    // Solving for Sun - Earth system - Verlet method
    System.write_file("Earth_Sun_Verlet.txt", init_text);
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            // Saves every 25 steps.
            System.write_file("Earth_Sun_Verlet.txt", "");
            plot_counter = 0;
        }
        solver.Verlet(System);
        plot_counter += 1;
    }

    // Solving for all planets
    //System.write_file("Celestial_positions.txt", init_text);

    int plot_counter = 100;
    ODEsolvers solver(dt);
    System.write_file("Celestial_positions.txt", "HEYAAAAAAA \n");
    for (int step=0; step<NumTimesteps; step++){
        if (plot_counter == 100){
            // Saves every 100 steps.
            System.write_file("Celestial_positions.txt", "");
            plot_counter = 0;
        }
        solver.Verlet(System);
        plot_counter += 1;
    }
    */
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
