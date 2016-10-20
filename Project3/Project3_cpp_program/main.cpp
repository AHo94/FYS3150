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

void solve_systems(SolarSystem &SolSys, int N, double dt, string filename, string method){
    clock_t start, finish;
    ODEsolvers solver(dt);
    int plot_counter = 0;
    string SNumsteps = to_string(N) + " ";
    string Sdt = to_string(dt) + "\n";
    SolSys.write_file(filename, SNumsteps, Sdt);
    if (method == "verlet"){
        cout << "Running Verlet method" << endl;
        start = clock();
        for (int step=0; step<N; step++){
            if (plot_counter == 20){
                // Saves every 100 steps.
                SolSys.write_file(filename, SNumsteps, Sdt);
                plot_counter = 0;
                SNumsteps = "";
                Sdt = "";
            }
            solver.Verlet(SolSys);
            plot_counter += 1;
        }
        finish = clock();
        cout << "Time elapsed for Verlet method: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    }
    else if (method == "euler"){
        cout << "Running Euler method" << endl;
        start = clock();
        for (int step=0; step<N; step++){
            if (plot_counter == 20){
                // Saves every 100 steps.
                SolSys.write_file(filename, SNumsteps, Sdt);
                plot_counter = 0;
                SNumsteps = "";
                Sdt = "";
            }
            solver.Euler_step(SolSys);
            plot_counter += 1;
        }
        finish = clock();
        cout << "Time elapsed for Euler method: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    }
    else if (method == "eulercromer"){
        cout << "Running Euler Cromer method" << endl;
        start = clock();
        for (int step=0; step<N; step++){
            if (plot_counter == 20){
                // Saves every 100 steps.
                SolSys.write_file(filename, SNumsteps, Sdt);
                plot_counter = 0;
                SNumsteps = "";
                Sdt = "";
            }
            solver.EulerCromer(SolSys);
            plot_counter += 1;
        }
        finish = clock();
        cout << "Time elapsed for Euler Cromer method: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    }
    else if (method == "verletGR"){
        cout << "Running Verlet method with GR adjustment" << endl;
        start = clock();
        for (int step=0; step<N; step++){
            if (plot_counter == 20){
                // Saves every 100 steps.
                SolSys.write_file(filename, SNumsteps, Sdt);
                plot_counter = 0;
                SNumsteps = "";
                Sdt = "";
            }
            solver.Verlet_GR(SolSys);
            plot_counter += 1;
        }
        finish = clock();
        cout << "Time elapsed for Verlet method w. GR adjustment: "
             << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    }
    else {
        cout << "Variable 'method' = " << method << " not properly set, try the following:" << endl;
        cout << "verlet \n" <<
                "euler \n" <<
                "eulercromer" <<
                "verletGR" << endl;
    }
    cout << "Data saved to: " << filename << "\n" << endl;
}

void New_system_and_solve(int N, double dt, string *names, double *masses, int NumCelestials
                          , string filename, string method){
    SolarSystem System;
    // Hardcoding position and velocity of the Sun. Assuming at origin and at rest.
    //vec3 SunPos(0,0,0);
    //vec3 SunVel(0,0,0);

    // Coding in the Sun
    //System.createCelestialBody(SunPos, SunVel, 1);

    for (int i=0; i<NumCelestials; i++){
        vec3 position, velocity;
        set_initial_cond(position, velocity, names[i]);
        System.createCelestialBody(position, velocity, masses[i]);
    }

    // Solve system
    solve_systems(System, N, dt, filename, method);
}

int main(){
    // Masses of the celestial bodies (excluding the Sun, as M_sun = 1) in the solar system
    double M_sun, M_earth, M_jupiter, M_mercury, M_venus, M_mars, M_saturn, M_uranus, M_neptune, M_pluto;
    double M_sun_real = 2*pow(10,30);
    M_sun = 1;
    M_earth = 6*pow(10,24)/M_sun_real;
    M_jupiter = 1.9*pow(10,27)/M_sun_real;
    M_mercury = 2.4*pow(10,23)/M_sun_real;
    M_venus = 4.9*pow(10,24)/M_sun_real;
    M_mars = 6.6*pow(10,23)/M_sun_real;
    M_saturn = 5.5*pow(10,26)/M_sun_real;
    M_uranus = 8.8*pow(10,25)/M_sun_real;
    M_neptune = 1.03*pow(10,26)/M_sun_real;
    M_pluto = 1.31*pow(10,22)/M_sun_real;

    // Earth-Sun system
    double dt = 0.0001;
    int NumTimesteps = 100000;
    string earth_sun_name[] = {"sun","earth"};
    double earth_sun_mass[] = {M_sun, M_earth};
    int NumCelestials = sizeof(earth_sun_mass)/sizeof(*earth_sun_mass);

    New_system_and_solve(NumTimesteps, dt, earth_sun_name, earth_sun_mass,
                         NumCelestials, "Earth_Sun_sys_euler.txt", "euler");
    New_system_and_solve(NumTimesteps, dt, earth_sun_name, earth_sun_mass,
                         NumCelestials, "Earth_Sun_sys_eulercromer.txt", "eulercromer");
    New_system_and_solve(NumTimesteps, dt, earth_sun_name, earth_sun_mass,
                         NumCelestials, "Earth_Sun_sys_verlet.txt", "verlet");

    // Increase dt as a comparison
    dt = 0.001;
    NumTimesteps = 10000;
    New_system_and_solve(NumTimesteps, dt, earth_sun_name, earth_sun_mass,
                         NumCelestials, "Earth_Sun_sys_euler_largerdt.txt", "euler");
    New_system_and_solve(NumTimesteps, dt, earth_sun_name, earth_sun_mass,
                         NumCelestials, "Earth_Sun_sys_verlet_largerdt.txt", "verlet");

    // Finding the escape velocity of a planet

    dt= 0.0001;
    NumTimesteps = 100000;
    vec3 PlanetPos(1,0,0);
    double EscapeVel = sqrt(8*acos(-1)*acos(-1)/(PlanetPos.length()));
    vec3 PlanetVel(0,EscapeVel,0);
    SolarSystem EscapeVelSystem;
    EscapeVelSystem.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1);
    EscapeVelSystem.createCelestialBody(PlanetPos, PlanetVel, M_earth);
    solve_systems(EscapeVelSystem, NumTimesteps, dt, "Escape_velocity_system.txt", "verlet");


    // Earth-Sun-Jupiter system
    dt = 0.0001;
    NumTimesteps = 150000;
    string ESJ_names[] = {"sun", "earth", "jupiter"};
    double ESJ_masses[] = {M_sun, M_earth, M_jupiter};
    NumCelestials = sizeof(ESJ_masses)/sizeof(*ESJ_masses);
    New_system_and_solve(NumTimesteps, dt, ESJ_names, ESJ_masses, NumCelestials, "ESJ_sys.txt", "verlet");

    // Stability test of verlet method in Earth-Sun-Jupiter system
    dt = 0.001;
    NumTimesteps = 15000;
    New_system_and_solve(NumTimesteps, dt, ESJ_names, ESJ_masses, NumCelestials, "ESJ_sys_largerdt.txt", "verlet");

    dt = 0.0001;
    NumTimesteps = 150000;
    // Jupiter mass now 10 times larger
    double ESJ_masses2[] = {M_sun, M_earth, 10*M_jupiter};
    New_system_and_solve(NumTimesteps, dt, ESJ_names, ESJ_masses2, NumCelestials, "ESJ_sys_10MJ.txt", "verlet");
    // Jupiter mass now 1000 times larger
    double ESJ_masses3[] = {M_sun, M_earth, 1000*M_jupiter};
    New_system_and_solve(NumTimesteps, dt, ESJ_names, ESJ_masses3, NumCelestials, "ESJ_sys_1000MJ.txt", "verlet");


    // Whole Solar system
    string Celestial_names[] = {"sun", "earth", "jupiter", "mercury", "venus", "mars", "saturn",
                            "uranus", "neptune", "pluto"};
    double Celestial_masses[] = {M_sun, M_earth, M_jupiter, M_mercury, M_venus, M_mars, M_saturn,
                             M_uranus, M_neptune, M_pluto};

    /*
    dt = 0.001;
    NumTimesteps = 250000;
    NumCelestials = sizeof(Celestial_masses)/sizeof(*Celestial_masses);

    New_system_and_solve(NumTimesteps, dt, Celestial_names, Celestial_masses, NumCelestials
                         , "SolarSys_All_planets.txt", "verlet");

    dt = 0.0001;
    NumTimesteps = 25000;
    New_system_and_solve(NumTimesteps, dt, Celestial_names, Celestial_masses, NumCelestials
                         , "SolarSys_All_planets_First4.txt", "verlet");
    */

    /*
    SolarSystem MercurySys;
    MercurySys.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1);
    vec3 Mercpos(0.3075, 0, 0);
    vec3 Mercvel(0, 12.44, 0);
    MercurySys.createCelestialBody(Mercpos, Mercvel, M_mercury);
    dt = 0.001;
    NumTimesteps = 100;
    for (int step=0; step<NumTimesteps; step++){
        MercurySys.write_file("Mercury_GR.txt");
        solver.Verlet_GR(MercurySys);
        plot_counter += 1;
    }
    */
    return 0;
}
