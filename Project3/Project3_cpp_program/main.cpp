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

    /*
    // Adding Earth
    vec3 Earthpos (9.890331046925951E-01, 1.768079890757788E-01, -1.738715302893284E-04);
    vec3 Earthvel (-3.268395786841218E-03, 1.689265025904021E-02, -9.889230545039174E-07);
    Earthvel /= YrstoDays;
    System.createCelestialBody(Earthpos, Earthvel, M_earth);

    // Adding Jupiter
    vec3 Jupiterpos (-5.433468170028908E+00, -3.819061221110369E-01, 1.231004384238452E-01);
    vec3 Jupitervel (4.425651679847022E-04, -7.171108917491057E-03, 1.992744446163222E-05);
    Jupitervel /= YrstoDays;
    System.createCelestialBody(Jupiterpos, Jupitervel, M_jupiter);

    // Adding Mercury
    vec3 Mercurypos (-1.388351215994794E-01, 2.874076124640064E-01, 3.611730762400382E-02);
    vec3 Mercuryvel (-3.081033504804020E-02, -1.153752302730325E-02, 1.883146626624065E-03);
    Mercuryvel /= YrstoDays;
    System.createCelestialBody(Mercurypos, Mercuryvel, M_mercury);

    // Adding Mars
    vec3 Marspos (1.074137801923908E+00, -8.751565791508180E-01, -4.484330241084649E-02);
    vec3 Marsvel (9.406114898320066E-03, 1.202410268499505E-02, 2.097435381751857E-05);
    Marsvel /= YrstoDays;
    System.createCelestialBody(Marspos, Marsvel, M_mars);

    // Adding Saturn
    vec3 Saturnpos (-2.318303159285024E+00, -9.761896118742531E+00, 2.619996914174937E-01);
    vec3 Saturnvel (5.122777078109817E-03, -1.306326757884738E-03, -1.812965626845902E-04);
    Saturnvel /= YrstoDays;
    System.createCelestialBody(Saturnpos, Saturnvel, M_saturn);

    // Adding Uranus
    vec3 Uranuspos (1.847838443295973E+01, 7.526847462019028E+00, -2.114361038013451E-01);
    vec3 Uranusvel (-1.512383759608680E-03, 3.459146288519939E-03, 3.242416050249801E-05);
    Uranusvel /= YrstoDays;
    System.createCelestialBody(Uranuspos, Uranusvel, M_uranus);

    // Adding Neptune
    vec3 Neptunepos (2.825072704568992E+01, -9.952093577677799E+00, -4.461218547546795E-01);
    vec3 Neptunevel (1.022588623946866E-03,  2.979574756810357E-03, -8.514325155122657E-05);
    Neptunevel /= YrstoDays;
    System.createCelestialBody(Neptunepos, Neptunevel, M_neptune);

    // Adding Pluto
    vec3 Plutopos (9.393096450667111E+00, -3.182064102580347E+01, 6.879522592437006E-01);
    vec3 Plutovel (3.065499934972441E-03, 2.293283900283695E-04, -9.119583887771224E-04);
    Plutovel /= YrstoDays;
    System.createCelestialBody(Plutopos, Plutovel, M_pluto);
    */
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
