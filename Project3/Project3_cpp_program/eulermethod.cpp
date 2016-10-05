#include "eulermethod.h"
#include <cmath>
#include <iostream>
using namespace std;
eulermethod::eulermethod()
{
}

void eulermethod::solve_euler_2D(double **velocity, double **position, int N){
    double h = 1.0/(N+1);
    double pi = atan(1)*4;  // Define pi
    for (int i=0; i<N-1; i++){
        double r_i = sqrt(pow(position[0][i],2) + pow(position[1][i],2));
        velocity[0][i+1] = velocity[0][i] - h*4*pow(pi,2)*position[0][i]/pow(r_i,3);
        velocity[1][i+1] = velocity[1][i] - h*4*pow(pi,2)*position[1][i]/pow(r_i,3);
        position[0][i+1] = position[0][i] + h*velocity[0][i];
        position[1][i+1] = position[1][i] + h*velocity[1][i];
    }
}
