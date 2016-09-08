#include <iostream>
#include <vector>
#include <math.h>   // For mathematical expression, i.e exp
#include <fstream>  // Writing to file
#include <iomanip>
using namespace std;
using std::setw;

double fill_arrays(double *x, double *a, double *b, double *c, double *f, int n)
{
    float L = 1;        // End point of x
    float h = L/(n-1);  // Value of h
    for (int i=0; i<n; i++)
    {
        x[i] = i*(L/(n-1));
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        f[i] = h*h*100*exp (x[i]);
    }
}

float new_v(float b, float c, float f,float v)
{
    // Function that calculates new v values
    return (1/b)*(f - c*v);
}
int main()
{
    int n = 10;         // number of gridpoints
    //float L = 1;        // End point of x
    //float h = L/(n-1);  // Value of h

    //double *x[n], *a[n], *b[n], *c[n], *f[n], *v[n];  // Vector with n elements
    double *x, *a, *b, *c, *f, *v;  // Pointer of the arrays
    // Creating new arrays
    x = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    v = new double[n];

    // Filling out gridpoints for x, based on number of points
    fill_arrays(x, a, b, c, f, n);
//    for (int i=0; i<n; i++)
//    {
//        x[i] = i*(L/(n-1));
//        a[i] = -1;
//        b[i] = 2;
//        c[i] = -1;
//        f[i] = h*h*100*exp (x[i]);
//    }

    // Calculate the temporary values of b and f
    for (int i=1; i<n; i++)
    {
        b[i] = b[i] - (c[i-1]*a[i])/b[i-1];
        f[i] = f[i] - f[i-1]*(a[i]/b[i-1]);
    }
    // Calculate v_i values
    v[n-1] = f[n-1]/b[n-1];
    for (int i=n-2; i>-1; i--)
    {
        //v[i] = 1/(b[i])*(f[i] - c[i]*v[i+1]);
        v[i] = new_v(b[i], c[i], f[i], v[i+1]);
    }

    // Writing data to file
    ofstream datafile;
    datafile.open("project_1_data.txt");
    datafile << "n = " << n << "\n";
    datafile << "x" << setw(15) << "v" << "\n";// << setw(15) << "f" << setw(15) << "b"  <<  "\n";
    for (int i=0; i < n; i++)
    {
        datafile << x[i] << setw(15) << v[i] << "\n";// << setw(15) << f[i] << setw(15) << b[i] << "\n";
    }
    datafile.close();

    cout << "Sucess!" << endl;
    return 0;
}
