#include <iostream>
#include <math.h>   // For mathematical expression, i.e exp
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
using namespace std;
using std::setw;

double fill_initial_arrays(double *x, double *a, double *b, double *c, double *f, int n)
{
    // Function filling the initial arrays
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

double forward_subt(double *a, double *b, double *c, double *f, int n)
{
    // Function that uses forward subtitution to calculate new b and f
    for (int i=1; i<n; i++)
    {
        b[i] = b[i] - (c[i-1]*a[i])/b[i-1];
        f[i] = f[i] - f[i-1]*(a[i]/b[i-1]);
    }
}

double backward_subt(double *b, double *c, double *f, double *v, int n)
{
    // Function that uses backward subtitution to find new v
    v[n-1] = f[n-1]/b[n-1];     // Initial (last) value of v
    for (int i=n-2; i>-1; i--)
    {
        v[i] = 1/(b[i])*(f[i] - c[i]*v[i+1]);
    }
}

void write_file(double *x, double *v, int n)
{
    // Function that writes data to a file
    ofstream datafile;
    datafile.open("project_1_data.txt");
    datafile << "n = " << n << "\n";
    datafile << "x" << setw(15) << "v" << "\n";
    for (int i=0; i < n; i++)
    {
        datafile << x[i] << setw(15) << v[i] << "\n";
    }
    datafile.close();
}

int main()
{
    int n = 10;       // number of gridpoints
    double *x, *a, *b, *c, *f, *v;  // Pointer of the arrays
    // Creating new arrays
    x = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    v = new double[n];

    // Solving the algorithms and write results of x and v to a file
    fill_initial_arrays(x, a, b, c, f, n);
    forward_subt(a, b, c, f, n);
    backward_subt(b, c, f, v, n);
    write_file(x, v, n);

    cout << "Sucess!" << endl;
    return 0;
}
