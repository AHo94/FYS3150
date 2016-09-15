#include <iostream>
#include <cmath>
#include <math.h>   // For mathematical expression, i.e exp
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
#include <time.h>
using namespace std;
using std::setw;

void fill_initial_arrays(double *x, double *a, double *b, double *c, double *f, int n)
{
    /* Function filling the initial arrays
    Will here assume that the values along the diagonal are the same
    to make sure that the algorithm works with the specialized case*/
    float L = 1;        // End point of x
    float h = L/(n-1);  // Steplength h
    for (int i=0; i<n; i++)
    {
        x[i] = i*(L/(n-1));
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        f[i] = h*h*100*exp(-10*x[i]);
    }
}

void forward_and_backward_subst(double *a, double *b, double *c, double *f, double *v, int n)
{
    /* Function solving forward and backward substitution
    Assuming different values along the diagonal of the matrix */
    for (int i=1; i<n; i++)
    {
        // Forward substitution
        b[i] = b[i] - (c[i-1]*a[i-1])/b[i-1];
        f[i] = (f[i] - (f[i-1]*a[i-1])/b[i-1]);
    }
    v[n-1] = f[n-1]/b[n-1];     // Initial (last) value of v
    for (int i=n-1; i>-1; i--)
    {
        // Backward substitution
        v[i-1] = (1/b[i-1])*(f[i-1] - c[i-1]*v[i]) ;
    }
}

void forward_simplified(double *x, double *b, double *f, double *v, int n)
{
    // Simplified algorithm for a special case where the values along the diagonal are the same
    float L = 1;            // End point of x
    float h = L/(n-1);      // Steplength h
    float float_converter = 1;  // Converts int to float value to prevent integer divison
    for (int i=0; i<n+1; i++)
    {
        x[i] = i*(L/(n-1));
        f[i] = h*h*100*exp(-10*x[i]);
        b[i] = float_converter*(i+1)/i;
    }

    for (int i=1; i<n; i++)
    {
        // Forward substitution
        //f[i] = f[i] + f[i-1]*float_converter*(i)/(i+1);
        f[i] = f[i] + f[i-1]/b[i];
    }
    //v[n-1] = f[n-1]/(L*n/(n-1));
    v[n-1] = f[n-1]/(b[n-1]);    // Initial (last) value of v
    for (int i=n-2; i>-1; i--)
    {
        // Backward substitution
        //v[i] = (float_converter*i/(i+1))*(f[i] + v[i+1]);
        v[i] = (f[i] + v[i+1])/b[i];
    }
}

void write_file(double *x, double *v, int n, string filename)
{
    // Function that writes data of x and v to a file
    ofstream datafile;
    datafile.open(filename);
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
    // TASK B)
    clock_t start, finish;
    start = clock();
    int n;                   // number of gridpoints
    double *x, *a, *b, *c, *f, *v;  // Pointer of the arrays
    string filename = "General_data_n";
    for (int i=0; i <= 3; i++)
    {
        /* For loop that solves the general method
        Uses values of n = 10, 100, 1000 */
        n = pow(10.0,i);
        x = new double[n];
        a = new double[n];
        b = new double[n];
        c = new double[n];
        f = new double[n];
        v = new double[n];

        string fileout = filename;
        string argument = to_string(n);
        fileout.append(argument);
        fileout.append(".txt");

        // Solving the algorithms and write results of x and v to a file
        fill_initial_arrays(x, a, b, c, f, n);
        forward_and_backward_subst(a, b, c, f, v, n);
        write_file(x, v, n, fileout);

    }

    // TASK C)
    // Freeing memory for next task
    delete[]a;
    delete[]c;
   // delete[]b;
    // Solving for specialized algorithm, with n = 10^6s
    /*
    n = pow(10,6);
    x = new double[n];
    f = new double[n];
    v = new double[n];
    b = new double[n];

    forward_simplified(x, b, f, v, n);
    write_file(x, v, n, "Project1c_data_simplified.txt");


//    // TASK D)
//    n = pow(10,6);
//    x = new double[n];
//    f = new double[n];
//    v = new double[n];
//    b = new double[n];

//    forward_simplified(x, f, v, n);
//    write_file(x, v, n, "Project1d_relative_error.txt");
    */
    cout << "Sucess!" << endl;
    delete [] x;
    delete [] f;
    delete [] v;
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << "s" << endl;;

    return 0;
}
