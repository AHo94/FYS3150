#include <iostream>
#include <cmath>
#include <math.h>   // For mathematical expression, i.e exp
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
using namespace std;
using std::setw;

void fill_initial_arrays(double *x, double *a, double *b, double *c, double *f, int n)
{
    /* Function filling the initial arrays
    Will here assume that the values along the diagonal are the same */
    float L = 1;        // End point of x
    float h = L/(n-1);  // Value of h
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
    float L = 1;
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

void forward_simplified(double *x, double *f, double *v, int n)
{
    // Simplified algorithm for a special case where the values along the diagonal are the same
    float L = 1;
    float h = L/(n-1);
    float float_converter = 1;  // Converts int to float value to prevent integer divison
    for (int i=0; i<n+1; i++)
    {
        x[i] = i*(L/(n-1));
        f[i] = h*h*100*exp(-10*x[i]);
    }

    for (int i=1; i<n; i++)
    {
        f[i] = f[i] + f[i-1]*float_converter*(i)/(i+1);
    }
    // Initial (last) value of v
    v[n-1] = f[n-1]*(L*n/(n-1));

    for (int i=n-2; i>-1; i--)
    {
       v[i] = (float_converter*i/(i+1))*(f[i] + v[i+1]);
    }
    cout << v[0] << endl;

}

void write_file(double *x, double *v, int n, string filename)
{
    // Function that writes data to a file
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
    //int n = 10;       // number of gridpoints
    double *x, *a, *b, *c, *f, *v;  // Pointer of the arrays
    int filename_index = 0;     // Used to select filenames
    // Array with filenames for the general case
    char *filenames[] = {"Project1_data_n10.txt", "Project1_data_n100.txt", "Project1_data_n1000.txt"};
    for (int n = 10; n < 1001; n = 10*n)
    {
        /* For loop runs through n = 10, 100, 100
        Creating new arrays for new n */
        x = new double[n];
        a = new double[n];
        b = new double[n];
        c = new double[n];
        f = new double[n];
        v = new double[n];

        // Solving the algorithms and write results of x and v to a file
        fill_initial_arrays(x, a, b, c, f, n);
        forward_and_backward_subst(a, b, c, f, v, n);
        write_file(x, v, n, filenames[filename_index]);
        filename_index = filename_index + 1;

    }
    /*
    // Creating new arrays
    x = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    v = new double[n];

    // Solving the algorithms and write results of x and v to a file
    fill_initial_arrays(x, a, b, c, f, n);
    forward_and_backward_subst(a, b, c, f, v, n);
    write_file(x, v, n, "Project1_data_n10.txt");

    // Increase number of points to 100 and do same calculations
    n = 100;
    x = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    v = new double[n];

    fill_initial_arrays(x, a, b, c, f, n);
    forward_and_backward_subst(a, b, c, f, v, n);
    write_file(x, v, n, "Project1_data_n100.txt");

    // Increase n to 1000
    n = 1000;
    x = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    v = new double[n];

    fill_initial_arrays(x, a, b, c, f, n);
    forward_and_backward_subst(a, b, c, f, v, n);
    write_file(x, v, n, "Project1_data_n1000.txt");
    */

    // TASK C)
    // Freeing memory for next task
    delete[]a;
    delete[]c;
    delete[]b;
    delete[]filenames;
    // Solving for specialized algorithm, with n = 10^6

    int n = pow(10,6);
    x = new double[n];
    f = new double[n];
    v = new double[n];

    forward_simplified(x, f, v, n);
    write_file(x, v, n, "Project1c_data_simplified.txt");

    cout << "Sucess!" << endl;
    delete [] x;
    delete [] f;
    delete [] v;


    return 0;
}
