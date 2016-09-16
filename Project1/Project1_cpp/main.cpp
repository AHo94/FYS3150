#include <iostream>
#include <cmath>
#include <math.h>   // For mathematical expression, i.e exp
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
#include <time.h>
#include <lib.h>
using namespace std;
using std::setw;

double f_tild(double x, float h){
    // Calculates the function f_tilde for a given x value
    return h*h*100*exp(-10*x);
}

void fill_initial_arrays(double *x, double *a, double *b, double *c, double *f, int n, float L){
    /* Function filling the initial arrays
    Will here assume that the values along the diagonal are the same
    to make sure that the algorithm works with the specialized case*/
    float h = L/(n+1);  // Steplength h
    for (int i=0; i<n; i++){
        x[i] = (i+1)*h;
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        //f[i] = h*h*100*exp(-10*x[i]);
        f[i] = f_tild(x[i], h);
    }
}

void forward_and_backward_subst(double *a, double *b, double *c, double *f, double *v, int n){
    /* Function solving forward and backward substitution
    Assuming different values along the diagonal of the matrix */
    for (int i=1; i<n; i++){
        // Forward substitution
        b[i] = b[i] - (c[i-1]*a[i-1])/b[i-1];
        f[i] = (f[i] - (f[i-1]*a[i-1])/b[i-1]);
    }
    v[n-1] = f[n-1]/b[n-1];     // Initial (last) value of v
    for (int i=n-1; i>0; i--){
        // Backward substitution
        v[i-1] = (1/b[i-1])*(f[i-1] - c[i-1]*v[i]) ;
    }
}

void simplified_algorithm(double *x, double *b, double *f, double *v, int n, float h){
    // Simplified algorithm for a special case where the values along the diagonal are the same
    //float h = L/(n+1);
    float float_converter = 1;  // Converts int to float value to prevent integer divison
    for (int i=0; i<n; i++){
        x[i] = (i+1)*h;
        //f[i] = h*h*100*exp(-10*x[i]);
        f[i] = f_tild(x[i], h);
        b[i] = float_converter*(i+1)/i;
    }

    for (int i=1; i<n; i++){
        // Forward substitution
        f[i] = f[i] + f[i-1]/b[i];
    }
    v[n-1] = f[n-1]/(b[n-1]);    // Initial (last) value of v
    for (int i=n-2; i>-1; i--){
        // Backward substitution
        v[i] = (f[i] + v[i+1])/b[i];
    }
}

void write_file(double *x, double *v, int n, string filename){
    // Function that writes data of x and v to a file
    ofstream datafile;
    datafile.open(filename);
    int max_points = 1000;  // Saving up to 1000 points
    datafile << "# First line is the value of n \n" ;
    datafile << n << "\n";
    int step;
    if(n > max_points){
        /*
        Creating an if-test to reduce the number of points saved.
        Saving at most 1000 points, if we have more points, we will do n/max_points jumps
        between the intervals. e.g: n = 10^4, then we save every 10 points.
        */
        step = n/max_points;
    }
    else{
        //  Saves every point if n < max_points
        step = 1;
    }
    for (int i=0; i < n; i=i+step){
        datafile << x[i] << setw(15)  << v[i] << "\n";
    }
    datafile.close();
}

void create_tridiagonal_matrix(double **A, int n){
    // Function that fills a given nxn matrix and fills the diagonals
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i==j) {
                A[i][j] = 2;
            } else if (i==j-1) {
                A[i][j] = -1;
            } else if (i==j+1) {
                A[i][j] = -1;
            } else{
                A[i][j] = 0;
            }
        }
    }
}


// Run main program
int main()
{

//    // TASK B) - General Algorithm
    clock_t start, finish;
    int n;                                // number of gridpoints
    float L = 1;                          // Endpoint of x
    double *x, *a, *b, *c, *f, *v;        // Pointers for array
    start = clock();
    string filename = "General_data_n";   // Filename of our general algorithm
    for (int i=1; i <= 3; i++){
        /* For loop that solves the general method
        Uses values of n = 10, 100, 1000 */
        n = (int) pow(10.0,i);
        x = new double[n];
        a = new double[n];
        b = new double[n];
        c = new double[n];
        f = new double[n];
        v = new double[n];

        // Adds something extra to the filename to distinguis between the files
        string fileout = filename;
        string argument = to_string(n);
        fileout.append(argument);
        fileout.append(".txt");

        // Solving the algorithms and write results of x and v to a file
        fill_initial_arrays(x, a, b, c, f, n, L);
        forward_and_backward_subst(a, b, c, f, v, n);
        write_file(x, v, n, fileout);
    }
    finish = clock();
    cout << "Time elapsed for general algorithm: " << ((finish-start)/(double)(CLOCKS_PER_SEC)/1000) << "s" << endl;

    // TASK C) - Simplified algorithm
    // Freeing memory for next task
    delete[]a;
    delete[]c;
    start = clock();
    string filename_simplified = "Simplified_data_n";   // Filename for simplified algorithm
    for (int i=1; i <= 6; i++){
        // For loop that runs through exponents from i=1 to i=6
        n = pow(10,i);
        x = new double[n];
        f = new double[n];
        v = new double[n];
        b = new double[n];

        // Adds something extra to the filename to distinguis between the files
        string fileout = filename_simplified;
        string argument = to_string(n);
        fileout.append(argument);
        fileout.append(".txt");

        if (i==3){
            // Stops clock for the specialized algorithm to compare the CPU time with the general algorithm.
            finish = clock();
            cout << "Time elapsed for specialized algorithm: "
                 << ((finish-start)/double(CLOCKS_PER_SEC)/1000) << "s" << endl;
        }

        // Solving the algorithms and write results of x and v to a file
        simplified_algorithm(x, b, f, v, n, L/(n+1));
        write_file(x, v, n, fileout);
    }
    //finish = clock();
    //cout << "Time elapsed for specialized algorithm: " << ((finish-start)/double(CLOCKS_PER_SEC)/1000) << "s" << endl;
    // TASK D) - Calculate relative error
    string filename_error = "Error_data_n";     // Filename for relative error data
    for (int i=1; i <= 7; i++){
        // For loop that runs through the exponents from i=1 to i=7
        n = pow(10,i);
        x = new double[n];
        f = new double[n];
        v = new double[n];
        b = new double[n];

        // Adds something extra to the filename to distinguis between the files
        string fileout = filename_error;
        string argument = to_string(n);
        fileout.append(argument);
        fileout.append(".txt");

        /* Solving the algorithms and write results of x and v to a file
           Using simplified algorithm */
        simplified_algorithm(x, b, f, v, n, L/(n+1));
        write_file(x, v, n, fileout);
    }

    // TASK E) - LU-decomposition
    double **A, d;
    start = clock();
    string filename_LUD = "LUDecomp_data_n";
    for (int i=1; i<=3; i++){
        n = pow(10,i);
        A = new double*[n];
        for (int i=0; i<n; i++) {
            A[i] = new double[n];
        }
        create_tridiagonal_matrix(A, n);
        int index[n];
        //double d;
        x = new double[n];
        f = new double[n];
        for (int i=0; i<n; i++){
            x[i] = (i+1)*(L/(n+1));
            f[i] = f_tild(x[i], L/(n+1));
        }

        /* Uses lib.cpp to compute LU-decompositiion
           Results are overwritten in the f variable */
        ludcmp(A,n,index,&d);
        lubksb(A,n,index, f);

        // Giving filename a specific name
        string fileout = filename_LUD;
        string argument = to_string(n);
        fileout.append(argument);
        fileout.append(".txt");

        write_file(x, f, n, fileout);
    }
    finish = clock();
    cout << "Time elapsed for LU-decomp: " << ((finish-start)/(CLOCKS_PER_SEC)) << "s" << endl;
    cout << "Program finished" << endl;
    delete [] x;
    delete [] f;
    delete [] v;
    delete [] b;
    delete [] A;
    return 0;
}
