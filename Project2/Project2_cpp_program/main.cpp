#include <iostream>
#include <cmath>
#include <math.h>
//#include <armadillo>


using namespace std;
//using namespace arma;

void initialize_matrix(double **A, double *d, double *e, double *rho, int n){
    /* This function initializes the initial matrix A for the project.
       Does not calculate the values along the off-diagonal as they are the same.
    */
    double h = (4.0/(n+1));
    for (int i=0; i<n; i++){
        // Creates the rho gridpoints
        rho[i+1] = rho[i] + i*h;
    }

    for (int i=0; i<n; i++){
        // Calculates the values along the diagonal
        d[i] = 2.0/(h*h) + pow(rho[i], 2);
    }

    // Starts filling the matrix elements.
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j) {
                A[i][j] = d[i];
            } else if (i==j-1) {
                A[i][j] = -1.0/(h*h);
            } else if (i==j+1) {
                A[i][j] = -1.0/(h*h);
            } else{
                A[i][j] = 0;
            }
        }
    }
}

double max_offdiag(double **A, int p, int q, int n){
    // Finds the max value along the off diagonal
    double max_value;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            double a_ij = fabs(A[i][j]);
            if (a_ij > max_value){
                max_value = a_ij; p = i, q = j;
            }
        }
    }
    return max_value;
}

void Jacobi_rotation(double **a, int k, int l, int n){
    double s, c; // Sine and cosine functions
    if (A[k][l] != 0.0) {
        double t, tau;
        tau = (A[l][l] - A[k][k]/(2*A[k][l]));

        if (tau >= 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }
}

int main(){
    //double A;
    //A = mat (1,2);
    double tolerance = 1.0e-8;
    double *d, *e, *rho, **A;
    int n = 5;
    d = new double[n];
    e = new double[n];
    rho = new double[n];
    A = new double*[n];
    for (int i=0; i<n; i++){
        A[i] = new double[n];
    }

    initialize_matrix(A, d, e, rho, n);

    for (int k=0; k<n; k++){
        cout << A[k][0] << " " << A[k][1] << " " << A[k][2] << " " << A[k][3] << " " << A[k][4] << endl;
    }
    cout << "pls" << endl;
    int p, q;
    double max_diag = max_offdiag(A, p, q, n);
    cout << max_diag << endl;
    //delete[]d;
    //delete[]e;
    //delete[]rho;
    //delete[]A;
    return 0;
}
