#include <iostream>
#include <cmath>
#include <math.h>
//#include <armadillo>


using namespace std;
//using namespace arma;

void initialize_matrix(double **A, double **R, double *d, double *rho, int n){
    /* This function initializes the initial matrix A for the project.
       Does not calculate the values along the off-diagonal as they are the same.
    */
    double h = (50.0/(n+1));
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
                R[i][j] = 1;
            } else if (i==j-1) {
                A[i][j] = -1.0/(h*h);
                R[i][j] = 0;
            } else if (i==j+1) {
                A[i][j] = -1.0/(h*h);
                R[i][j] = 0;
            } else{
                A[i][j] = 0;
                R[i][j] = 0;
            }
        }
    }
}

double max_offdiag(double **A, int p, int q, int n){
    // Finds the max value along the off diagonal
    double max_value;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            double a_ij = fabs(A[i][j]);
            if (a_ij > max_value){
                max_value = a_ij; p = i, q = j;
            }
        }
    }
    return max_value;
}

void Jacobi_rotation(double **A, double **R, int k, int l, int n){
    double s, c; // Sine and cosine functions
    if (A[k][l] != 0.0) {
        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);

        if (tau >= 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1+t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    A[k][k] = c*c*a_kk - 2*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0;
    A[l][k] = 0;
    for (int i=0; i<n; i++){
        if (i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    return;
}

int main(){
    //double A;
    //A = mat (1,2);
    double tolerance = 1.0e-8;
    double *d, *rho, **A, **R;
    int n = 5;
    d = new double[n];
    rho = new double[n];
    A = new double*[n];
    R = new double*[n];
    for (int i=0; i<n; i++){
        A[i] = new double[n];
        R[i] = new double[n];
    }

    initialize_matrix(A, R, d, rho, n);
    /*
    for (int k=0; k<n; k++){
        cout << A[k][0] << " " << A[k][1] << " " << A[k][2] << " " << A[k][3] << " " << A[k][4] << endl;
    }
    cout << "\n" << endl;

    for (int k=0; k<n; k++){
        cout << R[k][0] << " " << R[k][1] << " " << R[k][2] << " " << R[k][3] << " " << R[k][4] << endl;
    }
    cout << "pls" << endl;
    */
   // delete[]d;
   // delete[]rho;
    int p, q;
    double max_diag;
    int iterations = 0;
    int maxiter = 1000;
    max_diag = 1;

    while (max_diag > tolerance && iterations <= maxiter){
        max_diag = max_offdiag(A, p, q, n);
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (max_diag == fabs(A[i][j]))
                {
                    p=i, q=j;
                }
            }
        }
        /*
        for (int k=0; k<n; k++){
            cout << A[k][0] << " " << A[k][1] << " " << A[k][2] << " " << A[k][3] << " " << A[k][4] << endl;
        }
        cout << "\n" << endl;
        */
        //cout << p << " "<<q<<endl;
        Jacobi_rotation(A, R, p, q, n);
        iterations++;
    }
    cout <<"Number of iterations: " << iterations << endl;


    //delete[]d;
    //delete[]e;
    //delete[]rho;
    //delete[]A;
    return 0;
}
