#include <iostream>
#include <cmath>
#include <math.h>
#include <algorithm>
//#include <armadillo>


using namespace std;
//using namespace arma;

void initialize_matrix(double **A, double **R, double *d, double *rho, double rho_max, int n){
    /* This function initializes the initial matrix A and R for the project.
       Does not calculate the values along the off-diagonal as they are the same.
       Matrix R is the identity matrix
    */
    double h = (rho_max/(n+1));
    for (int i=0; i<n; i++){
        // Creates the rho gridpoints
        rho[i] = (i+1)*h;
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
    double max_value=0;
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

void orthogonal_test(double **R, int n){
    double *w1, *w2;
    w1 = new double[n];
    w2 = new double[n];
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            w1 = R[i];
            w2 = R[j];
            if (fabs(w1*w2 - 1) < 1e-10){
                cout << "Orthogonality conserved" << endl;
            }
            else{
                cout << "Orthogonality not conserved" << endl;
            }
        }
    }
    delete[]w1;
    delete[]w2;
}

void Jacobi_rotation(double **A, double **R, int k, int l, int n){
    // Jacobi rotation algorithm
    double s, c; // Sine and cosine functions
    if (A[k][l] != 0.0) {
        double t=0, tau=0;
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
    // Doing new rotation to get new matrix A
    double a_kk=0, a_ll=0, a_ik=0, a_il=0, r_ik=0, r_il=0;
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
        // Saving eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    orthogonal_test(R, n);
    return;
}


int main(){
    //double A;
    //A = mat (1,2);
    double tolerance = 1.0e-8;
    double *d, *rho, **A, **R;
    int n = 30;
    double rho_max = 40;
    d = new double[n];
    rho = new double[n];
    A = new double*[n];
    R = new double*[n];
    for (int i=0; i<n; i++){
        A[i] = new double[n];
        R[i] = new double[n];
    }

    initialize_matrix(A, R, d, rho, rho_max, n);

    // Deleting unused arrays
    delete[]d;
    delete[]rho;

    double max_diag = 1;
    int iterations = 0;
    int maxiter = 5000;
    while (max_diag > tolerance && iterations <= maxiter){
        int p = 0;
        int q = 0;
        max_diag = max_offdiag(A, p, q, n);
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (fabs(max_diag - fabs(A[i][j])) < tolerance){
                    p=i;
                    q=j;
                }
            }
        }
        Jacobi_rotation(A, R, p, q, n);
        iterations++;
    }

    double *lambda;
    lambda = new double[n];
    cout <<"Number of iterations: " << iterations << endl;
    for (int i=0; i<n; i++){
        lambda[i] = A[i][i];
    }
    // Sorting eigenvalues from lowest to highest
    std::sort(lambda, lambda+n);
    cout << "Lowest 5 eigen values is: " << endl;
    for (int i=0; i<5; i++){
        // Printing the 5 smallest eigenvalues
        cout << lambda[i] << endl;
    }

    return 0;
}
