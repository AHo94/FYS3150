#include <iostream>
#include <cmath>
#include <math.h>
#include <algorithm>    // Used to sort arrays
#include <numeric>
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
#include <sstream>
#include <time.h>

using namespace std;

void initialize_matrix(double **A, double **R, double *d, double *rho, double rho_max, int n, double omega=0){
    /* This function initializes the initial matrix A and R for the project.
     * Does not calculate the values along the off-diagonal as they are the same.
     * Matrix R is the identity matrix.
     * Omega is an optional argument. If omega = 0 we look at the non interacting case,
     * if omega != 0 we look at the interactng case.
    */
    double h = (rho_max/(n+1));
    for (int i=0; i<n; i++){
        // Creates the rho gridpoints
        rho[i] = (i+1)*h;
    }
    // Calculates the values along the diagonal, depending if we have interacting or non interacting case
    if (omega == 0){
        // For non interacting case, single electron.
        for (int i=0; i<n; i++){
            d[i] = 2.0/(h*h) + pow(rho[i], 2);
        }
    }
    else{
        // For interacting case, two electrons.
        for (int i=0; i<n; i++){
            d[i] = 2.0/(h*h) + omega*omega*pow(rho[i], 2) + 1.0/rho[i];
        }
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

double max_offdiag(double **A, int n){
    // Finds the max value along the off diagonal
    double max_value=0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            double a_ij = fabs(A[i][j]);
            if (a_ij > max_value){
                max_value = a_ij;
            }
        }
    }
    return max_value;
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
    return;
}

void orthogonal_test(double **R, int n){
    /* Function that checks if the final eigenvectors are orthogonal.
     * Has a small tolerance in case the values in the eigenvectors are very small.
     */
    double *w1, *w2;
    w1 = new double[n];
    w2 = new double[n];
    int N_orthogonal = 0;
    int N_non_orthogonal = 0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            w1 = R[i];
            w2 = R[j];
            if (fabs(std::inner_product(w1, w1+n, w2, 0)) < 1e-10){
                N_orthogonal++;
            }
            else{
                N_non_orthogonal++;
            }
        }
    }
    cout << "Possible number of orthogonal vectors: " << 0.5*n*(n-1) << endl;
    cout << "Number of orthogonal vectors: " << N_orthogonal << endl;
    cout << "Number of non-orthogonal vectors: " << N_non_orthogonal << "\n" <<endl;
    delete[]w1;
    delete[]w2;
}

void Smallest_eigenvector(double **A, double *lambda, int *vector_index, int n){
    /* This function finds the matrix index of the smallest eigenvector.
     * Compares the eigenvalues of the sorted eigenvalue array and the eigenvalue matrix.
     * Gives the index if the values are the same, or within a small tolerance.
    */
    double tolerance = 1.0e-10;
    for (int j=0; j<2; j++){
        for (int i=0; i<n; i++){
            if (fabs(lambda[j] - A[i][i]) < tolerance){
                vector_index[j] = i;
            }
        }
    }
}

void write_file(double **R, double *rho, double rho_max, double omega, int n, int *vector_index, string filename){
    // Function that writes eigenvector and rho data to an output file.
    ofstream datafile;
    datafile.open(filename);
    int min1 = vector_index[0];
    datafile << "# First row contains the n, rho_max and omega values respectively.";
    datafile << " The other rows contains the data to be plotted \n";
    datafile << "# First column contains rho values. Second column contains the eigenvector for the ground state \n";
    datafile << n << setw(15) << rho_max << setw(15) << omega << "\n";
    for (int i=0; i<n; i++){
        datafile << rho[i] << setw(15)
                 << pow(R[i][min1], 2) << "\n";
    }
    datafile.close();
}

int main(){
    clock_t start, finish;
    double *d, *rho, **A, **R;
    int n = 100;
    double rho_max = 10.0;

    cout << "Using a " << n << "x" << n << " matrix (n = "<< n << ")" << endl;
    cout << "with rho_max = " << rho_max << "\n" << endl;
    start = clock();
    d = new double[n];
    rho = new double[n];
    A = new double*[n];
    R = new double*[n];
    for (int i=0; i<n; i++){
        A[i] = new double[n];
        R[i] = new double[n];
    }

    initialize_matrix(A, R, d, rho, rho_max, n);
    double max_diag = 1;
    int iterations = 0;
    int maxiter = 200000;
    double tolerance = 1.0e-8;
    while (max_diag > tolerance && iterations <= maxiter){
        int p = 0;
        int q = 0;
        max_diag = max_offdiag(A, n);

        // For loop that finds the index of the max_diagonal value in the matrix A.
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
    finish = clock();
    cout << "Time elapsed for non interactive case: " << ((finish-start)/(double)(CLOCKS_PER_SEC)/1000) << "s" << endl;

    double *lambda;
    lambda = new double[n];
    cout <<"Number of iterations: " << iterations << endl;
    for (int i=0; i<n; i++){
        lambda[i] = A[i][i];
    }
    // Sorting eigenvalues from lowest to highest
    std::sort(lambda, lambda+n);
    cout << "Lowest 3 eigenvalues are: " << endl;
    for (int i=0; i<3; i++){
        // Printing the 3 smallest eigenvalues
        cout << lambda[i] << endl;
    }
    orthogonal_test(R, n);
    
    // --- Interacting case ---
    cout << "Calculating interacting case... \n" << endl;

    double omegas[] = {0.01, 0.5, 1, 5};
    int *vector_index;   // An array that stores the indices from the matrix R where the smallest eigenvector is
    vector_index = new int[1];
    string filename = "Eigenvector_data_omega_";
    for (int i=0; i<4; i++){
        cout << "Calculating for the case with omega = " << omegas[i] << endl;
        start = clock();
        d = new double[n];
        rho = new double[n];
        A = new double*[n];
        R = new double*[n];
        for (int i=0; i<n; i++){
            A[i] = new double[n];
            R[i] = new double[n];
        }
        initialize_matrix(A, R, d, rho, rho_max, n, omegas[i]);

        max_diag = 1;
        iterations = 0;
        while (max_diag > tolerance && iterations <= maxiter){
            int p = 0;
            int q = 0;
            max_diag = max_offdiag(A, n);
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
        finish = clock();
        cout << "Time elapsed: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

        lambda = new double[n];
        cout <<"Number of iterations: " << iterations << endl;
        for (int i=0; i<n; i++){
            lambda[i] = A[i][i];
        }

        // Sorting eigenvalues from lowest to highest
        std::sort(lambda, lambda+n);
        cout << "Lowest 3 eigenvalues are: " << endl;
        for (int i=0; i<3; i++){
            cout << lambda[i] << endl;
        }
        Smallest_eigenvector(A, lambda, vector_index, n);
        orthogonal_test(R, n);

        string fileout = filename;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        write_file(R, rho, rho_max, omegas[i], n, vector_index, fileout);
    }
    return 0;
}
