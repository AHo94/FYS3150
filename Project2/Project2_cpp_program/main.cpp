#include <iostream> // Printing purposes
//#include <cmath>
#include <math.h>   // Mathematics library
#include <algorithm>    // Used to sort arrays
//#include <numeric>
#include <fstream>  // Writing to file
#include <iomanip>  // setw identitation for output file
#include <sstream>  // Convert numbers to string within a set precision
#include <time.h>   // Clock
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

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

double max_offdiag(double **A, int n, double *max_diag_indices){
    // Finds the max value along the off diagonal
    double max_value=0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            double a_ij = fabs(A[i][j]);
            if (a_ij > max_value){
                max_value = a_ij;
                max_diag_indices[0] = i;
                max_diag_indices[1] = j;
            }
        }
    }
    return max_value;
}

void max_diag_testing(double **max_diag_test_matrix, int n, int i, int j){
    /* Unit test to see if the max_diag function gives the largest (absolute value) non diagonal matrix element.
     * Test matrix is a tridiagonal matrix with 2 along the diagonal and -1 along the non diagonal.
     * One of the non diagonal element is to be selected to be larger than |-1|, e.g -1000.
     */
    double *index_testing, *input_index_test;
    input_index_test = new double[2];
    index_testing = new double[2];
    input_index_test[0] = i;
    input_index_test[1] = j;

    // Fills the test matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j) {
                max_diag_test_matrix[i][j] = 2;
            } else if (i==j-1) {
                max_diag_test_matrix[i][j] = -1.0;
            } else if (i==j+1) {
                max_diag_test_matrix[i][j] = -1.0;
            } else{
                max_diag_test_matrix[i][j] = 0;
            }
        }
    }
    // One of the off diagonals to become the largest value
    max_diag_test_matrix[i][j] = -1000;
    cout << "Max (absolute) value in the test matrix: " << fabs(max_diag_test_matrix[i][j]) << endl;
    double test_max_diag = max_offdiag(max_diag_test_matrix, n, index_testing);
    if (fabs(fabs(test_max_diag) - fabs(max_diag_test_matrix[i][j])) < 1e-8){
        cout << "max_diag does give out the largest value" << endl;
        cout << "max_diag result: " << test_max_diag << " \n" << endl;
    }
    else{
        cout << "max_diag test did not work \n" << endl;
    }
    // Freeing memory
    delete[]index_testing;
    delete[]input_index_test;
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

void armadillo_test(double *d, double rho_max, int n){
    double h = (rho_max)/(n+1);
    mat Arma_test;
    Arma_test.zeros(n,n);
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j) {
                Arma_test(i,j) = d[i];
            } else if (i==j-1) {
                Arma_test(i,j) = -1.0/(h*h);
            } else if (i==j+1) {
                Arma_test(i,j) = -1.0/(h*h);
            } else{
                Arma_test(i,j) = 0;
            }
        }
    }
    vec eigval = eig_sym(Arma_test);
    cout << "Lowest 3 eigenvalues of the Armadillo simulation: " << endl;
    for (int i=0; i<3; i++)
        cout << eigval(i) << endl;
}

void write_file(double **R, double *rho, double *V, double lambda
                , double rho_max, double omega, int n, int *vector_index, string filename){
    // Function that writes eigenvector and rho data to an output file.
    ofstream datafile;
    datafile.open(filename);
    int min1 = vector_index[0];
    datafile << "# First row contains the n, rho_max and omega values respectively.";
    datafile << " The other rows contains the data to be plotted \n";
    datafile << "# First column contains rho values. Second column contains the eigenvector for the ground state"
                "third column contains the potential at point i. Fourth column contains the eigenvalues. \n";
    datafile << n << setw(15) << rho_max << setw(15) << omega << "\n";
    for (int i=0; i<n; i++){
        datafile << rho[i] << setw(15)
                 << pow(R[i][min1], 2) << setw(15)
                 << V[i] << setw(15) << lambda << "\n";
    }
    datafile.close();
}

int main(){
    clock_t start, finish;
    double *d, *rho, **A, **R, *max_diag_indices;
    double **max_diag_test_matrix;
    //double *index_testing, *input_index_test;
    int n = 400;
    double rho_max = 10;

    cout << "Using a " << n << "x" << n << " matrix (n = "<< n << ")" << endl;
    cout << "with rho_max = " << rho_max << "\n" << endl;
    // Initialize everything
    d = new double[n];
    rho = new double[n];
    A = new double*[n];
    R = new double*[n];
    max_diag_indices = new double[2];
    max_diag_test_matrix = new double*[n];
    for (int i=0; i<n; i++){
        A[i] = new double[n];
        R[i] = new double[n];
        max_diag_test_matrix[i] = new double[n];
    }
    // Unit test for max_offidag function
    max_diag_testing(max_diag_test_matrix, n, 1, 2);
    // Freeing memory
    for (int i=0; i<n; i++){
        delete[]max_diag_test_matrix[i];
    }
    delete[]max_diag_test_matrix;

    cout << "Starting Jacobi's algorithm for a single electron" << endl;
    // --- Non interacting case ---
    start = clock();
    initialize_matrix(A, R, d, rho, rho_max, n);
    double max_diag = 1;
    int iterations = 0;
    int maxiter = 400000;
    double tolerance = 1.0e-8;
    while (max_diag > tolerance && iterations <= maxiter){
        max_diag = max_offdiag(A, n, max_diag_indices);
        Jacobi_rotation(A, R, max_diag_indices[0], max_diag_indices[1], n);
        iterations++;
    }
    finish = clock();
    cout << "Time elapsed for non interactive case: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

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

    // --- Armadillo comparison ---
    cout << "Using Armadillo to compare the time\n" << endl;
    start = clock();
    armadillo_test(d, rho_max, n);
    finish = clock();
    cout << "Time elapsed for Armadillo's algorithm: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;
    //return 0;

    // --- Interacting case ---
    cout << "Calculating interacting case... \n" << endl;

    double omegas[] = {0.01, 0.5, 1, 5};
    int *vector_index;   // An array that stores the indices from the matrix R where the smallest eigenvector is
    double *V;
    double specified_rho_max[] = {40, 5, 4, 2};
    vector_index = new int[1];
    string filename = "Eigenvector_data_omega_";
    for (int i=0; i<4; i++){
        cout << "Calculating for the case with omega = " << omegas[i] << endl;
        start = clock();
        d = new double[n];
        rho = new double[n];
        A = new double*[n];
        R = new double*[n];
        V = new double[n];
        max_diag_indices = new double[n];
        for (int j=0; j<n; j++){
            A[j] = new double[n];
            R[j] = new double[n];
        }

        initialize_matrix(A, R, d, rho, specified_rho_max[i], n, omegas[i]);

        for (int j=0; j<n; j++){
            // Filling out the potential
            V[j] = pow(omegas[i]*rho[j],2) + 1.0/rho[j];
        }
        max_diag = 1;
        iterations = 0;
        while (max_diag > tolerance && iterations <= maxiter){
            max_diag = max_offdiag(A, n, max_diag_indices);
            Jacobi_rotation(A, R, max_diag_indices[0], max_diag_indices[1], n);
            iterations++;
        }
        finish = clock();
        cout << "Time elapsed: " << ((finish-start)/(double)(CLOCKS_PER_SEC)) << "s" << endl;

        lambda = new double[n];
        cout <<"Number of iterations: " << iterations << endl;
        for (int j=0; j<n; j++){
            lambda[j] = A[j][j];
        }

        // Sorting eigenvalues from lowest to highest
        std::sort(lambda, lambda+n);
        cout << "Lowest 3 eigenvalues are: " << endl;
        for (int j=0; j<3; j++){
            cout << lambda[j] << endl;
        }
        Smallest_eigenvector(A, lambda, vector_index, n);
        orthogonal_test(R, n);

        string fileout = filename;
        stringstream stream;
        stream << fixed << setprecision(2) << omegas[i];
        string argument = stream.str();
        fileout.append(argument);
        fileout.append(".txt");
        write_file(R, rho, V, lambda[0], specified_rho_max[i], omegas[i], n, vector_index, fileout);
    }
    return 0;
}
