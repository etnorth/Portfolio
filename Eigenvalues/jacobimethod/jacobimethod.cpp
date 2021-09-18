/*
Last edited by: Erlend T. North 10:21 01.10.2019
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <fstream>
// Tried catch.cpp at first, but did not have time to learn

using namespace std;
namespace ch = std::chrono;

void offdiag(arma::mat A, int *p, int *q, int n) {
   double max;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;  *p = i; *q = j;
           }
       }
   }
}

void Jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l, int n ) {
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }

    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
  A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
  A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
  A(l,k) = 0.0;  // same here
  for ( int i = 0; i < n; i++ ) {
    if ( i != k && i != l ) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
  return;
}

bool test_offdiag() {
    // Checks for largest offdiag-element in matrices A, B, C and D
    int n = 5;
    arma::mat A = arma::mat(n,n, arma::fill::zeros);
    arma::mat B = arma::mat(n,n, arma::fill::zeros);
    arma::mat C = arma::mat(n,n, arma::fill::zeros);
    arma::mat D = arma::mat(n,n, arma::fill::zeros);
    int p, q;

    A(0,0) =  1.0; A(0,1) =  4.0; A(0,2) =    7.0; A(0,3) =   5.0; A(0,4) =   22.0;
    A(1,0) =  4.0; A(1,1) =  5.0; A(1,2) =    8.0; A(1,3) =  57.0; A(1,4) =   62.0;
    A(2,0) =  7.0; A(2,1) =  8.0; A(2,2) =    9.0; A(2,3) = -52.0; A(2,4) = 1034.0;
    A(3,0) =  5.0; A(3,1) = 57.0; A(3,2) =  -52.0; A(3,3) =  25.0; A(3,4) =   83.0;
    A(4,0) = 22.0; A(4,1) = 62.0; A(4,2) = 1034.0; A(4,3) =  83.0; A(4,4) =  -38.0;
    offdiag(A, &p, &q, n);
    double A_offdiag = A(p,q);
    bool A_pass = A_offdiag == A(2,4);
    cout << "   Test for A returned: " << A_pass << endl;

    B(0,0) =    1.0; B(0,1) =      22.0; B(0,2) = 25.0; B(0,3) = 0.028; B(0,4) =    0.0083;
    B(1,0) =   22.0; B(1,1) =       6.0; B(1,2) =  8.0; B(1,3) = -87.0; B(1,4) = 1000000.0;
    B(2,0) =   25.0; B(2,1) =       8.0; B(2,2) =  5.0; B(2,3) =  17.0; B(2,4) =      94.0;
    B(3,0) =  0.028; B(3,1) =     -87.0; B(3,2) = 17.0; B(3,3) =  31.0; B(4,3) =      87.0;
    B(4,0) = 0.0083; B(4,1) = 1000000.0; B(4,2) = 94.0; B(4,3) =  87.0; B(4,4) =      31.0;
    offdiag(B, &p, &q, n);
    double B_offdiag = B(p,q);
    bool B_pass = B_offdiag == B(1,4);
    cout << "   Test for B returned: " << B_pass << endl;

    C(0,0) =   2.0; C(0,1) =  7.0; C(0,2) =   4.0; C(0,3) = 37.56; C(0,4) = 91.0;
    C(1,0) =   7.0; C(1,1) = 12.0; C(1,2) =   6.0; C(1,3) =  89.0; C(1,4) = 75.0;
    C(2,0) =   4.0; C(2,1) =  6.0; C(2,2) =   3.0; C(2,3) = 25.24; C(2,4) = 82.0;
    C(3,0) = 37.56; C(3,1) = 89.0; C(3,2) = 25.24; C(3,3) =  23.0; C(3,4) = 13.0;
    C(4,0) =  91.0; C(4,1) = 75.0; C(4,2) =  82.0; C(4,3) =  13.0; C(4,4) = 23.0;
    offdiag(C, &p, &q, n);
    double C_offdiag = C(p,q);
    bool C_pass = C_offdiag == C(0,4);
    cout << "   Test for C returned: " << C_pass << endl;

    D(0,0) =   72.0; D(0,1) = 3.1415; D(0,2) =   2.718; D(0,3) =    53.0; D(0,4) = 64.27;
    D(1,0) = 3.1415; D(1,1) =  188.0; D(1,2) =  3.1416; D(1,3) =  -183.0; D(1,4) =  13.0;
    D(2,0) =  2.718; D(2,1) = 3.1416; D(2,2) =    68.0; D(2,3) = -121.57; D(2,4) =  -3.0;
    D(3,0) =   53.0; D(3,1) = -183.0; D(3,2) = -121.57; D(3,3) =   1.119; D(3,4) =  87.0;
    D(4,0) =  64.27; D(4,1) =   13.0; D(4,2) =    -3.0; D(4,3) =    87.0; D(4,4) = 1.119;
    offdiag(D, &p, &q, n);
    double D_offdiag = D(p,q);
    bool D_pass = D_offdiag == D(1,3);
    cout << "   Test for D returned: " << D_pass << endl;

    return (A_pass && B_pass && C_pass && D_pass);
}

bool test_eig() {
    // Checks whether our algorithm returns the correct eigenvalues for a given matrix A
    double tolerance =1.0E-10;
    int maxiteration = 100000;
    int n = 5;
    double maxnondiag = 1;
    int iterations = 0;



    arma::mat A = arma::mat(n, n, arma::fill::zeros);
    arma::mat R = arma::mat(n, n, arma::fill::eye);

    for(int i = 0; i < n; i++) {
          for(int j = 0; j < n; j++) {
            if(i == j) {
              A(i,j) = 2;
            }
            else if(i == j+1) {
                A(i,j) = -1;
            }
            else if(i == j-1) {
                A(i,j) = -1;
            }
        }
    }

    while ( maxnondiag > tolerance && iterations <= maxiteration) {
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
    }

    arma::vec eigval = arma::vec(n);
    arma::vec A_diag = arma::vec(n);



    for(int j = 0; j < n; j++) {
        eigval(j) = 2 + 2*(-1)*std::cos( ((j+1)*M_PI) / (5+1) ); // Analytical expression for eigenvalues of A
        A_diag(j) = (A(j,j));
    }

    A_diag = arma::sort(A_diag);

    return arma::approx_equal(A_diag, eigval, "absdiff", tolerance);
}

int main(int argc, char *argv[]) { // argv[1]: dimension, argv[2]: bool for running tests
    //std::cout << std::scientific;

    fstream outfile; // Initializing outfile for later use

    //Tests
    if(atoi(argv[2])) {
        cout << endl << "Commencing tests..." << endl;

        cout << endl << "Testing for largest non-diagonal element:" << endl;
        if(test_offdiag()) {
            cout << "Largest non-diagonal element preserved!" << endl;
        }
        else {
            cout << "Largest non-diagonal element NOT preserved!" << endl;
        }

        cout << endl << "Testing for eigenvalues:" << endl;
        if(test_eig()) {
            cout << "Eigenvalues correct!" << endl;
        }
        else {
            cout << "Eigenvalues NOT correct!" << endl;
        }

        cout << endl;
    }


    int n = atoi(argv[1]);
    // We "move" h^2 outside of the matrix and "scale it" away

    //
    // Tridiagonal Toeplitz matrix
    //

    arma::mat A(n, n, arma::fill::zeros);
    arma::mat R(n, n, arma::fill::eye);

    for(int i = 0; i < n; i++) { // Filling A
          for(int j = 0; j < n; j++) {
            if(i == j) {
              A(i,j) = 2;
            }
            else if(i == j+1) {
                A(i,j) = -1;
            }
            else if(i == j-1) {
                A(i,j) = -1;
            }
        }
    }


    arma::vec eigval;
    arma::mat eigvec;

    ch::steady_clock::time_point start = ch::steady_clock::now();

    arma::eig_sym(eigval, eigvec, A);

    ch::steady_clock::time_point stop = ch::steady_clock::now();
    ch::duration<double> time_span_eig_sym = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by eig_sym = " << time_span_eig_sym.count()  << "s" << std::endl;


    double tolerance = 1.0E-10;
    int maxiteration = n*n*n;
    double maxnondiag = 1;

    int iterations = 0;
/*
    start = ch::steady_clock::now();

    while ( maxnondiag > tolerance && iterations <= maxiteration) {
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
    }

    stop = ch::steady_clock::now();
    if(iterations > maxiteration) {
        cout << "Reached max iterations" << endl;
    }
    else {
        cout << "Algorithm complete with " << iterations << " iterations!" << endl;
    }
    ch::duration<double> time_span_ours = ch::duration_cast<ch::nanoseconds>(stop - start);
    std::cout << "Time used by us = " << time_span_ours.count()  << "s" << std::endl;
*/

    /*
    outfile.open("../../stats.txt", std::fstream::out | std::ofstream::app);
    outfile << std::scientific ;
    outfile << n << ", " << iterations << ", " << time_span_eig_sym.count() << ", " << time_span_ours.count() << endl;
    outfile.close();
    */

    //
    // Harmonic Oscillator
    //
    cout << endl << "Entering the quantum domain!" << endl;

    double* rho = new double[n];
    rho[n-1] = 10;
    double h = (rho[n-1] - rho[0])/(n+1); // Feels wrong with (n+1), but Morten's EigvalueArma.cpp did it like this and gave seemingly better results

    for(int i = 0; i < n; i++) { // Filling rho
        rho[i] = rho[0] + (i+1)*h;
    }

    A = arma::mat(n, n, arma::fill::zeros); // Resetting A
    R = arma::mat(n, n, arma::fill::zeros); // Resetting R
    maxnondiag = 1;
    iterations = 0;
    for(int i = 0; i < n; i++) { // Filling A (harmonic)
          for(int j = 0; j < n; j++) {
            if(i == j) {
              A(i,j) = 2.0/(h*h) + (rho[i]*rho[i]); // Could move constants outside of loop, but not *that* necessarry
            }
            else if(i == j+1) {
                A(i,j) = -1.0/(h*h);
            }
            else if(i == j-1) {
                A(i,j) = -1.0/(h*h);
            }
        }
    }

    /* Armadillo's solver for harmonic
    eig_sym(eigval, eigvec, A);
    for(int i = 0; i < 5; i++) {
        cout << eigval[i] << endl;
    }
    */

/*
    while ( maxnondiag > tolerance && iterations <= maxiteration) { // main-loop for diagonalizing Quantum harmonic oscillator A
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
    }
    if(iterations > maxiteration) {
        cout << "Reached max iterations." << endl;
    }
    else {
        cout << "Algorithm complete with " << iterations << " iterations!" << endl;
    }

    eigval = sort(A.diag());
    cout << "eigval: " << endl;
    for(int i = 0; i < 5; i++) { // Prints five first eigenvalues
        cout << eigval[i] << endl;
    }


    outfile.open("../../harmonic.txt", std::fstream::out | std::ofstream::app);
    outfile << n << ", " << rho[n-1] << ", " << eigval[0] << ", " << eigval[1] << ", " << eigval[2] << ", " << eigval[3] << ", " << eigval[4] << endl;
    outfile.close();
*/

    //
    // Two-electron Quantum Dot
    //
    cout << endl << "Look out for Quantum dots!" << endl;

    double omega_r = 5;

    A = arma::mat(n, n, arma::fill::zeros); // Resetting A
    R = arma::mat(n, n, arma::fill::zeros); // Resetting R
    maxnondiag = 1;
    iterations = 0;
    for(int i = 0; i < n; i++) { // Filling A (qDot)
        for(int j = 0; j < n; j++) {
            if(i == j) {
                A(i,j) = 2.0/(h*h) + omega_r*omega_r*rho[i]*rho[i] + 1.0/rho[i] ;
            }
            else if(i == j+1) {
                A(i,j) = -1.0/(h*h);
            }
            else if(i == j-1) {
                A(i,j) = -1.0/(h*h);
            }
        }
    }

    /* Armadillo's solver for qDot
    eig_sym(eigval, eigvec, A);
    for(int i = 0; i < 5; i++) {
        cout << eigval[i] << endl;
    }
    */

    while ( maxnondiag > tolerance && iterations <= maxiteration) { // main-loop for diagonalizing Quantum dot A
       int p, q;
       offdiag(A, &p, &q, n);
       maxnondiag = fabs(A(p,q));
       Jacobi_rotate(A, R, p, q, n);
       iterations++;
    }
    if(iterations > maxiteration) {
        cout << "Reached max iterations." << endl;
    }
    else {
        cout << "Algorithm complete with " << iterations << " iterations!" << endl;
    }

    eigval = sort(A.diag());
    cout << "eigval: " << endl;
    for(int i = 0; i < 5; i++) { // Prints the five first eigenvalues
        cout << eigval[i] << endl;
    }

    outfile.open("../../qdot.txt", std::fstream::out | std::ofstream::app);
    outfile << n << ", " << rho[n-1] << ", " << omega_r << ", " << eigval[0] << ", " << eigval[1] << ", " << eigval[2] << ", " << eigval[3] << ", " << eigval[4] << endl;
    outfile.close();


    return 0;
}
