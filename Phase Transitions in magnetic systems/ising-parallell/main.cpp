/*
    Last edited: 17.11.2019 19:41 by Erlend T. North
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>
#include <armadillo>
#include <iomanip>
#include <sstream>
#include <mpi.h>

// kB = 1, J = 1

using namespace std;
namespace ch = std::chrono;

ofstream outfile;

long long seed = ch::high_resolution_clock::now().time_since_epoch().count(); //random_device did not work (always made the same sequence)
mt19937_64 mt_gen(seed);
uniform_real_distribution<double> rand_frac(0,1);

int periodic(int i, int limit, int add) { // Periodic Boundary Conditions
    return (i+limit+add) % (limit);
}

void initialize(int L, double temp, arma::Mat<int> &lattice, double& E, double& M, bool ordered) {  // Initialize energy and magnetization
    // Setup lattice and initial magnetization
    double spin;
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            if(ordered) {
                lattice(x,y) = 1;
            }
            else {
                spin = rand_frac(mt_gen);
                if(spin<0.5) {
                    lattice(x,y) = -1;
                }
                else {
                    lattice(x,y) = 1;
                }
            }

            /*
            if(temp < 1.5) {
                lattice(x,y) = 1; // Spin orientation for the ground state
            }
            */

            M += lattice(x,y);
        }
    }

    // Setup initial energy
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            E -= lattice(x,y) * (
                        lattice(x, periodic(y,L,-1)) +
                        lattice(periodic(x,L,-1), y)
                        );
        }
    }
}

void Metropolis(int L, arma::Mat<int> &lattice, double& E, double& M, double *w, int &acc, uniform_real_distribution<double> &rand_frac) { // The Metropolis algorithm
    acc = 0;

    // Loop over all spins
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            //Find random position
            int ix = rand_frac(mt_gen)*L; // If problem put "(int)" in front of rand_frac
            int iy = rand_frac(mt_gen)*L;
            int deltaE = 2*lattice(ix,iy) * (
                        lattice(ix,periodic(iy,L,-1)) +
                        lattice(periodic(ix,L,-1),iy) +
                        lattice(ix,periodic(iy,L, 1)) +
                        lattice(periodic(ix,L, 1),iy)
                        );

            // Here we perform the Metropolis test

            double r = rand_frac(mt_gen);

            double E_ = w[deltaE+8];

            if(r <= E_) {
                lattice(ix,iy) *= -1; // Flip one spin and accept new spin config

                acc += 1;

                // Update energy and magnetization
                M += 2*lattice(ix,iy);
                E += deltaE;
            }
        }
    }
}

void output(int L, int n, double temperature, double *average) { // Prints to file the results of the calculations
    double norm = 1.0/n; // Divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    // all expectation values are per spin, divide by 1/L/L
    double Evariance = (E2average- Eaverage*Eaverage)/L/L;
    double Mvariance = (M2average - Maverage*Maverage)/L/L;
    double Mabsvariance = (M2average - Mabsaverage*Mabsaverage)/L/L;
    outfile << setprecision(8);
    outfile << temperature << " , ";
    outfile << Eaverage/L/L << " , ";
    outfile << Evariance/temperature/temperature << " , ";
    // outfile << setw(15) << setprecision(8) << Maverage/L/L;
    outfile << Mabsvariance/temperature << " , ";
    outfile << Mabsaverage/L/L << endl;
}

int main(int argc, char *argv[]) { // Main function
    // Read in initial values
    double w[17], average[5], total_average[5], E, M, E2, M2;
    int L;
    int n;
    double initial_temp;
    double final_temp;
    double temp_step;
    bool ordered;
    int acc;
    string outfilename;
    long long seed;

    // MPI initializations

    int numprocs, my_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(my_rank == 0 && argc <= 1) {
        cout << "Bad usage: " << argv[0] << " read output file" << endl;
    }
    if(my_rank == 0 && argc > 1) {
        outfilename = argv[0];

        outfile.open("../../output/" + outfilename + ".txt");
        outfile << "T , <E> , Cv , X , <|M|>" << endl;
    }

    L = 20;
    n = 1000000;
    initial_temp = 1.0;
    final_temp = 2.4;
    temp_step = 0.1;

    /*
    Determine number of interval which are used by all processes
    myloop_begin gives the staring point on process my_rank
    myloop_end gives the end point of summation on process my_rank
    */
    int no_intervals = n/numprocs;
    int myloop_begin = my_rank*no_intervals + 1;
    int myloop_end = (my_rank+1)*no_intervals;
    if((my_rank == numprocs-1) && (myloop_end < n)) {
        myloop_end = n;
    }

    // Broadcast to all nodes common variables
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Initialize lattice
    arma::Mat<int> lattice(L,L,arma::fill::zeros);

    // Every node has its own seed for the random numbers, this is
    // important else if one starts with the same seed, one ends
    // with the same random numbers
    seed = ch::high_resolution_clock::now().time_since_epoch().count() + my_rank;
    mt19937_64 mt_gen(seed);
    uniform_real_distribution<double> rand_frac(0,1);



    for(double temp = initial_temp; temp <= final_temp; temp += temp_step) {
        // Initialize energy and magnetization
        E = M = 0;

        // Initialize lattice and expectation values
        initialize(L, temp, lattice, E, M, ordered);

        // setup array for possible energy changes
        for(int de = -8; de <= 8; de++) {
            w[de+8] = 0;
        }
        for(int de = -8; de <= 8; de+=4) {
            w[de+8] = exp(-de/temp);
        }

        // Initialize array for expectation values
        for(int i = 0; i < 5; i++) {
            average[i] = 0;
        }
        for(int i = 0; i < 5; i++) {
          total_average[i] = 0;
        }

        // Start Monte Carlo Computation
        for(int cycles = 1; cycles <= n; cycles++) {
            Metropolis(L, lattice, E, M, w, acc);

            // Update expectation values for local node
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }
        for(int i = 0; i < 5; i++) {
          MPI_Reduce(&average[i], &total_average[i], i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // Print results
        if(my_rank == 0){
          output(L, n, temp, average);
          cout << (temp-initial_temp)/temp_step << " out of " << (final_temp-initial_temp)/temp_step << " iterations complete" << endl;
        }
    }
    outfile.close();
    // End MPI
    MPI_Finalize();
    return 0;
}
