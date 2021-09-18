/*
    Last edited: 20.11.2019 00:01 by Alexandra Jahr Kolstad
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>
#include <armadillo>
#include <iomanip>
#include <sstream>

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

void initialize(int L, double temp, arma::Mat<int> &lattice, double& E, double& M) {  // Initialize energy and magnetization
    // Setup lattice and initial magnetization
    double spin;
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            if(temp < 1.5) {
                lattice(x,y) = 1; // Spin orientation for the ground state
            }
            else {

            lattice(x,y) = -1;

            }
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

void Metropolis(int L, arma::Mat<int> &lattice, double& E, double& M, double *w, int &acc) { // The Metropolis algorithm
    acc = 0;

    // Loop over all spins
    for(int x = 0; x < L; x++) {
        for(int y = 0; y < L; y++) {
            //Find random position
            int ix = rand_frac(mt_gen)*L; // If problem put "(int)" in front of rand_frac
            int iy = rand_frac(mt_gen)*L;
            //int ix = rand_frac(mt_gen) * (L-1) ;
            //int iy = rand_frac(mt_gen) * (L-1) ;
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

                //acc += 1;

                // Update energy and magnetization
                M += 2*lattice(ix,iy);
                E += deltaE;

                acc += 1 ;
            }
        }
    }
}

void output(int L, int n, double temperature, double *average) { // Prints to file the results of the calculations
    double norm = 1.0/n; // Divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    //double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    // all expectation values are per spin, divide by 1/L/L

    double Evariance = (E2average - Eaverage*Eaverage)/L/L ;
    double Mabsvariance = (M2average - Mabsaverage*Mabsaverage)/L/L ;
    double Cv = Evariance/temperature/temperature ;
    double susc = Mabsvariance/temperature ;

    outfile << setprecision(8);
    outfile << temperature << " , " ;
    outfile << Eaverage/L/L << " , " ;
    outfile << Cv << " , " ;
    outfile << susc << " , " ;
    outfile << Mabsaverage/L/L << endl ;
}

int main(int argc, char *argv[]) { // Main function
    // Read in initial values
    double w[17], average[5], E, M, E2, M2;
    int L = atoi(argv[1]);
    int n = atoi(argv[2]);
    double initial_temp = atof(argv[3]);
    double final_temp = atof(argv[4]);
    double temp_step = atof(argv[5]);
    int acc = 0;

    ch::steady_clock::time_point start = ch::steady_clock::now();

    // Initialize lattice
    arma::Mat<int> lattice(L,L,arma::fill::zeros);

    //string outfilename = "L" + to_string(L) + "_n" + to_string(n) + "_Tvar";
    string outfilename = argv[6] ;

    outfile.open("../../output/" + outfilename + ".txt", fstream::app);
    outfile << "T , <E> , Cv , X , <|M|>" << endl;
    outfile.close();

    for(double temp = initial_temp; temp <= final_temp; temp += temp_step) {
        // Initialize energy and magnetization
        E = M = 0;

        // setup array for possible energy changes
        //for(int de = -8; de <= 8; de++) {
        for (int de = -8 ; de <= 8 ; de+=4){
            w[de+8] = 0;
        }
        //for(int de = -8; de <= 8; de+=4) {
        for (int de = -8 ; de <= 8 ; de+=4) {
            w[de+8] = exp(-de/temp);
        }

        // Initialize array for expectation values
        for(int i = 0; i < 5; i++) {
            average[i] = 0;
        }
        initialize(L, temp, lattice, E, M);

        // Output files

        /*
        stringstream tempstream;
        tempstream << fixed << setprecision(1) << temp;
        string tempstring = tempstream.str();
        */

        //string outfilename = "L" + to_string(L) + "_n" + to_string(n) + "_T" + tempstring + "_ord" + to_string(ordered);
        outfile.open("../../output/" + outfilename + ".txt", fstream::app);

        //outfile.open("../../output/" + outfilename + ".txt");
        // WAS HERE

        /*
        ofstream energyfile;
        energyfile.open("../../energy/" + outfilename + ".txt");
        energyfile << "E" << endl;

        ofstream magnetfile;
        magnetfile.open("../../magnet/" + outfilename + ".txt");
        magnetfile << "M" << endl;

        ofstream acceptfile;
        acceptfile.open("../../accept/" + outfilename + ".txt");
        acceptfile << "Acceptance" << endl;

        // Initial values to files
        energyfile << setprecision(8) << E/L/L << endl;
        magnetfile << setprecision(8) << M/L/L << endl;
        acceptfile << setprecision(8) << acc/L/L << endl;
        */

        // Start Monte Carlo Computation
        for(int cycles = 1; cycles <= n; cycles++) {
            Metropolis(L, lattice, E, M, w, acc);

            // Update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
            /*
            energyfile << E/L/L << endl;
            magnetfile << M/L/L << endl;
            acceptfile << acc/L/L << endl;
            */
        }

        // Print results
        output(L, n, temp, average);
        outfile.close();
        /*
        energyfile.close();
        magnetfile.close();
        acceptfile.close();
        */
        cout << (temp-initial_temp)/temp_step << " out of " << (final_temp-initial_temp)/temp_step << " iterations complete" << endl;
    }

    ch::steady_clock::time_point stop = ch::steady_clock::now();
    ch::duration<double> time_span = ch::duration_cast<ch::nanoseconds>(stop - start);
    cout << "Program used " << time_span.count() << " seconds" << endl;

    return 0;
}
