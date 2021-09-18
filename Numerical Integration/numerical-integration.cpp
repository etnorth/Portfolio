/*
  Last edited: 21.10.2019 20:50 by Alexandra Jahr Kolstad
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <random>
#include <mpi.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;
namespace ch = std::chrono;

// Functions
double* linspace(double start,double stop, int n) { // Creates linspaced dynamic array
    double h = (stop - start)/(n-1);
    double* arr = new double[n];
    for(int i = 0; i < n; i++) {
      arr[i] = start + stop*i*h;
    }
    return arr;
}

double psi(double x1, double y1, double z1, double x2, double y2, double z2, double alpha = 2) { // Defines the function to integrate
    double value = exp(-2*alpha*(sqrt(x1 * x1 + y1 * y1 + z1 * z1) + sqrt(x2 * x2 + y2 * y2 + z2 * z2)));
    double length = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
      if(length < ZERO) {
          return 0;
        }
    return value / length ;
}

double psi_sphere(double r1, double r2, double t1, double t2, double p1, double p2, double alpha = 2) { // Defines the function to integrate in spherical coordinates
    double cosb = cos(t1) * cos(t2) + sin(t1) * sin(t2) * cos(p1-p2);
    double value = exp(-3 * (r1+r2) )* r1 * r1 * r2 * r2 * sin(t1) * sin(t2);
    double length = r1*r1 + r2*r2 - 2 * r1 * r2 * cosb;
      if(length < ZERO) {
          return 0;
        }
      else {
          return (value) / sqrt(length) ;
  }
}

double psi_sphere_MC(double r1, double r2, double t1, double t2, double p1, double p2, double alpha = 2) { // Defines the function to integrate in spherical coordinates
    double cosb = cos(t1) * cos(t2) + sin(t1) * sin(t2) * cos(p1-p2);
    double value = r1 * r1 * r2 * r2 * sin(t1) * sin(t2);
    double length = r1*r1 + r2*r2 - 2 * r1 * r2 * cosb;
    if(length < ZERO) {
        return 0;
    }
    else {
        return (value) / sqrt(length) ;
    }
}

void gauleg(double x1, double x2, double x[], double w[], int n) { // Morten's code
    /*
    This  function takes the lower and upper limits of integration x1, x2
    calculates and returns the abcissas in x[0, ..., n-1] and the weights
    in w[0, ..., n-1] of length n of the Gauss-Legendre n-point quadrature formulae.
    */
    int         m,j,i;
    double      z1,z,xm,xl,pp,p3,p2,p1;
    double      const  pi = 3.14159265359;
    double      *x_low, *x_high, *w_low, *w_high;

    m  = (n + 1)/2; // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low  = x; // pointer initialization
    x_high = x + n - 1;
    w_low  = w;
    w_high = w + n - 1;

    for(i = 1; i <= m; i++) { // loops over desired roots
        z = cos(pi * (i - 0.25)/(n + 0.5));
        /*
        Starting with the above approximation to the ith root
        we enter the mani loop of refinement bt Newtons method.
        */

        do {
            p1 =1.0;
            p2 =0.0;

            // Loop up recurrence relation to get the Legendre polynomial evaluated at x

            for(j = 1; j <= n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
        }

            /*
    p1 is now the desired Legrendre polynomial. Next compute
    ppp its derivative by standard relation involving also p2,
    polynomial of one lower order.
    */

            pp = n * (z * p1 - p2)/(z * z - 1.0);
            z1 = z;
            z  = z1 - p1/pp;                   // Newton's method
        }

        while(fabs(z - z1) > ZERO);
        /*
        Scale the root to the desired interval and put in its symmetric
        counterpart. Compute the weight and its symmetric counterpart
        */

        *(x_low++)  = xm - xl * z;
        *(x_high--) = xm + xl * z;
        *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
        *(w_high--) = *(w_low++);
    }
}

double gammln( double xx){ //Morten's code
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void gaulag(double *x, double *w, int n, double alf){ //Morten's code
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln(n))/(pp*n*p2);
	}
}

void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &std){ //Morten's code with minor adjustments
    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double * x = new double [n];
    double x1, x2, y1, y2, z1, z2, f;
    double mc = 0.0;
    double sigma = 0.0;
    int i;
    double jacob = pow((b-a),6);

    for (i = 0; i < n; i++){
        x1 = distribution(generator)*(b-a)+a;
        x2 = distribution(generator)*(b-a)+a;
        y1 = distribution(generator)*(b-a)+a;
        y2 = distribution(generator)*(b-a)+a;
        z1 = distribution(generator)*(b-a)+a;
        z2 = distribution(generator)*(b-a)+a;
        f = psi(x1, y1, z1, x2, y2, z2);
        mc += f;
        x[i] = f;
    }
    mc = mc/n;
    for (i = 0; i < n; i++){
        sigma += (x[i] - mc)*(x[i] - mc);
    }
    sigma = sigma*jacob/n;
    std = sqrt(sigma)/sqrt(n);
    integral = mc*jacob;
    delete [] x;
}

void Polar_MonteCarlo_Importance(int n, double  &integral, double  &std){ //Morten's code with minor adjustments
    std::random_device rd;
    std::mt19937_64 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double * x = new double [n];
    double r1, r2, t1, t2, p1, p2, f,rr1,rr2;
    double mc = 0.0;
    double sigma = 0.0;
    double jacob = 4*pow(M_PI,4)/16;
    int i;

    for (i = 0; i < n; i++){
        rr1 = distribution(generator);
        rr2 = distribution(generator);
        r1 = -0.25*log(1-rr1);
        r2 = -0.25*log(1-rr2);
        t1 = distribution(generator)*M_PI;
        t2 = distribution(generator)*M_PI;
        p1 = distribution(generator)*2*M_PI;
        p2 = distribution(generator)*2*M_PI;
        f = psi_sphere_MC(r1, r2, t1, t2, p1, p2);
        mc += f;
        x[i] = f;
    }
    mc = mc/n;
    for (i = 0; i < n; i++){
        sigma += (x[i] - mc)*(x[i] - mc);
    }
    sigma = sigma*jacob/n;
    std = sqrt(sigma)/sqrt(n);
    integral = mc*jacob;
    delete [] x;
}

int main(int argc, char *argv[]) {

    int n = atoi(argv[1]);      //Command line arguments
    double la = atof(argv[2]);

    double exact = (5*M_PI*M_PI)/(16*16);

    //---------------------------------------------------------------------------------------------

    // Gauss-Legendre method

    double *w = new double[n];
    double *x = new double[n];

    //gauleg(-la, la, x, w, n);

    ch::steady_clock::time_point start = ch::steady_clock::now();

    double legendre_sum = 0.0;
    // Calculating the approximation of the integral with six dimensions for Gauss-Legendre
    /*
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                for(int l = 0; l < n; l++) {
                    for(int o = 0; o < n; o++) {
                        for(int p = 0; p < n; p++) {
                          legendre_sum += (w[i] * w[j] * w[k] * w[l] * w[o] * w[p]) * psi(x[i], x[j], x[k], x[l], x[o], x[p]);
                        }
                    }
                }
            }
        }
    }*/

    ch::steady_clock::time_point stop = ch::steady_clock::now();
    ch::duration<double> time_span_gauss_legendre = ch::duration_cast<ch::nanoseconds>(stop - start);

    delete [] x;
    delete [] w;



    //-------------------------------------------------------------------------------------------

    // Gauss-Laguerre method

    double *r = new double[n+1];
    double *the = new double[n];
    double *phi = new double[n];
    double *wr = new double[n+1];
    double *wthe = new double[n];
    double *wphi = new double[n];

    //gaulag(r, wr, n+1, 0);
    //gauleg(0, M_PI, the, wthe, n);
    //gauleg(0, 2*M_PI, phi, wphi, n);

    start = ch::steady_clock::now();


    double laguerre_sum = 0.0;
    // Calculating the approximation of the integral with six dimensions for Gauss-Laguerre
    /*
    for(int i = 1; i < n+1; i++) {
        for(int j = 1; j < n+1; j++) {
            for(int k = 0; k < n; k++) {
                for(int l = 0; l < n; l++) {
                    for(int o = 0; o < n; o++) {
                        for(int p = 0; p < n; p++) {
                          laguerre_sum += (wr[i] * wr[j] * wthe[k] * wthe[l] * wphi[o] * wphi[p]) * psi_sphere(r[i], r[j], the[k], the[l], phi[o], phi[p]);
                        }
                    }
                }
            }
        }
    }*/

    stop = ch::steady_clock::now();
    ch::duration<double> time_span_gauss_laguerre = ch::duration_cast<ch::nanoseconds>(stop - start);

    delete [] r;
    delete [] the;
    delete [] phi;
    delete [] wr;
    delete [] wthe;
    delete [] wphi;



    //-------------------------------------------------------------------------------------------

    //Brute Force Monte Carlo

    double BMC_sum;
    double BMC_std;

    start = ch::steady_clock::now();

    Brute_MonteCarlo(n, -la, la, BMC_sum, BMC_std);

    stop = ch::steady_clock::now();
    ch::duration<double> time_span_BMC = ch::duration_cast<ch::nanoseconds>(stop - start);

    //-------------------------------------------------------------------------------------------

    //Spherical Monte Carlo

    double SMC_sum;
    double SMC_std;

    start = ch::steady_clock::now();

    Polar_MonteCarlo_Importance(n, SMC_sum, SMC_std);

    stop = ch::steady_clock::now();
    ch::duration<double> time_span_SMC = ch::duration_cast<ch::nanoseconds>(stop - start);

    //-------------------------------------------------------------------------------------------

    //Parallelized Spherical Monte Carlo

    //  MPI initializations
    int local_n, numprocs, my_rank;
    double PSMC_sum, PSMC_std, local_sum, local_std;
    double time_start, time_end, time_span_PSMC;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    if(n%numprocs==0) {   // To test if n is an odd number or an even number
      // This is good
      local_n = n/numprocs;
    }
    else {
      // This is bad
      local_n = (n-1)/numprocs;
    }

    time_start = MPI_Wtime();
    PSMC_sum = 0.0;
    PSMC_std = 0.0;
    Polar_MonteCarlo_Importance(local_n, local_sum, local_std);
    MPI_Reduce(&local_sum, &PSMC_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_std, &PSMC_std, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    time_end = MPI_Wtime();
    PSMC_sum = PSMC_sum/numprocs;
    PSMC_std = PSMC_std/numprocs;
    time_span_PSMC = time_end-time_start;
    if ( my_rank == 0) {

        fstream outfile;

    /*
        outfile.open("lambda.txt", std::fstream::out | std::ofstream::app);
        outfile << n << " , " << la << " , " << fabs(exact-legendre_sum) << " , " << fabs(exact-laguerre_sum) << " , " << time_span_gauss_legendre.count() << " , " << time_span_gauss_laguerre.count() << endl;
        outfile.close();

        outfile.open("integrationpoints.txt", std::fstream::out | std::ofstream::app);
        outfile << n << " , " << la << " , " << fabs(exact-legendre_sum) << " , " << fabs(exact-laguerre_sum) << " , " << time_span_gauss_legendre.count() << " , " << time_span_gauss_laguerre.count() << endl;
        outfile.close();

        outfile.open("montecarlo.txt", std::fstream::out | std::ofstream::app);
        outfile << n << " , " << la << " , " << fabs(exact-BMC_sum) << " , " << fabs(exact-SMC_sum) << " , " << fabs(exact-PSMC_sum) << " , " << time_span_BMC.count() << " , "  << time_span_SMC.count() << " , "  << time_span_PSMC << endl;
        outfile.close();

        outfile.open("montecarlooptimized.txt", std::fstream::out | std::ofstream::app);
        outfile << n << " , " << la << " , " << fabs(exact-BMC_sum) << " , " << fabs(exact-SMC_sum) << " , " << fabs(exact-PSMC_sum) << " , " << time_span_BMC.count() << " , "  << time_span_SMC.count() << " , "  << time_span_PSMC << endl;
        outfile.close();

        outfile.open("variance.txt", std::fstream::out | std::ofstream::app);
        outfile << n << " , " << la << " , " << BMC_std*BMC_std*n << " , " << SMC_std*SMC_std*n << " , " << PSMC_std*PSMC_std*n << " , " << time_span_BMC.count() << " , "  << time_span_SMC.count() << " , "  << time_span_PSMC << endl;
        outfile.close();
    */

        // Final output
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << " " << "\n" ;
        cout << "Gauleg = " << setw(30) << setprecision(15)  << defaultfloat << legendre_sum << endl;
        cout << "Exact answer = " << setw(20) << setprecision(15) << defaultfloat << exact << endl;
        cout << "Error = " << setw(30) << setprecision(15) << defaultfloat << fabs(exact-legendre_sum) << endl;
        std::cout << "Time used by Gauleg = " << scientific << time_span_gauss_legendre.count()  << "s" << std::endl;
        cout << " " << "\n" ;
        cout << "GauLag = " << setw(30) << setprecision(15)  << defaultfloat << laguerre_sum << endl;
        cout << "Exact answer = " << setw(20) << setprecision(15) << defaultfloat << exact << endl;
        cout << "Error = " << setw(30) << setprecision(15) << defaultfloat << fabs(exact-laguerre_sum) << endl;
        std::cout << "Time used by Gaulag = " << scientific << time_span_gauss_laguerre.count()  << "s" << std::endl;
        cout << " " << "\n" ;
        cout << "BMC = " << setw(30) << setprecision(15) << defaultfloat << BMC_sum << endl;
        cout << "Exact answer = " << setw(20) << setprecision(15) << defaultfloat << exact << endl;
        cout << "Error = " << setw(30) << setprecision(15) << defaultfloat << fabs(exact-BMC_sum) << endl;
        std::cout << "Time used by BMC = " << scientific << time_span_BMC.count()  << "s" << std::endl;
        cout << " " << "\n" ;
        cout << "SMC = " << setw(30) << setprecision(15)  << defaultfloat << SMC_sum << endl;
        cout << "Exact answer = " << setw(20) << setprecision(15) << defaultfloat << exact << endl;
        cout << "Error = " << setw(30) << setprecision(15) << defaultfloat << fabs(exact-SMC_sum) << endl;
        std::cout << "Time used by SMC = " << scientific << time_span_SMC.count()  << "s" << std::endl;
        cout << " " << "\n" ;
        cout << "PSMC = " <<  setw(30) << setprecision(15) << defaultfloat << PSMC_sum << endl;
        cout << "Exact answer = " << setw(20) << setprecision(15) << defaultfloat << exact << endl;
        cout << "Error = " << setw(30) << setprecision(15) << defaultfloat << fabs(exact-PSMC_sum) << endl;
        cout << "Time used by PSMC  = " << scientific << time_span_PSMC << "s" << endl ;
        cout << "on number of processors: " << defaultfloat << numprocs << endl;
        cout << " " << "\n" ;
        cout << "Variance BMC = " << defaultfloat << BMC_std*BMC_std*n << "\n" ;
        cout << "Variance SMC = " << defaultfloat << SMC_std*SMC_std*n << "\n" ;
        cout << "Variance PSMC = " << defaultfloat << PSMC_std*PSMC_std*n << endl;
        cout << " " << "\n" ;
        cout << "Standard deviation BMC = " << defaultfloat << BMC_std << "\n" ;
        cout << "Standard deviation SMC = " << defaultfloat << SMC_std << "\n" ;
        cout << "Standard deviation PSMC = " << defaultfloat << PSMC_std << endl;
        cout << " " << "\n" ;

    }
    MPI_Finalize();


  return 0;
}
