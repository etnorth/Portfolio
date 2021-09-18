/*
    Last edited: 12.12.2019 14:05 by Erlend T. North
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <armadillo>
#include <cmath>
#include <body.h>
#include <vector>
#include <chrono>

using namespace std;
namespace ch = chrono;

// Functions
double mod(arma::rowvec vec) { // Calculates the modulus of a vector
    double val = 0;
    for(int i = 0; i < vec.n_elem; i++) {
        val += vec(i)*vec(i);
    }
    return sqrt(val);
}

void forwardEuler(Body &body, vector<Body> system, int i, double dt) { // Performs 1 iteration of Forward Euler
    body.acceleration(system, i);
    body.new_vel(body.get_vel() + body.get_acc()*dt);
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt, i+1); // This is actually Euler-Cromer ;)
    //body.write_data(i);

}

void velocityVerlet(Body &body, vector<Body> system, int i, double dt, double dt_pos, double dt_vel) {
    body.acceleration(system, i); // Cannot re-use new_acc as other bodies have moved
    arma::rowvec acc_old = body.get_acc();
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt + acc_old*dt_pos, i+1);
    body.acceleration(system, i+1);
    body.new_vel(body.get_vel() + (acc_old+body.get_acc())*dt_vel);
}

void rel_forwardEuler(Body &body, vector<Body> system, int i, double dt) { // Performs 1 iteration of Forward Euler
    body.rel_acceleration(system, i);
    body.new_vel(body.get_vel() + body.get_acc()*dt);
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt, i+1); // This is actually Euler-Cromer ;)
    //body.write_data(i);

}

void rel_velocityVerlet(Body &body, vector<Body> system, int i, double dt, double dt_pos, double dt_vel) {
    body.rel_acceleration(system, i); // Cannot re-use new_acc as other bodies have moved
    arma::rowvec acc_old = body.get_acc();
    body.new_pos(body.get_pos().row(i) + body.get_vel()*dt + acc_old*dt_pos, i+1);
    body.acceleration(system, i+1);
    body.new_vel(body.get_vel() + (acc_old+body.get_acc())*dt_vel);
}

int main() {
    ifstream infile("../../input.txt");

    int n, yr;
    double dt, dt_pos, dt_vel;

    string trash, planet;
    vector<string> planetList;
    if(infile.is_open()) { // Read input-file
        infile >> trash >> trash >> trash; // "n" >> " , " >> "yr"
        infile >> n >> trash >> yr; // 100 >> " , " >> 100
        infile >> trash; // "planets"
        while(getline(infile, planet))  {
            planetList.push_back(planet);
        }
        planetList.erase(planetList.begin()); // 1st line is white-space

    }
    else { // File not readable
        cout << "Unable to open file" << endl;
        cout << "Quitting program..." << endl;
        return 1;
    }

    dt = (1.0*(double)yr)/(double)n;
    dt_pos = (dt*dt)/2.0;
    dt_vel = dt/2.0;

    vector<Body> system;
    vector<ofstream> systemdata;

    ch::steady_clock::time_point start;
    ch::steady_clock::time_point stop;

    //---------------------------------------------------------------------------------------------

    // Forward Euler
    for(int bodyCount = 0; bodyCount < planetList.size(); bodyCount++) { // Construct bodies and output-files
        system.push_back(Body(planetList[bodyCount], n));
        systemdata.push_back(ofstream("../../../data/ForwardEuler/" + system[bodyCount].get_name() + ".txt"));
        systemdata[bodyCount] << "x , y , z" << endl;
        systemdata[bodyCount] << system[bodyCount].get_pos()(0,0) << " , " << system[bodyCount].get_pos()(0,1) << " , " << system[bodyCount].get_pos()(0,2) << endl;
    }

    start = ch::steady_clock::now();

    // Main loop
    for(int i = 0; i < n-1; i++) { // Perform Forward Euler
        for(int bodyCount = 0; bodyCount < system.size(); bodyCount++) {
            if(system[bodyCount].get_name() == "Sun") { // Static Sun
            }
            else {
                forwardEuler(system[bodyCount], system, i, dt);
                //rel_forwardEuler(system[bodyCount], system, i, dt); // Sun-Mercury
            }
            systemdata[bodyCount] << system[bodyCount].get_pos()(i+1,0) << " , " << system[bodyCount].get_pos()(i+1,1) << " , " << system[bodyCount].get_pos()(i+1,2) << endl;
        }
    }

    stop = ch::steady_clock::now();
    ch::duration<double> time_span_forwardEuler = ch::duration_cast<ch::nanoseconds>(stop - start);
    cout << "Forward Euler complete!" << endl;


    for(int bodyCount = 0; bodyCount < systemdata.size(); bodyCount++) { // Close output-files
        systemdata[bodyCount].close();
    }
    systemdata.clear();
    system.clear();

    cout << "Time used by Forward Euler = " << time_span_forwardEuler.count()  << "s" << std::endl;

    //---------------------------------------------------------------------------------------------

    // Velocity Verlet
    for(int bodyCount = 0; bodyCount < planetList.size(); bodyCount++) { // Construct output-files
        system.push_back(Body(planetList[bodyCount], n));
        systemdata.push_back(ofstream("../../../data/VelocityVerlet/" + system[bodyCount].get_name() + ".txt"));
        systemdata[bodyCount] << "x , y , z" << endl;
        systemdata[bodyCount] << system[bodyCount].get_pos()(0,0) << " , " << system[bodyCount].get_pos()(0,1) << " , " << system[bodyCount].get_pos()(0,2) << endl;
    }

    start = ch::steady_clock::now();

    // Main loop
    for(int i = 0; i < n-1; i++) { // Perform Velocity Verlet
        for(int bodyCount = 0; bodyCount < system.size(); bodyCount++) {
            if(system[bodyCount].get_name() == "Sun") { // Static Sun
            }
            else {
                velocityVerlet(system[bodyCount], system, i, dt, dt_pos, dt_vel);
                //rel_velocityVerlet(system[bodyCount], system, i, dt, dt_pos, dt_vel); // Sun-Mercury
            }
            systemdata[bodyCount] << system[bodyCount].get_pos()(i+1,0) << " , " << system[bodyCount].get_pos()(i+1,1) << " , " << system[bodyCount].get_pos()(i+1,2) << endl;
        }
    }

    stop = ch::steady_clock::now();
    ch::duration<double> time_span_velocityVerlet = ch::duration_cast<ch::nanoseconds>(stop - start);
    cout << "Velocity Verlet complete!" << endl;

    for(int bodyCount = 0; bodyCount < systemdata.size(); bodyCount++) { // Close output-files
        systemdata[bodyCount].close();
    }
    systemdata.clear();
    system.clear();

    cout << "Time used by Velocity Verlet = " << time_span_velocityVerlet.count()  << "s" << std::endl;

    //system("cd ..");
    //system("cd ..");
    //system("python ../../plot.py");

    return 0;
}
