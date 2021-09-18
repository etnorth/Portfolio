/*
    Last edited: 12.12.2019 14:05 by Erlend T. North
*/

#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

// "Mention" functions from other files for use here (body.h is now dependent on main.cpp)
double mod(arma::rowvec vec);

class Body {
    // Attributes
    string name;
    double mass;
    arma::mat pos;
    arma::rowvec vel;
    arma::rowvec acc;

public:
    // Methods
    Body(string new_name, int n) { // Constructs object
        name = new_name;
        pos = arma::mat(n,3,arma::fill::zeros);
        vel = arma::rowvec(3,arma::fill::zeros);
        acc = arma::rowvec(3,arma::fill::zeros);
        ifstream infile("../../BodyInput/" + name + ".txt");

        string trash;
        if(infile.is_open()) {
            infile >> name >> trash >> mass;
            infile >> pos(0,0) >> trash >> pos(0,1) >> trash >> pos(0,2);
            infile >> vel(0) >> trash >> vel(1) >> trash >> vel(2);
            infile.close();


            //vel = vel*365.242199; // Converts from AU/day to AU/yr

            /*
            ofstream data;
            data.open("../../bodyOutput/" + name + ".txt");
            data << "rx , ry , rz" << endl;
            data << setprecision(8) << pos(0,0) << " , " << pos(0,1) << " , " << pos(0,2) << endl;
            data.close();
            */
        }
        else {
            cout << "Unable to open file" << endl;
        }

    }

    string get_name() { // Fetches body.name
        return name;
    }

    double get_mass() { // Fetches body.mass
        return mass;
    }

    arma::mat get_pos() { // Fetches body.pos
        return pos;
    }

    arma::mat get_vel() { // Fetches body.vel
        return vel;
    }

    arma::mat get_acc() { // Fetches body.acc
        return acc;
    }

    void new_pos(arma::rowvec new_pos, int i) { // Assigns new position at iteration i
        pos.row(i) = new_pos;
    }

    void new_vel(arma::rowvec new_vel) { // Assigns new velocity
        vel = new_vel;
    }

    void new_acc(arma::rowvec new_acc) { // Assigns new acceleration
        acc = new_acc;
    }

    void acceleration(vector<Body> system, int i, double GMs=4*M_PI*M_PI) {
        arma::rowvec vec;
        acc = arma::rowvec(3,arma::fill::zeros);
        for(int bodyCount = 0; bodyCount < system.size(); bodyCount++) {
            if(system[bodyCount].get_name() != name) {
                vec = system[bodyCount].get_pos().row(i) - pos.row(i);
                acc += (GMs*system[bodyCount].get_mass()) / pow(mod(vec), 3)*vec;
            }
        }
    }

    void rel_acceleration(vector<Body> system, int i, double GMs=4*M_PI*M_PI, double c=63239.7263) {
        arma::rowvec vec, l;
        double r;
        acc = arma::rowvec(3,arma::fill::zeros);
        for(int bodyCount = 0; bodyCount < system.size(); bodyCount++) {
            if(system[bodyCount].get_name() != name) {
                vec = system[bodyCount].get_pos().row(i) - pos.row(i);
                r = mod(vec);
                l = mod(cross(pos.row(i),vec));
                acc += ((GMs*system[bodyCount].get_mass()) / pow(r, 3) * (1 + ((3*l*l)/(r*r*c*c))))*vec;
            }
        }
    }

    void write_data(int i) {
        ofstream data;
        data.open("../../bodyOutput/" + name + ".txt", fstream::app);
        data << pos(i,0) << " , " << pos(i,1) << " , " << pos(i,2) << endl;
        data.close();
    }

};
