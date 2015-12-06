#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include "planet.h"
#include "gaussian_random.h"
#include <armadillo>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using std::vector;
using namespace arma;

class solarsystem
{
public:
    // Constructor
    solarsystem();
    //constants and arrays used for Energy-analysis
    double* vecKinpot;
    double* finalEnergy;//total Energy of i-th particle at the last timestep
    vector<double> average_kin; // <k>*N
    vector<double> average_pot; // <V>*N
    //constants and vectors used for the method
    double R0=20;
    double averagemass=10;
    int numplanets;
    double G;
    double epsilon=0.002;
    vector<Planet> planets;
    Mat<double> A;
    Mat<double> dA;
    //random objects initialisation
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    gaussian_random object;
    //general functions
    void addplanet(Planet Planet1);
    void addrandomplanet();
    void setmatrices();
    //Velocity-Verlet
    void VelocityVerlet(double dt,int n);
    void calculateForces();
    // RK-4
    void RungeKuttamethod(double dt,int n);
    mat derivate( Mat<double> B);
    // funktions for Energy-analysis
    void kinPotEnergy();
    void kinPotEnergy2();
    void virial_output(); //outputs the virial-Koefficient and the total Energy of the i-th particle
};

#endif // SOLARSYSTEM_H
