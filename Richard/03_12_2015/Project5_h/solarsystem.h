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
    double M_PI = 3.14159265359;
    void kinPotEnergy();
    void kinPotEnergy2();
    void radialDensity();//find the density vs radius
    double* centermass;
    double* vecKinpot;
    double* finalEnergy;//total Energy of i-th particle at the last timestep
    vector<double> average_kin; // <k>*N
    vector<double> average_pot; // <V>*N
    double R0=20;
    double averagemass=10;
    int numplanets;
    int numplanetsInSystem;
    double G;
    double epsilon=0.0005;
    solarsystem();
    vector<Planet> planets;
    Mat<double> A;
    Mat<double> dA;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    gaussian_random object;

    void addplanet(Planet Planet1);
    void setmatrices();
    void virial_output(); //outputs the virial-Koefficient and the total Energy of the i-th particle
    void VelocityVerlet(double dt,int n);
    void RungeKuttamethod(double dt,int n);
    mat derivate( Mat<double> B);
    void calculateForces();
    void addrandomplanet(double R_0);
};

#endif // SOLARSYSTEM_H
