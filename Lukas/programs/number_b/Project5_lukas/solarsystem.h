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
    void kinPotEnergy();
    void kinPotEnergy2();
    double* vecKinpot;
    double R0=20;
    double averagemass=10;
    int numplanets;
    double G;
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
    void VelocityVerlet(double dt,int n);
    void RungeKuttamethod(double dt,int n);
    mat derivate( Mat<double> B);
    void calculateForces();
    void addrandomplanet(double R_0);
};

#endif // SOLARSYSTEM_H
