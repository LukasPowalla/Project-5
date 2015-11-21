#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include "planet.h"
#include <armadillo>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <vector>

using std::vector;
using namespace arma;

class solarsystem
{
public:
    solarsystem();
    vector<Planet> planets;
    Mat<double> A;
    Mat<double> dA;

    void addplanet(Planet Planet1);
    void setmatrices();
    void VelocityVerlet(double dt,int n);
    void RungeKuttamethod(double dt,int n);
    mat derivate( Mat<double> B);
    void calculateForces();
};

#endif // SOLARSYSTEM_H
