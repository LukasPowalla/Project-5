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
    int number_of_planets=0;
    vector<Planet> vectorofplanets;
    void addplanet(Planet Planet1);
    void VelocityVerlet(vector<Planet> Planetvec,double dT,double t_max);
    double ForceonPlanet(Planet stellar);
};

#endif // SOLARSYSTEM_H
