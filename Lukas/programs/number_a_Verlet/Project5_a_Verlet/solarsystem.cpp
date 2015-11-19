#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <iomanip>

using namespace arma;
using namespace std;

solarsystem::solarsystem()
{

}

void solarsystem::addplanet(Planet Planet1){
    number_of_planets+=1;               // the number of planets increases (+1)
    vectorofplanets.push_back(Planet1); // I add the planet to the vector of planets
}

void solarsystem::VelocityVerlet(vector<Planet> Planetvec,double dT,double t_max){


    for(int i=0;i<3;i++){
         Planetvec[0].position[i]+=Planetvec[0].velocity[i]*dT + ForceonPlanet(Planetvec[0])*0.5*dT*dT; //unsure with ForeceonmPlanet
    }

}

double solarsystem::ForceonPlanet(Planet stellar){


}
