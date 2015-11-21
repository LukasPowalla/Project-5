/*
 * This will be a stupid straight forward Verlet solver for two bodies
 *
 * */

#include <iostream>
#include <fstream>
#include <cmath>
#include "planet.h"
#include "solarsystem.h"
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    solarsystem mysystem;
    Planet p1(0,0,0,0,0,0,0,0,0,1) ;
    Planet p2(1,0,0,0,6,0,0,0,0,0.000001);
    mysystem.addplanet(p1);
    mysystem.addplanet(p2);
    mysystem.VelocityVerlet(0.001,1000);
    //mysystem.RungeKuttamethod(0.001,1000);
}




//two bodys at the moment -old coldefragment - probabely used later
//this is the core of the Verlet algorithm
// alllda ... schreib ne FKT für die nächsten 3 zeilen; argumente sind ein ausgewählter planet und ein planet, der auf ihn wirkt
/*double forcevector_i[3];
for(int i=0;i<3;i++){
     forcevector_i[i]=(Planetvec[1].velocity[i]-Planetvec[0].velocity[i])*ForceonPlanet(Planetvec[0],Planetvec[1]);
     Planetvec[0].position[i]+=Planetvec[0].velocity[i]*dT +forcevector_i[i]*0.5*dT*dT; //unsure with ForeceonmPlanet
     Planetvec[0].velocity[i]+=0.5*(forcevector_i[i]+(Planetvec[1].velocity[i]-Planetvec[0].velocity[i])*ForceonPlanet(Planetvec[0],Planetvec[1]))*dT;
}*/
