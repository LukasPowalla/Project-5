/*
 *
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
    int numberofplanets=100;
    for(int i=0;i<numberofplanets;i++){
        mysystem.addrandomplanet(20);
    }

    /*Planet p1(0,0,0,0,0,0,0,0,0,1) ;
    Planet p2(1,0,0,0,6,0,0,0,0,0.000001);
    mysystem.addplanet(p1);
    mysystem.addplanet(p2);*/
    mysystem.VelocityVerlet(0.001,3000);
    //mysystem.RungeKuttamethod(0.001,10000);

    /*for(int k=0;k<numberofplanets;k++){
        cout<<mysystem.planets[k].position[0]<<endl;
        cout<<mysystem.planets[k].position[1]<<endl;
        cout<<mysystem.planets[k].position[2]<<endl;
        cout<<mysystem.planets[k].m<<endl<<endl;
    }*/
}



