#include <iostream>
#include <fstream>
#include <cmath>
#include "planet.h"
#include <armadillo>

using namespace arma;
using namespace std;

double calculateforces(Planet p1,Planet p2,int i){
    double G=4*M_PI*M_PI;
    double dx=p1.position[0]-p2.position[0];
    double dy=p1.position[1]-p2.position[1];
    double dz=p1.position[2]-p2.position[2];
    double dr2=dx*dx+dy*dy+dz*dz;
        p1.force[i]=G*p1.m*p2.m*(p2.position[i]-p1.position[i])/pow(dr2,1.5);
        return(p1.force[i]);
}


double verletalgorithm(Planet p1,Planet p2,double dt, int n){
    ofstream verletposition;
    verletposition.open("verletposition.txt");
    for(int i=0;i<3;i++){
        p1.force[i]=calculateforces(p1,p2,i);
        p2.force[i]=calculateforces(p2,p1,i);
        p1.velocity[i]=p1.velocity[i]+0.5*p1.force[i]*dt/p1.m;
        p2.velocity[i]=p2.velocity[i]+0.5*p2.force[i]*dt/p2.m;
    }

        for(int k=0;k<n;k++){
            p1.position[0]=p1.position[0]+p1.velocity[0]*dt;
            p2.position[0]=p2.position[0]+p2.velocity[0]*dt;
            p1.position[1]=p1.position[1]+p1.velocity[1]*dt;
            p2.position[1]=p2.position[1]+p2.velocity[1]*dt;
            p1.position[2]=p1.position[2]+p1.velocity[2]*dt;
            p2.position[2]=p2.position[2]+p2.velocity[2]*dt;
            p1.force[0]=calculateforces(p1,p2,0);
            p2.force[0]=calculateforces(p2,p1,0);
            p1.force[1]=calculateforces(p1,p2,1);
            p2.force[1]=calculateforces(p2,p1,1);
            p1.force[2]=calculateforces(p1,p2,2);
            p2.force[2]=calculateforces(p2,p1,2);
            p1.velocity[0]=p1.velocity[0]+p1.force[0]*dt/p1.m;
            p2.velocity[0]=p2.velocity[0]+p2.force[0]*dt/p2.m;
            p1.velocity[1]=p1.velocity[1]+p1.force[1]*dt/p1.m;
            p2.velocity[1]=p2.velocity[1]+p2.force[1]*dt/p2.m;
            p1.velocity[2]=p1.velocity[2]+p1.force[2]*dt/p1.m;
            p2.velocity[2]=p2.velocity[2]+p2.force[2]*dt/p2.m;
            verletposition<<p2.position[0]<<"  "<<p2.position[1]<<"  "<<p1.position[0]<<"  "<<p1.position[1]<<endl;
        }
    verletposition.close();

}



int main()
{
    Planet p1(-1,0,0,0,-2,0,0,0,0,1) ;
    Planet p2(1,0,0,0,2,0,0,0,0,1);
    verletalgorithm(p1,p2,0.001,1000);

}
