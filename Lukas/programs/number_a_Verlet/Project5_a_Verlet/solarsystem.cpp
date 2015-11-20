#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <iomanip>
#include <cmath>

using namespace arma;
using namespace std;

solarsystem::solarsystem()
{

}

void solarsystem::addplanet(Planet Planet1){
    planets.push_back(Planet1); // I add the planet to the vector of planets
}

void solarsystem::VelocityVerlet(double dt,int n){
    setmatrices();
    cout << A;
    calculateForces(); //this might be shit
    A.col(1)=A.col(1)+0.5*dt*dA.col(1);
    for(int k=0;k<n;k++){
        A.col(0)=A.col(0)+A.col(1)*dt;      // r(t+dt)
        calculateForces();                  // a(t+dt)
        A.col(1)=A.col(1)+dt*dA.col(1); // v(t+0.5dt) -> v(t+3/2dt)
        cout<<endl<<A(3,0)<<"  "<<A(4,0)<<endl;
        for(int i=0;i<planets.size();i+=1){
            for(int j=0;j<3;j++){
                planets[i].position[j]=A(3*i+j,0);
                planets[i].velocity[j]=A(3*i+j,1);
            }
        }
    }

}
void solarsystem::setmatrices(){
    A=Mat<double>(planets.size()*3,2);
    dA=Mat<double>(planets.size()*3,2);

    for(int i=0;i<planets.size();i+=1){
        for(int j=0;j<3;j++){
            A(3*i+j,0)=planets[i].position[j];
            A(3*i+j,1)=planets[i].velocity[j];
        }
    }

    for(int i=0;i<planets.size();i+=1){
        for(int j=0;j<3;j++){
            dA(3*i+j,0)=planets[i].velocity[j];
            dA(3*i+j,1)=planets[i].force[j]/planets[i].m;
        }
    }
}


void solarsystem::calculateForces(){
    double G=4*M_PI*M_PI; // referred to the earths solar system

    for(int k=0;k<planets.size();k++){
        for(int l=0;l<3;l++){
            planets[k].force[l]=0;
        }
    }
    for(int i=0;i<planets.size();i++){
        for(int j=i+1;j<planets.size();j++){
            double dx=planets[i].position[0]-planets[j].position[0];
            double dy=planets[i].position[1]-planets[j].position[1];
            double dz=planets[i].position[2]-planets[j].position[2];
            double dr2=dx*dx+dy*dy+dz*dz;
            //calculate forces
            double Fx=(G*(planets[i].m)*(planets[j].m)*dx)/pow(dr2,1.5);
            double Fy=(G*(planets[i].m)*(planets[j].m)*dy)/pow(dr2,1.5);
            double Fz=(G*(planets[i].m)*(planets[j].m)*dz)/pow(dr2,1.5);
            //update planet properties (not necessary)
            planets[i].force[0]-=Fx;
            planets[j].force[0]+=Fx;
            planets[i].force[1]-=Fy;
            planets[j].force[1]+=Fy;
            planets[i].force[2]-=Fz;
            planets[j].force[2]+=Fz;
        }
    }

    for(int i=0;i<planets.size();i+=1){
        for(int j=0;j<3;j++){
            dA(3*i+j,0)=planets[i].velocity[j];
            dA(3*i+j,1)=planets[i].force[j]/planets[i].m;
        }
    }

}
