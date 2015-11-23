#include "random.h"
#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <fstream>
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
    //print file
    //ofstream planetpositionVerlet;
    // planetpositionVerlet.open ("planetposition_verlet.txt");
    // planetpositionVerlet << "r_x          r_y         r_z         t ";
    setmatrices();
    calculateForces();                               //explonation of the algorithm
    A.col(1)=A.col(1)+0.5*dt*dA.col(1);              // v(t)->v(t+0.5dt)
    for(int k=0;k<n;k++){
        A.col(0)=A.col(0)+A.col(1)*dt;               // r(t+dt)
        calculateForces();                           // a(t+dt)
        A.col(1)=A.col(1)+dt*dA.col(1);              // v(t+0.5dt) -> v(t+3/2dt)
        //cout<<endl<<A(3,0)<<"  "<<A(4,0)<<"  "<<k<<endl;
        for(int i=0;i<planets.size();i+=1){
            //planetpositionVerlet<<endl;
            for(int j=0;j<3;j++){
                planets[i].position[j]=A(3*i+j,0);
                planets[i].velocity[j]=A(3*i+j,1);
                //planetpositionVerlet<<A(3*i+j,0)<<"  ";
            }
        }
    }
    //planetpositionVerlet.close();
}

// the function setmatrices initialises the matrix A and dA, which contains all necessary information about all of the planets
// dA=(d/dt)A -> this means that A contains: A=[positioncomponents of all planets;velocitycomponents of all planets]
// and dA contains: dA=[velocitycomponents of all planets;accelerationcomponents of all planets]
// This is very useful for RK4; However I also use it for Velocityverlet
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

// the function calculateForces recalculates the Forces of the actual constellation of planets
// the function calculateForces dates up the matrix dA as well as the planet objects
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
            //update planet properties
            planets[i].force[0]-=Fx;
            planets[j].force[0]+=Fx;
            planets[i].force[1]-=Fy;
            planets[j].force[1]+=Fy;
            planets[i].force[2]-=Fz;
            planets[j].force[2]+=Fz;
        }
    }
    //update the matrix dA
    for(int i=0;i<planets.size();i+=1){
        for(int j=0;j<3;j++){
            dA(3*i+j,0)=planets[i].velocity[j];
            dA(3*i+j,1)=planets[i].force[j]/planets[i].m;
        }
    }

}

//This will be useful by RK 4 algorithm
// the function derivate returns the derivation of B
mat solarsystem::derivate( Mat<double> B){
    double G=4*M_PI*M_PI; // referred to the earths solar system
    int n=planets.size();
    Mat<double> dC;
    dC=Mat<double>(n*3,2,fill::zeros);
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            double dx=B(3*i,0)-B(3*j,0);
            double dy=B(3*i+1,0)-B(3*j+1,0);
            double dz=B(3*i+2,0)-B(3*j+2,0);
            double dr2=dx*dx+dy*dy+dz*dz;
            //calculate forces
            double Fx=(G*(planets[i].m)*(planets[j].m)*dx)/pow(dr2,1.5);
            double Fy=(G*(planets[i].m)*(planets[j].m)*dy)/pow(dr2,1.5);
            double Fz=(G*(planets[i].m)*(planets[j].m)*dz)/pow(dr2,1.5);//until here B
            //update planet properties (not necessary)
            dC(3*i,1)-=Fx;
            dC(3*j,1)+=Fx;
            dC(3*i+1,1)-=Fy;
            dC(3*j+1,1)+=Fy;
            dC(3*i+2,1)-=Fz;
            dC(3*j+2,1)+=Fz;//returen dC with v from B and new acceleration
        }
    }
    for(int i=0;i<n;i+=1){
        for(int j=0;j<3;j++){
            dC(3*i+j,1)=dC(3*i+j,1)/planets[i].m;
            // cout << dC/planets[i].m<<endl;
        }
    }
    dC.col(0)=B.col(1);
    return(dC);
}

void solarsystem::RungeKuttamethod(double dt,int n){
    setmatrices();
    int size=planets.size();
    Mat<double> k1,k2,k3,k4;
    k1=Mat<double>(size*3,2);
    k2=Mat<double>(size*3,2);
    k3=Mat<double>(size*3,2);
    k4=Mat<double>(size*3,2);
    //ofstream RungeKutta_position;
    //RungeKutta_position.open("Rungekuttaposition.txt");
    for(int i=0;i<n;i++){
        k1=derivate(A);
        k2=derivate(A+k1*(dt/2));
        k3=derivate(A+k2*(dt/2));
        k4=derivate(A+k3*dt);
        A=A+(1.0/6.0)*(k1+2*k2+2*k3+k4)*dt;
        //RungeKutta_position << A(3,0)<<"  "<<A(4,0)<<"  "<<i<<endl;
    }
    //RungeKutta_position.close();
}

void addrandomplanet(double R_0){
    Planet randomplanet;
    // assume random_1 to be uniform distributed between 0 and 1
    double random_1,random_mass_gaussian;
    randomplanet.m=random_mass_gaussian;

    randomplanet.velocity[0]=0;
    randomplanet.velocity[1]=0;
    randomplanet.velocity[2]=0;

    randomplanet.position[0]=R_0*pow(random_1,(1.0/3.0))*pow((1-2*random_1)*(1-2*random_1),0.5)*cos(2*M_PI*random_1);
    randomplanet.position[0]=R_0*pow(random_1,(1.0/3.0))*pow((1-2*random_1)*(1-2*random_1),0.5)*sin(2*M_PI*random_1);
    randomplanet.position[0]=R_0*pow(random_1,(1.0/3.0))*(1-2*random_1);

}
