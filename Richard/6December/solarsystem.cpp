#include "solarsystem.h"
#include "planet.h"
#include "gaussian_random.h"
#include <iostream>
#include <fstream>
#include <armadillo>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <random>


using namespace arma;
using namespace std;

solarsystem::solarsystem() : gen(this->rd()), dis(0,1)
{
    this->numplanets=0;
}


void solarsystem::addplanet(Planet Planet1){
    planets.push_back(Planet1); // I add the planet to the vector of planets
    average_kin.push_back(0);
    average_pot.push_back(0);
    this->numplanets+=1;
    G=4*M_PI*M_PI*R0*R0*R0/(32*this->numplanets*averagemass);
}

void solarsystem::VelocityVerlet(double dt,int n){
     //--------------------------------print file-----------------------------
    ofstream planetpositionVerlet;
    planetpositionVerlet.open ("planetposition_verlet.txt");
    ofstream kinPotFile;
    kinPotFile.open ("kinPotVerlet.txt");
    //----------------------------------------------------------------------
    setmatrices();
    calculateForces();                               //explonation of the algorithm
    A.col(1)=A.col(1)+0.5*dt*dA.col(1);              // v(t)->v(t+0.5dt)
    for(int k=0;k<n;k++){
        A.col(0)=A.col(0)+A.col(1)*dt;               // r(t+dt)
        calculateForces();                           // a(t+dt)
        A.col(1)=A.col(1)+dt*dA.col(1);              // v(t+0.5dt) -> v(t+3/2dt)
        for(int i=0;i<this->numplanets;i+=1){
            for(int j=0;j<3;j++){
                planets[i].position[j]=A(3*i+j,0);
                planets[i].velocity[j]=A(3*i+j,1);
                //outputfile-------------------------------------------------
                planetpositionVerlet<<A(3*i+j,0)<<"         ";
                //-----------------------------------------------------------
            }
        }
        //-------outpufile---------------------------------------------------
        planetpositionVerlet<<endl;  

            kinPotEnergy();
            kinPotFile<<numplanetsInSystem<<"       "<<theTotalEnergy;
            for (int i = 0; i < (this->numplanets); i++){
                kinPotFile<<setw(30)<<"         "<<potentialEnergy[i]<<"            "<<kineticEnergy[i];
            }
            kinPotFile<<endl;
    }
    kinPotFile.close();
    virial_output();
    planetpositionVerlet.close();
}

// the function setmatrices initialises the matrix A and dA, which contains all necessary information about all of the planets
// dA=(d/dt)A -> this means that A contains: A=[positioncomponents of all planets;velocitycomponents of all planets]
// and dA contains: dA=[velocitycomponents of all planets;accelerationcomponents of all planets]
// This is very useful for RK4; However I also use it for Velocityverlet
void solarsystem::setmatrices(){
    A=Mat<double>(this->numplanets*3,2);
    dA=Mat<double>(this->numplanets*3,2);
    
    for(int i=0;i<this->numplanets;i+=1){
        for(int j=0;j<3;j++){
            A(3*i+j,0)=planets[i].position[j];
            A(3*i+j,1)=planets[i].velocity[j];
        }
    }
    
    for(int i=0;i<this->numplanets;i+=1){
        for(int j=0;j<3;j++){
            dA(3*i+j,0)=planets[i].velocity[j];
            dA(3*i+j,1)=planets[i].force[j]/planets[i].m;
        }
    }
}

// the function calculateForces recalculates the Forces of the actual constellation of planets
// the function calculateForces dates up the matrix dA as well as the planet objects
void solarsystem::calculateForces(){

    
    for(int k=0;k<this->numplanets;k++){
        for(int l=0;l<3;l++){
            planets[k].force[l]=0;
        }
    }
    for(int i=0;i<this->numplanets;i++){
        for(int j=i+1;j<this->numplanets;j++){
            double dx=A(3*i,0)-A(3*j,0);
            double dy=A(3*i+1,0)-A(3*j+1,0);
            double dz=A(3*i+2,0)-A(3*j+2,0);
            double dr2=dx*dx+dy*dy+dz*dz+epsilon*epsilon;
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
    for(int i=0;i<this->numplanets;i+=1){
        for(int j=0;j<3;j++){
            dA(3*i+j,0)=planets[i].velocity[j];
            dA(3*i+j,1)=planets[i].force[j]/planets[i].m;
        }
    }
    
}

//This will be useful by RK 4 algorithm
// the function derivate returns the derivation of B
mat solarsystem::derivate( Mat<double> B){
    Mat<double> dC;
    dC=Mat<double>(this->numplanets*3,2,fill::zeros);
    for(int i=0;i<this->numplanets;i++){
        for(int j=i+1;j<this->numplanets;j++){
            double dx=B(3*i,0)-B(3*j,0);
            double dy=B(3*i+1,0)-B(3*j+1,0);
            double dz=B(3*i+2,0)-B(3*j+2,0);
            double dr2=dx*dx+dy*dy+dz*dz+epsilon*epsilon;
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
    for(int i=0;i<this->numplanets;i+=1){
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
    Mat<double> k1,k2,k3,k4;
    k1=Mat<double>(this->numplanets*3,2);
    k2=Mat<double>(this->numplanets*3,2);
    k3=Mat<double>(this->numplanets*3,2);
    k4=Mat<double>(this->numplanets*3,2);
    //ofstream RungeKutta_position;
    //RungeKutta_position.open("Rungekuttaposition.txt");
    //-----------------------------------------------------------------------
    ofstream kinPotFile2;
    kinPotFile2.open ("kinPotFileRK4.txt");
    kinPotFile2<<setprecision(10)<<setw(30) <<"totalE"<<setprecision(10)<<setw(30) <<"totalK"<<setprecision(10)<<setw(30) <<"totalp";
    for (int i = 0; i < this->numplanets; i++){
        kinPotFile2<< setprecision(10)<<setw(30)<<"k"<<i<<setprecision(10)<<setw(30) <<"p"<<i;
    }
    kinPotFile2<<endl;
    //-----------------------------------------------------------------------
    for(int i=0;i<n;i++){
        k1=derivate(A);
        k2=derivate(A+k1*(dt/2));
        k3=derivate(A+k2*(dt/2));
        k4=derivate(A+k3*dt);
        A=A+(1.0/6.0)*(k1+2*k2+2*k3+k4)*dt;
        //RungeKutta_position << A(3,0)<<"  "<<A(4,0)<<"  "<<i<<endl;
        //--------------------------------------------------------
            kinPotEnergy();
            kinPotFile2<<numplanetsInSystem<<"       "<<theTotalEnergy;
            for (int k = 0; k < this->numplanets; k++){
                kinPotFile2<<setw(30)<<"         "<<potentialEnergy[k]<<"            "<<kineticEnergy[k];
            }
            kinPotFile2<<endl;

        //------------------------------------------------------------
    }
    kinPotFile2.close();
    virial_output();
    //RungeKutta_position.close();
}

void solarsystem::addrandomplanet(){
    Planet randomplanet;
    randomplanet.m=object.generateGaussianNoise(10,1);
    
    randomplanet.velocity[0]=0;
    randomplanet.velocity[1]=0;
    randomplanet.velocity[2]=0;
    
    double rand_r=dis(gen);
    double rand_theta=dis(gen);
    double rand_phi=dis(gen);

    randomplanet.position[0]=R0*pow(rand_r,(1.0/3.0))*sqrt((1-(pow(1-2*rand_theta,2))))*cos(2*M_PI*rand_phi);
    randomplanet.position[1]=R0*pow(rand_r,(1.0/3.0))*sqrt((1-(pow(1-2*rand_theta,2))))*sin(2*M_PI*rand_phi);
    randomplanet.position[2]=R0*pow(rand_r,(1.0/3.0))*(1-(2*rand_theta));
    addplanet(randomplanet);
}

//here we calc the potenial energy and the kintic energy;
void solarsystem::kinPotEnergy(){

    kineticEnergy = new double[ this->numplanets];

    double kinetic = 0.0;
    double potenial = 0.0;
    double totalKinetic = 0.0;
    double totalPotenial = 0.0;
    double velocitySquared= 0.0;

    //kineticEnergy = 1/2 * m * v^2
    for (int i = 0; i < this->numplanets;i++){
        velocitySquared = A(3*i,1)*A(3*i,1)+A(3*i+1,1)*A(3*i+1,1)+A(3*i+2,1)*A(3*i+2,1);
        kineticEnergy[i] = 0.5*planets[i].m*velocitySquared;

        totalKinetic += kineticEnergy[i];
    }

    //potineal energy U = -G Mm/r----------------
    potentialEnergy = new double[ this->numplanets];
    double r;
    Mat<double> p = Mat<double>(this->numplanets,this->numplanets,fill::zeros); //matrix p will include the potenail form all planets.

    for (int i = 0;i < this->numplanets; i++ ){
        for (int j = i+1; j < this->numplanets; j++){
            r = (A(3*i,0)- A(3*j,0))*(A(3*i,0)- A(3*j,0)) + (A(3*i+1,0)- A(3*j+1,0))*(A(3*i+1,0)- A(3*j+1,0))+(A(3*i+2,0)- A(3*j+2,0))*(A(3*i+2,0)- A(3*j+2,0));
            r = sqrt(r);
            p(i,j)=-planets[i].m*planets[j].m*G/r;
            p(j,i)=p(i,j);
            totalPotenial += p(j,i);
        }
    }
    //calcuting the potenial energy for each planet from the marix P
    for(int i = 0; i < this->numplanets; i++){
        for(int j = 0; j < this->numplanets; j++){
            potenial += p(i,j);
        }
        potentialEnergy[i] = potenial;
    }
    theTotalEnergy = totalKinetic+totalPotenial;

    //Virial analysis
     numplanetsInSystem = 0;

    //clac the energy of the bound system ! virial !!!
    for (int i = 0;i < this->numplanets; i++ ){
           if (( kineticEnergy[i]+potentialEnergy[i])<0.0){
               numplanetsInSystem += 1;
               average_kin[i] += kineticEnergy[i];
               average_pot[i]+= potentialEnergy[i]*0.5;
           }

       }

}
//---------------------------------------------------------------

void solarsystem::virial_output(){
    ofstream virialanalysis;
    virialanalysis.open("Virialkoefficient.txt");
    for(int i=0;i<this->numplanets;i++){
       // virialanalysis<<i<<"   "<< vecKinpot[2+2*i]+vecKinpot[3+2*i]<<"    ";
       // virialanalysis<< average_pot[i]/average_kin[i]<<endl;
    }
    virialanalysis.close();
}
