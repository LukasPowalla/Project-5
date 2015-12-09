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
    this->numplanets+=1;
    G=4*M_PI*M_PI*R0*R0*R0/(32*this->numplanets*averagemass);
}

void solarsystem::VelocityVerlet(double dt,int n){
     //--------------------------------print file-----------------------------
    ofstream planetpositionVerlet;
    ofstream kinPotFile;
    ofstream ergodicFile;

    planetpositionVerlet.open ("planetposition_verlet.txt");
    kinPotFile.open ("kinPotVerlet.txt");
    ergodicFile.open ("ergodic.txt");
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

            }
        }
            //-------outpufile---------------------------------------------------
            kinPotEnergy();
            centermassfunction();
            //keeping the track for the cetnermass.
            planetpositionVerlet<<centermass[0]<<"         "<<centermass[1]<<"         "<<centermass[2]<<"         "<<endl;
            kinPotFile<<numplanetsInSystem<<"       "<<theTotalEnergy;
            for (int i = 0; i < (this->numplanets); i++){
                kinPotFile<<setw(30)<<"         "<<potentialEnergy[i]<<"            "<<kineticEnergy[i];
            }
            ergodicFile<<ergodic<<"         "<<k<<endl;
            kinPotFile<<endl;

    }
    kinPotFile.close();
    planetpositionVerlet.close();
    ergodicFile.close();
    radialDensity();//calc the radial density in the laste step where we have a equilibrium position ;
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
    //-----------------------------------------------------------------------
    ofstream RungeKutta_position;
    ofstream kinPotFile2;
    ofstream ergodicFile;
    RungeKutta_position.open("Rungekuttaposition.txt");
    kinPotFile2.open ("kinPotFileRK4.txt");
    ergodicFile.open ("ergodic.txt");
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
            kinPotEnergy();            
            centermassfunction();
            //keeping the track for the cetnermass.
            for(int z = 0; z < this->numplanets; z++){
            RungeKutta_position<< A(3*z,0)<<"         "<<A(3*z+1,0)<<"         "<<A(3*z+2,0)<<"         ";
            }
            RungeKutta_position<<centermass[0]<<"         "<<centermass[1]<<"         "<<centermass[2]<<"         "<<endl;
            kinPotFile2<<numplanetsInSystem<<"       "<<theTotalEnergy;
            for (int k = 0; k < this->numplanets; k++){
                kinPotFile2<<setw(30)<<"         "<<potentialEnergy[k]<<"            "<<kineticEnergy[k];
            }
            kinPotFile2<<endl;
            ergodicFile<<ergodic<<"         "<<i<<endl;

        //------------------------------------------------------------
    }
    kinPotFile2.close();
    RungeKutta_position.close();
    ergodicFile.close ();
}

void solarsystem::addrandomplanet(){
    Planet randomplanet;
    randomplanet.m=object.generateGaussianNoise(100,1);
    
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
    average_kin = 0.0;
    average_pot = 0.0;

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
        potentialEnergy[i] = potenial/2;
        potenial = 0;
    }
    theTotalEnergy = totalKinetic+totalPotenial;

    //Virial analysis
     numplanetsInSystem = 0;

    //clac the energy of the bound system ! virial !!!
    for (int i = 0;i < this->numplanets; i++ ){
        if (( kineticEnergy[i]+potentialEnergy[i])<0.0){
            numplanetsInSystem += 1;
            average_kin += kineticEnergy[i];
            average_pot += potentialEnergy[i];
        }

    }

    ergodic = average_pot/average_kin;
}

void solarsystem::centermassfunction(){
    centermass = new double[3];//(centrMassX,centerMassY,centerMassZ)
    double centermassX = 0.0;
    double centermassY = 0.0;
    double centermassZ = 0.0;
    double totalmass = 0.0;
    for (int i = 0;i<this->numplanets;i++){

     if (( kineticEnergy[i]+potentialEnergy[i])<0.0){
        totalmass  += planets[i].m;
        centermassX += A(3*i,0)*planets[i].m;
        centermassY += A(3*i+1,0)*planets[i].m;
        centermassZ += A(3*i+2,0)*planets[i].m;
      }
  }
    centermassX = centermassX/totalmass;
    centermassY = centermassY/totalmass;
    centermassZ = centermassZ/totalmass;

    centermass[0] = centermassX;
    centermass[1] = centermassY;
    centermass[2] = centermassZ;
}
void solarsystem::radialDensity(){

     double* radialDistance = new double[this->numplanets];
     double x = 0;
     double y = 0;
     double z = 0;
     double averageDistance = 0.0;
     double Nplanets = 0;
     double standardDeviation = 0.0;

     //looping though all bonded planets and clac the radial distance
     for (int i = 0;i<this->numplanets;i++){
        if (kineticEnergy[i]+potentialEnergy[i]<0){

            x = centermass[0] - A(3*i,0);
            y = centermass[1] - A(3*i+1,0);
            z = centermass[2] - A(3*i+2,0);
            double k = x*x+y*y+z*z;

            radialDistance[i] = sqrt(k);
            averageDistance += radialDistance[i];
            Nplanets += 1.0;

        }else{
            radialDistance[i] = 200; // just a big number to avoid worng calc(this plant is not in the system anyway).
        }
    }

    averageDistance = averageDistance/Nplanets;

     //standard deviation
     for (int i = 0;i<this->numplanets;i++){
        if (radialDistance[i] < 200){
            standardDeviation += (radialDistance[i] - averageDistance)*(radialDistance[i] - averageDistance);
        }
    }
     standardDeviation = sqrt(standardDeviation/Nplanets);

     cout<<"number of planet in the system : "<<Nplanets<<endl;
     cout<<"average distance" <<averageDistance<<endl;
     cout<<"standard deviation" <<standardDeviation<<endl;
     double m= 0.0;
     for (int i=0; i < this->numplanets;i++){
     m += planets[i].m;
    }
     cout<<m<<endl;

   double maxRadius = 100;
   double step = 0.001;
   double N = 0.0; // number of planet inside a givin radius.
   double volume = 0.0;
   int steps =  ((int)maxRadius/step + 1);
   double* density = new double[steps];
   int loopintger = 0;// just an intger for the loop;


     for (double r = 0.0; r < maxRadius; r += step){
       for (int i = 0;i < numplanetsInSystem; i++){
           if(radialDistance[i]<r){
               N +=1;
           }
         }
       volume = 4/3 * r*r*r*M_PI;
       density[loopintger]= N/volume;
       N=0;
       loopintger++;
     }

    //write data to a folder
    ofstream radialDensityTxt;
    radialDensityTxt.open ("radialDensity.txt");
    loopintger=0;
    for (double r = 0.0; r < maxRadius; r += step){
       radialDensityTxt<<r<<"         "<<density[loopintger]<<endl;
       loopintger++;
    }

     radialDensityTxt.close();
}
