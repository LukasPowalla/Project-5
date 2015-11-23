/*
 * This program includes the VerletVelocity-solver and Runge Kutta method (4th order) for a
 * two body system;
 *
 * For the Velocityverletalgorith, we use the position and velocity of the planet Objects.
 * functions needed for Velocity Verlet:
 * calculateforces -
 * verletalgorithm -
 *
 *
 * For RK4, we use a matrix representation; Planet-object positions and velocities are only needeed to initialize
 * the position and velocity sits in the constructed matrix A=(position, velocity).
 * Functions needed for RK4 method:
 * derivate -
 * rungeKuttaMethod -
 */


#include <iostream>
#include <fstream>
#include <cmath>
#include "planet.h"
#include <armadillo>

using namespace arma;
using namespace std;

// calculateforces calculates the forces on planet p1 and p2 and returns the force on the planet, which stands on the first place of the
// argument
double calculateforces(Planet p1,Planet p2,int i){
    double G=4*M_PI*M_PI;
    double dx=p1.position[0]-p2.position[0];
    double dy=p1.position[1]-p2.position[1];
    double dz=p1.position[2]-p2.position[2];
    double dr2=dx*dx+dy*dy+dz*dz;
    p1.force[i]=G*p1.m*p2.m*(p2.position[i]-p1.position[i])/pow(dr2,1.5);
    return(p1.force[i]);
}

// verletalgorithm performs the Velocity Verlet algorithm for endtime T= n* dt;
// only two body problem p1,p2
void verletalgorithm(Planet p1,Planet p2,double dt, int n){
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

// the function derivate returns a matrix, which is the derivative of the matrix B in the argument;
// the planet in the argument is only needed because we want to use the masses of the planets p1 and p2
mat derivate( Mat<double> B,Planet p1,Planet p2){
    double G=4*M_PI*M_PI; // referred to the earths solar system
    Mat<double> dC;
    dC=Mat<double>(6,2,fill::zeros);
    for(int i=0;i<2;i++){
        for(int j=i+1;j<2;j++){
            double dx=B(3*i,0)-B(3*j,0);
            double dy=B(3*i+1,0)-B(3*j+1,0);
            double dz=B(3*i+2,0)-B(3*j+2,0);
            double dr2=dx*dx+dy*dy+dz*dz;
            //calculate forces
            double Fx=(G*(p1.m)*(p2.m)*dx)/pow(dr2,1.5);
            double Fy=(G*(p1.m)*(p2.m)*dy)/pow(dr2,1.5);
            double Fz=(G*(p1.m)*(p2.m)*dz)/pow(dr2,1.5);//until here B
            //update planet properties
            dC(3*i,1)-=Fx;
            dC(3*j,1)+=Fx;
            dC(3*i+1,1)-=Fy;
            dC(3*j+1,1)+=Fy;
            dC(3*i+2,1)-=Fz;
            dC(3*j+2,1)+=Fz;//returen dC with v from B and new acceleration
        }
    }
    for(int j=0;j<3;j++){
        dC(j,1)=dC(j,1)/p1.m;
        dC(3+j,1)=dC(3+j,1)/p2.m;
    }
    dC.col(0)=B.col(1);
    return(dC);
}

// rungeKuttaMethod performs the Runge Kutta (4th order ) algorithm for a two body problem;
// we use the matrix representation A=(position;velocity)
void rungeKuttaMethod(Planet p1,Planet p2, double dt, int n){
    // Matrix initializing
    Mat <double>A;
    A= Mat <double>(6,2);

    for(int k=0;k<3;k++){
        A(k,0)=p1.position[k];
        A(k,1)=p1.velocity[k];
        A(3+k,0)=p2.position[k];
        A(3+k,1)=p2.velocity[k];
    }
    // end of Matrix initialising
    //beginn of the RK4 calculations; calculate the derivatives k1,k2,k3,k4
    Mat<double> k1,k2,k3,k4;
    k1=Mat<double>(6,2);
    k2=Mat<double>(6,2);
    k3=Mat<double>(6,2);
    k4=Mat<double>(6,2);
    ofstream RungeKutta_position;
    RungeKutta_position.open("Rungekuttaposition.txt");
    for(int i=0;i<n;i++){
        k1=derivate(A,p1,p2);
        k2=derivate(A+k1*(dt/2),p1,p2);
        k3=derivate(A+k2*(dt/2),p1,p2);
        k4=derivate(A+k3*dt,p1,p2);
        A=A+(1.0/6.0)*(k1+2*k2+2*k3+k4)*dt;
        RungeKutta_position <<A(3,0)<<"  "<<A(4,0)<<"  "<<A(0,0)<<"  "<<A(1,0)<<endl;
    }
    RungeKutta_position.close();
}
// end of RK4


int main()
{
    // Automatical data export included
    Planet p1(-1,0,0,0,-2,0,0,0,0,1) ;
    Planet p2(1,0,0,0,2,0,0,0,0,1);
    //verletalgorithm(p1,p2,0.001,1000);
    rungeKuttaMethod(p1,p2,0.001,1000);
}
