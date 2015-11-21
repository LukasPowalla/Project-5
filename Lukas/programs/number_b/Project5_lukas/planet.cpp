#include "planet.h"

Planet::Planet()
{

}

Planet::Planet(double rx, double ry, double rz,double v_x, double v_y, double v_z,double f_x,double f_y, double f_z, double mas){
    position[0]=rx;
    position[1]=ry;
    position[2]=rz;
    velocity[0]=v_x;
    velocity[1]=v_y;
    velocity[2]=v_z;
    force[0]=f_x;
    force[1]=f_y;
    force[2]=f_z;
    m=mas;
}
