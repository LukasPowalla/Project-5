#ifndef PLANET_H
#define PLANET_H



class Planet
{
  public:
    Planet(double rx, double ry, double rz,double v_x, double v_y, double v_z, double mas);
    Planet();
    double position[3];
    double velocity[3];
    double force[3];
    double m;
};

#endif // PLANET_H
