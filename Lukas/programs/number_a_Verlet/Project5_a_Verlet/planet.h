#ifndef PLANET_H
#define PLANET_H



class Planet
{
  public:
    Planet(double x, double y, double z,double vx, double vy, double vz, double m);
    Planet();
    double position[3];
    double velocity[3];
    double m;
};

#endif // PLANET_H
