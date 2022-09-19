#ifndef EOS_HPP_
#define EOS_HPP_
#include "../../defs.hpp"

#include <cmath>


class mesh;


double thermalspeed(double &rho, double &p, double &vthcoeff);

#pragma acc routine worker
double soundspeed(double &rho, double &p, double &gamma);

#pragma acc routine worker
double pres(double &dens, double &ene, double &m1, double &m2, double &m3, double &gamma);


double ene(double &rho, double &pres, double &v1, double &v2, double &v3, double &gamma);


void cons_to_prim(mesh &m);

#pragma acc routine worker
double energy_from_temperature_protection(double &dens, double &ene, double &m1, double &m2, double &m3, double &minTemp, double &gamma);


#endif // EOS_HPP_
