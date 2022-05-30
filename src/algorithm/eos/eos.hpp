#ifndef EOS_HPP_
#define EOS_HPP_
#include "../../defs.hpp"

#include <cmath>


class mesh;


double thermalspeed(double &rho, double &p, double &vthcoeff);


double soundspeed(double &rho, double &p, double &gamma);


double pres(double &dens, double &ene, double &m1, double &m2, double &m3, double &gamma);


double ene(double &rho, double &pres, double &v1, double &v2, double &v3, double &gamma);


void cons_to_prim(mesh &m);


#ifdef ENABLE_DUSTFLUID
void cons_to_prim_dust(mesh &m);
#endif // ENABLE_DUSTFLUID


#endif // EOS_HPP_
