#include <cmath>
#include "eos.hpp"


double soundspeed(double &rho, double &p, double &gamma){
    return sqrt(gamma * p / rho);
}


double pres(double &dens, double &ene, double &m1, double &m2, double &m3, double &gamma){
    return (ene - 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / dens) * (gamma - 1.);
}


double ene(double &rho, double &pres, double &v1, double &v2, double &v3, double &gamma){
    return pres / (gamma - 1.) + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3);
}
