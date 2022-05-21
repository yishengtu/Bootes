#ifndef EOS_HPP_
#define EOS_HPP_

#include <cmath>

float soundspeed(float &rho, float &p, float &gamma){
    return sqrt(gamma * p / rho);
}


float pres(float &dens, float &ene, float &m1, float &m2, float &m3, float &gamma){
    return (ene - 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / dens) * (gamma - 1.);
}


float ene(float &rho, float &pres, float &v1, float &v2, float &v3, float &gamma){
    return pres / (gamma - 1.) + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3);
}
#endif // EOS_HPP_
