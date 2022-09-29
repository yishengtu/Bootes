#include <cmath>
#include "eos.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"


double thermalspeed(double &rho, double &p, double &vthcoeff){
    return sqrt(vthcoeff * p / rho);
}


double soundspeed(double &rho, double &p, double &gamma){
    return sqrt(gamma * p / rho);
}


double pres(double &dens, double &ene, double &m1, double &m2, double &m3, double &gamma){
    return (ene - 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / dens) * (gamma - 1.);
}


double ene(double &rho, double &pres, double &v1, double &v2, double &v3, double &gamma){
    return pres / (gamma - 1.) + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3);
}


void cons_to_prim(mesh &m){
    #ifdef GPU
    #pragma acc parallel loop collapse (3)
    #else
    #pragma omp parallel for collapse (3) schedule (static)
    #endif
    for (int kk = m.x3s; kk < m.x3l ; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double v1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                double v2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                double v3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                m.prim(IDN, kk, jj, ii) = m.cons(IDN, kk, jj, ii);
                m.prim(IV1, kk, jj, ii) = v1;
                m.prim(IV2, kk, jj, ii) = v2;
                m.prim(IV3, kk, jj, ii) = v3;
                m.prim(IPN, kk, jj, ii) = pres(m.cons(IDN, kk, jj, ii), m.cons(IEN, kk, jj, ii),
                                                 m.cons(IM1, kk, jj, ii), m.cons(IM2, kk, jj, ii), m.cons(IM3, kk, jj, ii), m.hydro_gamma);
            }
        }
    }
}

double energy_from_temperature_protection(double &dens, double &ene, double &m1, double &m2, double &m3, double &minTemp, double &gamma){
    double KE = 0.5 * (m1 * m1 + m2 * m2 + m3 * m3) / dens;
    double eint = ene - KE;
    if (std::isnan(ene)){
        return dens * minTemp / (gamma - 1.) + KE;
    }
    return std::max(ene, dens * minTemp / (gamma - 1.) + KE);
}




