#include "../eos/eos.hpp"
#include "../index_def.hpp"
#include "hll.hpp"
#include <cmath>


void hll( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma){

    // calculate pressure and velocities
    double pL = pres(valsL[IDN], valsL[IEN], valsL[IM1], valsL[IM2], valsL[IM3], gamma);
    double pR = pres(valsR[IDN], valsR[IEN], valsR[IM1], valsR[IM2], valsR[IM3], gamma);
    double uL = valsL[IMP] / valsL[IDN];
    double uR = valsR[IMP] / valsR[IDN];
    // step 1: wave speed estimates, Toro 10.5.1
    double aL = soundspeed(valsL[IDN], pL, gamma);
    double aR = soundspeed(valsR[IDN], pR, gamma);
    double sL = uL - aL;
    double sR = uR + aR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < 5; val_ind ++){
        double flux_L = valsL[val_ind] * uL;
        double flux_R = valsR[val_ind] * uR;
        if (val_ind == IMP){
            flux_L += pL;
            flux_R += pR;
        }
        if (val_ind == IEN){
            flux_L += pL * uL;
            flux_R += pR * uR;
        }
        if      (0 <= sL){ fluxs[val_ind] = flux_L; }
        else if (0 >= sR ){ fluxs[val_ind] = flux_R; }
        else{
            fluxs[val_ind] = (sR * flux_L - sL * flux_R + sL * sR * (valsR[val_ind] - valsL[val_ind])) / (sR - sL);
        }
    }
}


void hll_grav( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &PhiL,
          double &PhiR,
          double &gamma){

    // calculate pressure and velocities
    double pL = pres(valsL[IDN], valsL[IEN], valsL[IM1], valsL[IM2], valsL[IM3], gamma);
    double pR = pres(valsR[IDN], valsR[IEN], valsR[IM1], valsR[IM2], valsR[IM3], gamma);
    double uL = valsL[IMP] / valsL[IDN];
    double uR = valsR[IMP] / valsR[IDN];
    // step 1: wave speed estimates, Toro 10.5.1
    double aL = soundspeed(valsL[IDN], pL, gamma);
    double aR = soundspeed(valsR[IDN], pR, gamma);
    double sL = uL - aL;
    double sR = uR + aR;
    // step 2: calculate flux
    double p_a_grhoL = pL + valsL[IDN] * PhiL;
    double p_a_grhoR = pR + valsR[IDN] * PhiR;
    double p_m_grhoL = pL - valsL[IDN] * PhiL;
    double p_m_grhoR = pR - valsR[IDN] * PhiR;
    for (int val_ind = 0; val_ind < 5; val_ind ++){
        double flux_L = valsL[val_ind] * uL;
        double flux_R = valsR[val_ind] * uR;
        if (val_ind == IMP){
            flux_L += p_m_grhoL;
            flux_R += p_m_grhoR;
        }
        if (val_ind == IEN){
            flux_L += p_a_grhoL * uL;
            flux_R += p_a_grhoR * uR;
        }
        if      (0 <= sL){ fluxs[val_ind] = flux_L; }
        else if (0 >= sR ){ fluxs[val_ind] = flux_R; }
        else{
            fluxs[val_ind] = (sR * flux_L - sL * flux_R + sL * sR * (valsR[val_ind] - valsL[val_ind])) / (sR - sL);
        }
    }
}


